#!/usr/bin/env python3
"""Compare original image data in two motion-control-stack savedData trees.

The default discovery rules understand the files written by
apps/python_fire_server/handle_data.py: the reference volume NRRD and the
per-slice NRRDs for subsequent volumes.  Derived registration products are
not part of the primary comparison.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import itertools
import json
import math
import os
import re
import sys
import textwrap
import warnings
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import SimpleITK as sitk


SUPPORTED_SUFFIXES = {".nhdr", ".nrrd", ".mha", ".mhd", ".nii", ".gz"}
DERIVED_MARKERS = (
    "_upsampled",
    "registered",
    "registration",
    "resampled",
    "difference",
    "_mask",
)
SLICE_ID_RE = re.compile(r"(?:^|_)volume_(\d+)_slice_(\d+)(?:_|\.|$)", re.I)
VOLUME_ID_RE = re.compile(
    r"(?:^|_)volume_(\d+)(?:_\d{8}T\d{6})?(?:\.(?:nhdr|nrrd|mha|mhd|nii(?:\.gz)?))$",
    re.I,
)

AXIS_RELATION_TOLERANCE = 1e-4
EXACT_INDEX_TOLERANCE = 1e-9
EXACT_MATRIX_TOLERANCE = 1e-9
EXACT_PHYSICAL_TOLERANCE_MM = 1e-6
SPACING_RELATIVE_TOLERANCE = 1e-6
SPACING_ABSOLUTE_TOLERANCE_MM = 1e-9
DEFAULT_SPATIAL_TOLERANCE_MM = 0.001
DEFAULT_PIXEL_CORRELATION_THRESHOLD = 0.99
INSUFFICIENT_OVERLAP_PERCENTAGE = 5.0
NRRD_FALLBACK_DTYPES = {
    "char": "i1", "int8": "i1", "signed char": "i1",
    "uchar": "u1", "uint8": "u1", "unsigned char": "u1",
    "short": "i2", "int16": "i2", "ushort": "u2", "uint16": "u2",
    "int": "i4", "int32": "i4", "uint": "u4", "uint32": "u4",
    "longlong": "i8", "int64": "i8", "ulonglong": "u8", "uint64": "u8",
    "float": "f4", "double": "f8",
}


@dataclass(frozen=True)
class ImageRecord:
    identifier: str
    path: Path
    relative_path: str
    identity_source: str


@dataclass
class LoadedImage:
    array: np.ndarray
    geometry: dict[str, Any]
    image: sitk.Image | None


def _json_value(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, tuple):
        return list(value)
    return value


def image_identifier(path: Path) -> tuple[str | None, str]:
    """Return a stable FIRE acquisition key, independent of protocol prefix."""
    name = path.name
    match = SLICE_ID_RE.search(name)
    if match:
        return f"volume:{int(match.group(1)):04d}/slice:{int(match.group(2)):04d}", "filename volume/slice indices"
    match = VOLUME_ID_RE.search(name)
    if match and "_slice_" not in name.lower():
        return f"volume:{int(match.group(1)):04d}", "filename volume index"
    return None, "no reliable acquisition identifier"


def _is_supported(path: Path) -> bool:
    lower = path.name.lower()
    return path.suffix.lower() in SUPPORTED_SUFFIXES and (not lower.endswith(".gz") or lower.endswith(".nii.gz"))


def discover_images(
    directory: Path, image_pattern: str | None = None
) -> tuple[dict[str, list[ImageRecord]], list[dict[str, str]]]:
    records: dict[str, list[ImageRecord]] = defaultdict(list)
    excluded: list[dict[str, str]] = []
    for path in sorted(directory.rglob("*")):
        if not path.is_file() or not _is_supported(path):
            continue
        relative = path.relative_to(directory).as_posix()
        if image_pattern and not (
            fnmatch.fnmatch(relative, image_pattern) or fnmatch.fnmatch(path.name, image_pattern)
        ):
            excluded.append({"path": relative, "reason": "does not match --image-pattern"})
            continue
        lower = path.name.lower()
        marker = next((item for item in DERIVED_MARKERS if item in lower), None)
        if marker:
            excluded.append({"path": relative, "reason": f"derived image marker: {marker}"})
            continue
        identifier, source = image_identifier(path)
        if identifier is None:
            excluded.append({"path": relative, "reason": source})
            continue
        records[identifier].append(ImageRecord(identifier, path, relative, source))
    return dict(records), excluded


def match_images(
    images_a: dict[str, list[ImageRecord]], images_b: dict[str, list[ImageRecord]]
) -> dict[str, Any]:
    duplicate_a = {key: [item.relative_path for item in value] for key, value in images_a.items() if len(value) > 1}
    duplicate_b = {key: [item.relative_path for item in value] for key, value in images_b.items() if len(value) > 1}
    ambiguous = set(duplicate_a) | set(duplicate_b)
    keys_a, keys_b = set(images_a), set(images_b)
    matched_keys = sorted((keys_a & keys_b) - ambiguous)
    return {
        "pairs": [(key, images_a[key][0], images_b[key][0]) for key in matched_keys],
        "missing_from_a": sorted(keys_b - keys_a),
        "missing_from_b": sorted(keys_a - keys_b),
        "duplicates_a": duplicate_a,
        "duplicates_b": duplicate_b,
        "ambiguous_identifiers": sorted(ambiguous),
    }


def _orientation_from_directions(directions: Sequence[list[float] | None]) -> str | None:
    spatial = [value for value in directions if value is not None]
    if len(spatial) != 3 or any(len(value) != 3 for value in spatial):
        return None
    matrix = np.asarray(spatial, dtype=np.float64).T
    norms = np.linalg.norm(matrix, axis=0)
    if np.any(norms == 0):
        return None
    unit = matrix / norms
    positive = ("L", "P", "S")
    negative = ("R", "A", "I")
    labels = []
    used: set[int] = set()
    for column in range(3):
        axis = int(np.argmax(np.abs(unit[:, column])))
        if axis in used:
            return None
        used.add(axis)
        labels.append(positive[axis] if unit[axis, column] >= 0 else negative[axis])
    return "".join(labels)


def _physical_points(image: sitk.Image, continuous_bounds: bool = False) -> list[list[float]]:
    """Return all image corner points in SimpleITK x,y,z index convention."""
    size = np.asarray(image.GetSize(), dtype=float)
    low = np.full(image.GetDimension(), -0.5 if continuous_bounds else 0.0)
    high = size - (0.5 if continuous_bounds else 1.0)
    points = []
    for choices in itertools.product((0, 1), repeat=image.GetDimension()):
        index = tuple(np.where(np.asarray(choices, dtype=bool), high, low).tolist())
        points.append(list(image.TransformContinuousIndexToPhysicalPoint(index)))
    return points


def _metadata_dict(image: sitk.Image) -> dict[str, str]:
    metadata = {}
    for key in image.GetMetaDataKeys():
        try:
            metadata[key] = image.GetMetaData(key)
        except RuntimeError:
            metadata[key] = "<unreadable>"
    return metadata


def _sitk_geometry(image: sitk.Image, array_shape: Sequence[int]) -> dict[str, Any]:
    dimension = image.GetDimension()
    direction = np.asarray(image.GetDirection(), dtype=float).reshape(dimension, dimension)
    spacing = list(image.GetSpacing())
    vectors = [(direction[:, index] * spacing[index]).tolist() for index in range(dimension)]
    center_index = tuple(((np.asarray(image.GetSize(), dtype=float) - 1.0) / 2.0).tolist())
    center_corners = _physical_points(image)
    support_corners = _physical_points(image, continuous_bounds=True)
    center_corner_array = np.asarray(center_corners)
    support_corner_array = np.asarray(support_corners)
    if dimension == 3:
        orientation = sitk.DICOMOrientImageFilter_GetOrientationFromDirectionCosines(image.GetDirection())
    else:
        orientation = _orientation_from_directions(vectors)
    return {
        "dimension": dimension,
        "dimensions": list(image.GetSize()),
        "numpy_shape": list(array_shape),
        "index_order": "SimpleITK (x,y,z,...); NumPy array axes are reversed (...,z,y,x)",
        "pixel_type": image.GetPixelIDTypeAsString(),
        "spacing": spacing,
        "origin": list(image.GetOrigin()),
        "direction_matrix": direction.tolist(),
        "directions": vectors,
        "orientation": orientation,
        "coordinate_space": "LPS (SimpleITK normalized)",
        "center_continuous_index": list(center_index),
        "center_physical_point": list(image.TransformContinuousIndexToPhysicalPoint(center_index)),
        "center_corner_physical_points": center_corners,
        "support_corner_physical_points": support_corners,
        "physical_bounding_box_min": np.min(support_corner_array, axis=0).tolist(),
        "physical_bounding_box_max": np.max(support_corner_array, axis=0).tolist(),
        "center_bounding_box_min": np.min(center_corner_array, axis=0).tolist(),
        "center_bounding_box_max": np.max(center_corner_array, axis=0).tolist(),
        "physical_extent_along_image_axes_mm": (np.asarray(image.GetSize()) * np.asarray(spacing)).tolist(),
        "metadata": _metadata_dict(image),
    }


def _load_detached_nrrd_fallback(path: Path, sitk_error: RuntimeError) -> LoadedImage:
    """Retain direct-array support when ITK cannot represent a nonspatial axis."""
    fields: dict[str, str] = {}
    filenames: list[str] = []
    list_mode = False
    for raw_line in path.read_text(encoding="ascii").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or line.startswith("NRRD"):
            continue
        if list_mode:
            filenames.append(line)
            continue
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        fields[key.strip().lower()] = value.strip()
        if key.strip().lower() == "data file" and value.strip().upper().startswith("LIST"):
            list_mode = True
    if not list_mode and fields.get("data file"):
        filenames = [fields["data file"]]
    type_name = fields.get("type", "").lower()
    if fields.get("encoding", "raw").lower() != "raw" or type_name not in NRRD_FALLBACK_DTYPES:
        raise sitk_error
    dtype = np.dtype(NRRD_FALLBACK_DTYPES[type_name]).newbyteorder(
        "<" if fields.get("endian", "little").lower() == "little" else ">"
    )
    sizes = tuple(int(item) for item in fields["sizes"].split())
    arrays = [np.fromfile(path.parent / filename, dtype=dtype) for filename in filenames]
    flat = arrays[0] if len(arrays) == 1 else np.concatenate(arrays)
    if flat.size != math.prod(sizes):
        raise ValueError(f"NRRD expected {math.prod(sizes)} values but found {flat.size}: {path}")
    array = flat.reshape(tuple(reversed(sizes)))
    geometry = {
        "dimension": int(fields.get("dimension", len(sizes))),
        "dimensions": list(sizes),
        "numpy_shape": list(array.shape),
        "index_order": "NRRD fastest-to-slowest axes; NumPy axes reversed",
        "pixel_type": str(dtype),
        "spacing": None,
        "origin": None,
        "direction_matrix": None,
        "directions": None,
        "orientation": None,
        "coordinate_space": fields.get("space"),
        "metadata": fields,
        "simpleitk_read_error": str(sitk_error),
    }
    return LoadedImage(array, geometry, None)


def load_image(path: Path) -> LoadedImage:
    """Read an image with SimpleITK and copy its reversed-axis NumPy view."""
    try:
        image = sitk.ReadImage(str(path))
    except RuntimeError as error:
        if path.suffix.lower() == ".nhdr":
            return _load_detached_nrrd_fallback(path, error)
        raise
    array = sitk.GetArrayFromImage(image)
    return LoadedImage(array, _sitk_geometry(image, array.shape), image)


def compare_geometry(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    components = {}
    for key in ("dimension", "dimensions", "spacing", "origin", "directions", "orientation", "coordinate_space"):
        components[key] = a.get(key) == b.get(key)
    return {
        "equal": all(components.values()),
        "status": "IDENTICAL" if all(components.values()) else "DIFFERENT",
        "components_equal": components,
        "a": a,
        "b": b,
    }


def _pearson(x: np.ndarray, y: np.ndarray) -> float | None:
    if x.size < 2:
        return None
    x_centered = x - np.mean(x)
    y_centered = y - np.mean(y)
    denominator = float(np.sqrt(np.sum(x_centered * x_centered) * np.sum(y_centered * y_centered)))
    if denominator == 0 or not math.isfinite(denominator):
        return None
    return float(np.sum(x_centered * y_centered) / denominator)


def _nonfinite_counts(array: np.ndarray) -> dict[str, int]:
    if not np.issubdtype(array.dtype, np.inexact):
        return {"nan": 0, "positive_inf": 0, "negative_inf": 0}
    return {
        "nan": int(np.isnan(array).sum()),
        "positive_inf": int(np.isposinf(array).sum()),
        "negative_inf": int(np.isneginf(array).sum()),
    }


def _integer_valued_count(array: np.ndarray) -> int:
    if not (np.issubdtype(array.dtype, np.integer) or np.issubdtype(array.dtype, np.floating)):
        return 0
    finite = np.isfinite(array)
    return int(np.count_nonzero(finite & (array == np.rint(array))))


def compare_arrays(a: np.ndarray, b: np.ndarray) -> dict[str, Any]:
    same_shape = a.shape == b.shape
    same_dtype = a.dtype == b.dtype
    result: dict[str, Any] = {
        "shape_a": list(a.shape), "shape_b": list(b.shape),
        "dtype_a": a.dtype.str, "dtype_b": b.dtype.str,
        "shape_equal": same_shape, "dtype_equal": same_dtype,
        "voxel_count": int(a.size) if same_shape else None,
        "unequal_voxels": None, "unequal_percentage": None,
        "max_absolute_difference": None, "mean_absolute_error": None,
        "rmse": None, "mean_signed_difference": None,
        "pearson_correlation": None, "byte_equal": False,
        "exact_voxel_equal": False,
        "nonfinite_a": _nonfinite_counts(a), "nonfinite_b": _nonfinite_counts(b),
        "integer_valued_count_a": _integer_valued_count(a),
        "integer_valued_count_b": _integer_valued_count(b),
    }
    if not same_shape:
        return result

    if np.issubdtype(a.dtype, np.inexact) or np.issubdtype(b.dtype, np.inexact):
        nan_equal = np.isnan(a) & np.isnan(b)
    else:
        nan_equal = np.zeros(a.shape, dtype=bool)
    equal = np.equal(a, b) | nan_equal
    unequal = int(equal.size - np.count_nonzero(equal))
    result["unequal_voxels"] = unequal
    result["unequal_percentage"] = 100.0 * unequal / equal.size if equal.size else 0.0
    result["exact_voxel_equal"] = bool(same_dtype and unequal == 0)
    result["byte_equal"] = bool(same_dtype and a.tobytes(order="C") == b.tobytes(order="C"))

    finite = np.isfinite(a) & np.isfinite(b)
    if np.any(finite):
        x = np.asarray(a[finite], dtype=np.float64)
        y = np.asarray(b[finite], dtype=np.float64)
        difference = x - y
        absolute = np.abs(difference)
        result["max_absolute_difference"] = float(np.max(absolute))
        result["mean_absolute_error"] = float(np.mean(absolute))
        result["rmse"] = float(np.sqrt(np.mean(difference * difference)))
        result["mean_signed_difference"] = float(np.mean(difference))
        result["pearson_correlation"] = _pearson(x, y)
        result["finite_pairs_for_statistics"] = int(x.size)
        result["_finite_sum_a"] = float(np.sum(x))
        result["_finite_sum_b"] = float(np.sum(y))
        result["_finite_sum_squares_a"] = float(np.sum(x * x))
        result["_finite_sum_squares_b"] = float(np.sum(y * y))
        result["_finite_sum_products"] = float(np.sum(x * y))
    else:
        result["finite_pairs_for_statistics"] = 0
    return result


def _index_to_physical_matrix(image: sitk.Image) -> np.ndarray:
    direction = np.asarray(image.GetDirection(), dtype=np.float64).reshape(
        image.GetDimension(), image.GetDimension()
    )
    return direction @ np.diag(np.asarray(image.GetSpacing(), dtype=np.float64))


def _index_mapping(source: sitk.Image, target: sitk.Image) -> tuple[np.ndarray, np.ndarray]:
    """Map source continuous indices to target continuous indices."""
    target_matrix = _index_to_physical_matrix(target)
    source_matrix = _index_to_physical_matrix(source)
    matrix = np.linalg.solve(target_matrix, source_matrix)
    offset = np.linalg.solve(
        target_matrix,
        np.asarray(source.GetOrigin(), dtype=np.float64)
        - np.asarray(target.GetOrigin(), dtype=np.float64),
    )
    return matrix, offset


def _signed_axis_relation(
    source: sitk.Image, target: sitk.Image, spatial_tolerance_mm: float
) -> dict[str, Any]:
    """Detect a source-index to target-index signed permutation."""
    if source.GetDimension() != target.GetDimension():
        return {"detected": False, "reason": "different dimensions"}
    dimension = source.GetDimension()
    matrix, offset = _index_mapping(source, target)
    permutation: list[int] = []
    signs: list[int] = []
    used: set[int] = set()
    for source_axis in range(dimension):
        target_axis = int(np.argmax(np.abs(matrix[:, source_axis])))
        if target_axis in used:
            return {
                "detected": False,
                "reason": "index mapping is not a one-to-one axis permutation",
                "index_matrix_source_to_target": matrix.tolist(),
                "index_offset_source_to_target": offset.tolist(),
            }
        used.add(target_axis)
        permutation.append(target_axis)
        signs.append(1 if matrix[target_axis, source_axis] >= 0 else -1)
    signed_permutation = np.zeros_like(matrix)
    for source_axis, target_axis in enumerate(permutation):
        signed_permutation[target_axis, source_axis] = signs[source_axis]
    matrix_error = float(np.max(np.abs(matrix - signed_permutation)))
    detected = matrix_error <= AXIS_RELATION_TOLERANCE
    if not detected:
        return {
            "detected": False,
            "reason": "axis mapping includes rotation, shear, or unequal spacing",
            "matrix_max_error": matrix_error,
            "index_matrix_source_to_target": matrix.tolist(),
            "index_offset_source_to_target": offset.tolist(),
        }

    source_size = np.asarray(source.GetSize(), dtype=int)
    target_size = np.asarray(target.GetSize(), dtype=int)
    size_compatible = all(
        source_size[source_axis] == target_size[target_axis]
        for source_axis, target_axis in enumerate(permutation)
    )
    ideal_offset = np.zeros(dimension, dtype=float)
    for source_axis, target_axis in enumerate(permutation):
        if signs[source_axis] < 0:
            ideal_offset[target_axis] = target_size[target_axis] - 1
    offset_error = offset - ideal_offset
    spacing_differences = []
    spacing_compatible = True
    for source_axis, target_axis in enumerate(permutation):
        source_spacing = source.GetSpacing()[source_axis]
        target_spacing = target.GetSpacing()[target_axis]
        spacing_differences.append(float(source_spacing - target_spacing))
        spacing_compatible &= bool(np.isclose(
            source_spacing,
            target_spacing,
            rtol=SPACING_RELATIVE_TOLERANCE,
            atol=SPACING_ABSOLUTE_TOLERANCE_MM,
        ))

    sample_discrepancies = []
    for choices in itertools.product((0, 1), repeat=dimension):
        source_index = np.where(
            np.asarray(choices, dtype=bool), source_size - 1, np.zeros(dimension)
        ).astype(float)
        target_index = signed_permutation @ source_index + ideal_offset
        source_point = np.asarray(
            source.TransformContinuousIndexToPhysicalPoint(tuple(source_index.tolist()))
        )
        target_point = np.asarray(
            target.TransformContinuousIndexToPhysicalPoint(tuple(target_index.tolist()))
        )
        sample_discrepancies.append(float(np.linalg.norm(source_point - target_point)))
    center_source = (source_size - 1) / 2.0
    center_target = signed_permutation @ center_source + ideal_offset
    center_discrepancy = float(
        np.linalg.norm(
            np.asarray(source.TransformContinuousIndexToPhysicalPoint(tuple(center_source.tolist())))
            - np.asarray(target.TransformContinuousIndexToPhysicalPoint(tuple(center_target.tolist())))
        )
    )
    sample_discrepancies.append(center_discrepancy)
    max_discrepancy = max(sample_discrepancies)
    exact_discrete_grid = bool(
        size_compatible
        and spacing_compatible
        and matrix_error <= EXACT_MATRIX_TOLERANCE
        and np.max(np.abs(offset_error)) <= EXACT_INDEX_TOLERANCE
        and max_discrepancy <= EXACT_PHYSICAL_TOLERANCE_MM
    )
    tolerance_discrete_grid = bool(
        size_compatible
        and spacing_compatible
        and matrix_error <= AXIS_RELATION_TOLERANCE
        and max_discrepancy <= spatial_tolerance_mm
    )
    return {
        "detected": True,
        "source_axis_to_target_axis": permutation,
        "source_axis_signs": signs,
        "reflected_source_axes": [index for index, sign in enumerate(signs) if sign < 0],
        "index_matrix_source_to_target": matrix.tolist(),
        "signed_permutation_matrix": signed_permutation.tolist(),
        "index_offset_source_to_target": offset.tolist(),
        "ideal_integer_offset_source_to_target": ideal_offset.tolist(),
        "index_offset_residual": offset_error.tolist(),
        "matrix_max_error": matrix_error,
        "size_compatible": size_compatible,
        "spacing_compatible": spacing_compatible,
        "spacing_differences_source_minus_target_mm": spacing_differences,
        "exact_discrete_grid_equivalent": exact_discrete_grid,
        "within_spatial_tolerance": tolerance_discrete_grid,
        "lossless_reorientation_feasible": tolerance_discrete_grid,
        "spatial_tolerance_mm": spatial_tolerance_mm,
        "corresponding_center_discrepancy_mm": center_discrepancy,
        "corresponding_sample_mean_discrepancy_mm": float(np.mean(sample_discrepancies)),
        "corresponding_sample_max_discrepancy_mm": max_discrepancy,
    }


def _corner_set_distances(points_a: Sequence[Sequence[float]], points_b: Sequence[Sequence[float]]) -> dict[str, float]:
    a = np.asarray(points_a, dtype=np.float64)
    b = np.asarray(points_b, dtype=np.float64)
    distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    nearest = np.concatenate((np.min(distances, axis=1), np.min(distances, axis=0)))
    return {"mean_mm": float(np.mean(nearest)), "max_mm": float(np.max(nearest))}


def _physical_overlap(source: sitk.Image, target: sitk.Image, relation: dict[str, Any]) -> dict[str, Any]:
    """Compare physical support extents, exactly for signed-axis relations."""
    if source.GetDimension() != target.GetDimension():
        return {"method": "unavailable", "percentage_of_a": 0.0, "percentage_of_b": 0.0}
    dimension = source.GetDimension()
    if relation.get("detected"):
        matrix = np.asarray(relation["index_matrix_source_to_target"], dtype=float)
        offset = np.asarray(relation["index_offset_source_to_target"], dtype=float)
        source_size = np.asarray(source.GetSize(), dtype=float)
        mapped = []
        for choices in itertools.product((0, 1), repeat=dimension):
            index = np.where(np.asarray(choices, dtype=bool), source_size - 0.5, -0.5)
            mapped.append(matrix @ index + offset)
        mapped = np.asarray(mapped)
        source_low, source_high = np.min(mapped, axis=0), np.max(mapped, axis=0)
        target_low = np.full(dimension, -0.5)
        target_high = np.asarray(target.GetSize(), dtype=float) - 0.5
        intersection_lengths = np.maximum(
            0.0, np.minimum(source_high, target_high) - np.maximum(source_low, target_low)
        )
        target_lengths = target_high - target_low
        source_lengths = source_high - source_low
        intersection = float(np.prod(intersection_lengths))
        volume_a = float(np.prod(target_lengths))
        volume_b = float(np.prod(source_lengths))
        union = volume_a + volume_b - intersection
        return {
            "method": "signed-axis support intersection in target continuous-index space",
            "percentage_of_a": 100.0 * intersection / volume_a if volume_a else 0.0,
            "percentage_of_b": 100.0 * intersection / volume_b if volume_b else 0.0,
            "intersection_over_union_percentage": 100.0 * intersection / union if union else 0.0,
            "source_support_bounds_in_target_index": [source_low.tolist(), source_high.tolist()],
            "target_support_bounds_in_target_index": [target_low.tolist(), target_high.tolist()],
        }

    # General fallback: conservative axis-aligned physical bounding-box overlap.
    corners_a = np.asarray(_physical_points(target, continuous_bounds=True))
    corners_b = np.asarray(_physical_points(source, continuous_bounds=True))
    low_a, high_a = np.min(corners_a, axis=0), np.max(corners_a, axis=0)
    low_b, high_b = np.min(corners_b, axis=0), np.max(corners_b, axis=0)
    lengths = np.maximum(0.0, np.minimum(high_a, high_b) - np.maximum(low_a, low_b))
    intersection = float(np.prod(lengths))
    volume_a = float(np.prod(high_a - low_a))
    volume_b = float(np.prod(high_b - low_b))
    return {
        "method": "axis-aligned physical bounding-box approximation",
        "percentage_of_a": 100.0 * intersection / volume_a if volume_a else 0.0,
        "percentage_of_b": 100.0 * intersection / volume_b if volume_b else 0.0,
        "intersection_over_union_percentage": (
            100.0 * intersection / (volume_a + volume_b - intersection)
            if volume_a + volume_b > intersection else 0.0
        ),
    }


def _apply_lossless_reorientation(array: np.ndarray, relation: dict[str, Any]) -> np.ndarray:
    """Reindex a NumPy array using a physically verified SimpleITK index map."""
    dimension = array.ndim
    source_index_order = np.transpose(array, axes=tuple(reversed(range(dimension))))
    permutation = relation["source_axis_to_target_axis"]
    source_axis_for_target = [permutation.index(target_axis) for target_axis in range(dimension)]
    aligned = np.transpose(source_index_order, axes=source_axis_for_target)
    signs = relation["source_axis_signs"]
    for target_axis, source_axis in enumerate(source_axis_for_target):
        if signs[source_axis] < 0:
            aligned = np.flip(aligned, axis=target_axis)
    return np.transpose(aligned, axes=tuple(reversed(range(dimension))))


def _target_sampling_overlap_mask(source: sitk.Image, target: sitk.Image) -> np.ndarray:
    """Mask target centers lying inside source voxel support, in NumPy order."""
    matrix, offset = _index_mapping(target, source)
    target_shape_numpy = tuple(reversed(target.GetSize()))
    numpy_indices = np.indices(target_shape_numpy, dtype=np.float64)
    target_indices = [numpy_indices[target.GetDimension() - 1 - axis] for axis in range(target.GetDimension())]
    valid = np.ones(target_shape_numpy, dtype=bool)
    source_size = np.asarray(source.GetSize(), dtype=float)
    for source_axis in range(source.GetDimension()):
        coordinate = np.full(target_shape_numpy, offset[source_axis], dtype=np.float64)
        for target_axis in range(target.GetDimension()):
            coordinate += matrix[source_axis, target_axis] * target_indices[target_axis]
        valid &= coordinate >= -0.5
        valid &= coordinate < source_size[source_axis] - 0.5
    return valid


def _intensity_distribution(values: np.ndarray) -> dict[str, Any]:
    values = np.asarray(values, dtype=np.float64)
    if not values.size:
        return {"count": 0}
    return {
        "count": int(values.size),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "mean": float(np.mean(values)),
        "standard_deviation": float(np.std(values)),
        "median": float(np.median(values)),
        "integer_valued_count": int(np.count_nonzero(values == np.rint(values))),
        "integer_valued_percentage": float(100.0 * np.count_nonzero(values == np.rint(values)) / values.size),
    }


def _diagnostic_affine_intensity_fit(native: np.ndarray, converted: np.ndarray) -> dict[str, Any]:
    native = np.asarray(native, dtype=np.float64)
    converted = np.asarray(converted, dtype=np.float64)
    if native.size < 2 or np.var(converted) == 0:
        return {"available": False, "reason": "insufficient converted intensity variance"}
    design = np.column_stack((converted, np.ones(converted.size)))
    scale, offset = np.linalg.lstsq(design, native, rcond=None)[0]
    residual = native - (scale * converted + offset)
    return {
        "available": True,
        "model": "native = scale * converted + offset",
        "scale": float(scale),
        "offset": float(offset),
        "residual_mae": float(np.mean(np.abs(residual))),
        "residual_rmse": float(np.sqrt(np.mean(residual * residual))),
    }


def compare_physical_space(
    native: LoadedImage,
    converted: LoadedImage,
    keep_arrays: bool = False,
    spatial_tolerance_mm: float = DEFAULT_SPATIAL_TOLERANCE_MM,
    pixel_correlation_threshold: float = DEFAULT_PIXEL_CORRELATION_THRESHOLD,
    direct_metrics: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Independently assess spatial geometry and pixel correspondence."""
    if direct_metrics is None:
        direct_metrics = compare_arrays(native.array, converted.array)
    image_a, image_b = native.image, converted.image
    if image_a is None or image_b is None:
        pixel_status = "PIXEL_PASS_EXACT" if direct_metrics.get("exact_voxel_equal") else (
            "PIXEL_SIMILAR"
            if direct_metrics.get("pearson_correlation") is not None
            and direct_metrics["pearson_correlation"] >= pixel_correlation_threshold
            else "PIXEL_DIFFERENT"
        )
        return {
            "coordinate_convention": "unavailable",
            "alignment_method": "none",
            "interpolated": False,
            "lossless_reorientation_feasible": False,
            "physical_overlap": {"percentage_of_a": 0.0, "percentage_of_b": 0.0},
            "sampling_overlap_voxels": 0,
            "sampling_overlap_percentage_of_a": 0.0,
            "aligned_intensity_metrics": None,
            "diagnostic_affine_intensity_fit": {"available": False, "reason": "SimpleITK could not represent this image geometry"},
            "spatial_statuses": ["SPATIAL_INDETERMINATE"],
            "spatial_status": "SPATIAL_INDETERMINATE",
            "pixel_status": pixel_status,
            "pixel_comparison_basis": "direct_array",
            "statuses": ["SPATIAL_INDETERMINATE", pixel_status],
            "reason": "SimpleITK physical geometry unavailable; direct-array comparison retained",
        }
    relation = _signed_axis_relation(image_b, image_a, spatial_tolerance_mm)
    overlap = _physical_overlap(image_b, image_a, relation)
    center_delta = np.asarray(native.geometry["center_physical_point"]) - np.asarray(
        converted.geometry["center_physical_point"]
    )
    bbox_min_delta = np.asarray(native.geometry["physical_bounding_box_min"]) - np.asarray(
        converted.geometry["physical_bounding_box_min"]
    )
    bbox_max_delta = np.asarray(native.geometry["physical_bounding_box_max"]) - np.asarray(
        converted.geometry["physical_bounding_box_max"]
    )
    corner_centers = _corner_set_distances(
        native.geometry["center_corner_physical_points"],
        converted.geometry["center_corner_physical_points"],
    )
    corner_support = _corner_set_distances(
        native.geometry["support_corner_physical_points"],
        converted.geometry["support_corner_physical_points"],
    )
    result: dict[str, Any] = {
        "coordinate_convention": "SimpleITK physical LPS millimeters",
        "coordinate_frame_assumption": (
            "Identity physical transform is used because both images are interpreted by SimpleITK in LPS. "
            "The NRRD headers alone cannot independently prove that their origins reference the same scanner frame."
        ),
        "index_to_physical_formula": "p = origin + direction @ diag(spacing) @ index",
        "orientation_a": native.geometry["orientation"],
        "orientation_b": converted.geometry["orientation"],
        "spatial_tolerance_mm": spatial_tolerance_mm,
        "pixel_correlation_threshold": pixel_correlation_threshold,
        "axis_relation": relation,
        "lossless_reorientation_feasible": bool(
            relation.get("lossless_reorientation_feasible", False)
            and corner_support["max_mm"] <= spatial_tolerance_mm
        ),
        "physical_overlap": overlap,
        "center_difference_vector_mm": center_delta.tolist(),
        "center_difference_mm": float(np.linalg.norm(center_delta)),
        "center_corner_set_discrepancy": corner_centers,
        "support_corner_set_discrepancy": corner_support,
        "bounding_box_min_difference_vector_mm": bbox_min_delta.tolist(),
        "bounding_box_max_difference_vector_mm": bbox_max_delta.tolist(),
    }
    if relation.get("detected") and relation.get("size_compatible"):
        ideal_reoriented = _apply_lossless_reorientation(converted.array, relation)
        result["axis_reorientation_diagnostic_metrics"] = compare_arrays(
            native.array.reshape(-1), ideal_reoriented.reshape(-1)
        )
        result["axis_reorientation_diagnostic_warning"] = (
            "Index-only diagnostic using the detected signed permutation; this is not a valid physical alignment "
            "unless lossless_reorientation_feasible is true."
        )
    else:
        result["axis_reorientation_diagnostic_metrics"] = None

    aligned_array = None
    overlap_mask = None
    if result["lossless_reorientation_feasible"]:
        aligned_array = _apply_lossless_reorientation(converted.array, relation)
        overlap_mask = np.ones(native.array.shape, dtype=bool)
        method = "lossless_discrete_axis_reorientation"
        interpolated = False
    elif image_a.GetDimension() == image_b.GetDimension() and overlap["percentage_of_a"] > 0:
        resampler = sitk.ResampleImageFilter()
        resampler.SetReferenceImage(image_a)
        resampler.SetTransform(sitk.Transform(image_a.GetDimension(), sitk.sitkIdentity))
        resampler.SetInterpolator(sitk.sitkLinear)
        resampler.SetDefaultPixelValue(float("nan"))
        aligned_image = resampler.Execute(sitk.Cast(image_b, sitk.sitkFloat64))
        aligned_array = sitk.GetArrayFromImage(aligned_image)
        overlap_mask = _target_sampling_overlap_mask(image_b, image_a) & np.isfinite(aligned_array)
        method = "identity_physical_resampling_linear"
        interpolated = True
    else:
        method = "none"
        interpolated = False

    result["alignment_method"] = method
    result["interpolated"] = interpolated
    exact_spatial = bool(
        relation.get("exact_discrete_grid_equivalent", False)
        and corner_support["max_mm"] <= EXACT_PHYSICAL_TOLERANCE_MM
    )
    tolerance_spatial = bool(
        relation.get("within_spatial_tolerance", False)
        and corner_support["max_mm"] <= spatial_tolerance_mm
    )
    spatial_statuses = []
    if exact_spatial:
        spatial_statuses.append("SPATIAL_PASS_EXACT")
        spatial_status = "SPATIAL_PASS_EXACT"
    elif tolerance_spatial:
        spatial_statuses.append("SPATIAL_PASS_TOLERANCE")
        spatial_status = "SPATIAL_PASS_TOLERANCE"
    else:
        grid_failure = not relation.get("within_spatial_tolerance", False)
        coverage_failure = (
            corner_support["max_mm"] > spatial_tolerance_mm
            or min(overlap.get("percentage_of_a", 0.0), overlap.get("percentage_of_b", 0.0))
            < INSUFFICIENT_OVERLAP_PERCENTAGE
        )
        if grid_failure:
            spatial_statuses.append("SPATIAL_FAIL_GRID")
        if coverage_failure:
            spatial_statuses.append("SPATIAL_FAIL_COVERAGE")
        if not spatial_statuses:
            spatial_statuses.append("SPATIAL_INDETERMINATE")
        spatial_status = spatial_statuses[0]

    if aligned_array is not None and overlap_mask is not None:
        finite = overlap_mask & np.isfinite(native.array) & np.isfinite(aligned_array)
        result["sampling_overlap_voxels"] = int(np.count_nonzero(finite))
        result["sampling_overlap_percentage_of_a"] = float(100.0 * np.count_nonzero(finite) / finite.size)
        if np.any(finite):
            native_values = np.asarray(native.array[finite], dtype=np.float64)
            converted_values = np.asarray(aligned_array[finite], dtype=np.float64)
            aligned_metrics = compare_arrays(native_values, converted_values)
            result["aligned_intensity_metrics"] = aligned_metrics
            result["intensity_distribution_a"] = _intensity_distribution(native_values)
            result["intensity_distribution_b"] = _intensity_distribution(converted_values)
            fit = _diagnostic_affine_intensity_fit(native_values, converted_values)
            result["diagnostic_affine_intensity_fit"] = fit
            primary_rmse = aligned_metrics.get("rmse")
            if fit.get("available") and (
                abs(fit["scale"] - 1.0) > 1e-3
                or abs(fit["offset"]) > 1e-3
                or (
                    primary_rmse is not None
                    and primary_rmse > 0
                    and fit["residual_rmse"] < 0.99 * primary_rmse
                )
            ):
                result["intensity_scaling_diagnostic"] = "DIFFERENT_INTENSITY_SCALING"
        else:
            result["aligned_intensity_metrics"] = None
            result["diagnostic_affine_intensity_fit"] = {"available": False, "reason": "no overlapping samples"}
    else:
        result["sampling_overlap_voxels"] = 0
        result["sampling_overlap_percentage_of_a"] = 0.0
        result["aligned_intensity_metrics"] = None
        result["diagnostic_affine_intensity_fit"] = {"available": False, "reason": "alignment unavailable"}
    aligned_metrics = result.get("aligned_intensity_metrics")
    exact_direct = bool(direct_metrics.get("exact_voxel_equal"))
    exact_lossless_aligned = bool(
        not interpolated and aligned_metrics and aligned_metrics.get("exact_voxel_equal")
    )
    if exact_direct:
        pixel_status = "PIXEL_PASS_EXACT"
        pixel_basis = "direct_array"
    elif exact_lossless_aligned:
        pixel_status = "PIXEL_PASS_EXACT"
        pixel_basis = "lossless_physical_reorientation"
    else:
        comparison_metrics = aligned_metrics or direct_metrics
        pixel_basis = method if aligned_metrics else "direct_array"
        correlation = comparison_metrics.get("pearson_correlation") if comparison_metrics else None
        if correlation is None:
            pixel_status = "PIXEL_INCOMPARABLE"
        elif correlation >= pixel_correlation_threshold:
            pixel_status = "PIXEL_SIMILAR"
        else:
            pixel_status = "PIXEL_DIFFERENT"
    result["spatial_statuses"] = spatial_statuses
    result["spatial_status"] = spatial_status
    result["spatial_exact_discrete_grid_equivalence"] = exact_spatial
    result["spatial_tolerance_equivalence"] = tolerance_spatial
    result["pixel_status"] = pixel_status
    result["pixel_comparison_basis"] = pixel_basis
    result["statuses"] = spatial_statuses + [pixel_status]
    if result.get("intensity_scaling_diagnostic"):
        result["statuses"].append(result["intensity_scaling_diagnostic"])
    if keep_arrays:
        result["_aligned_array"] = aligned_array
        result["_overlap_mask"] = overlap_mask
    return result


def compare_pair(
    identifier: str,
    record_a: ImageRecord,
    record_b: ImageRecord,
    spatial_tolerance_mm: float = DEFAULT_SPATIAL_TOLERANCE_MM,
    pixel_correlation_threshold: float = DEFAULT_PIXEL_CORRELATION_THRESHOLD,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "image_identifier": identifier,
        "file_a": str(record_a.path),
        "file_b": str(record_b.path),
    }
    try:
        loaded_a = load_image(record_a.path)
        loaded_b = load_image(record_b.path)
        row.update(compare_arrays(loaded_a.array, loaded_b.array))
        geometry = compare_geometry(loaded_a.geometry, loaded_b.geometry)
        row["geometry"] = geometry
        row["geometry_status"] = geometry["status"]
        physical = compare_physical_space(
            loaded_a,
            loaded_b,
            spatial_tolerance_mm=spatial_tolerance_mm,
            pixel_correlation_threshold=pixel_correlation_threshold,
            direct_metrics=row,
        )
        row["physical_comparison"] = physical
        row["physical_statuses"] = physical["statuses"]
        row["spatial_status"] = physical["spatial_status"]
        row["spatial_statuses"] = physical["spatial_statuses"]
        row["pixel_status"] = physical["pixel_status"]
        if not row["shape_equal"]:
            row["comparison_result"] = "FAIL_SHAPE"
        elif not row["dtype_equal"]:
            row["comparison_result"] = "FAIL_DTYPE"
        elif not row["exact_voxel_equal"]:
            row["comparison_result"] = "FAIL_VOXELS"
        elif not geometry["equal"]:
            row["comparison_result"] = "PASS_VOXELS_ONLY"
        else:
            row["comparison_result"] = "PASS_EXACT"
    except Exception as error:
        row.update({
            "shape_equal": False, "dtype_equal": False, "exact_voxel_equal": False,
            "byte_equal": False, "voxel_count": None, "unequal_voxels": None,
            "unequal_percentage": None, "max_absolute_difference": None,
            "mean_absolute_error": None, "rmse": None, "mean_signed_difference": None,
            "pearson_correlation": None, "geometry_status": "NOT_COMPARED",
            "comparison_result": "ERROR", "error": f"{type(error).__name__}: {error}",
        })
    return row


def _weighted_metric(rows: Sequence[dict[str, Any]], key: str) -> float | None:
    values = [(row.get(key), row.get("finite_pairs_for_statistics", 0)) for row in rows]
    usable = [(float(value), int(weight)) for value, weight in values if value is not None and weight]
    total = sum(weight for _, weight in usable)
    return sum(value * weight for value, weight in usable) / total if total else None


def aggregate_rows(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    compatible = [row for row in rows if row.get("voxel_count") is not None]
    total_voxels = sum(int(row["voxel_count"]) for row in compatible)
    unequal_voxels = sum(int(row["unequal_voxels"] or 0) for row in compatible)
    squared_sum = sum(
        (float(row["rmse"]) ** 2) * int(row.get("finite_pairs_for_statistics", 0))
        for row in compatible if row.get("rmse") is not None
    )
    finite_count = sum(int(row.get("finite_pairs_for_statistics", 0)) for row in compatible)
    max_values = [float(row["max_absolute_difference"]) for row in compatible if row.get("max_absolute_difference") is not None]
    unequal_percentages = [float(row["unequal_percentage"]) for row in compatible if row.get("unequal_percentage") is not None]
    sum_a = sum(float(row.get("_finite_sum_a", 0.0)) for row in compatible)
    sum_b = sum(float(row.get("_finite_sum_b", 0.0)) for row in compatible)
    sum_squares_a = sum(float(row.get("_finite_sum_squares_a", 0.0)) for row in compatible)
    sum_squares_b = sum(float(row.get("_finite_sum_squares_b", 0.0)) for row in compatible)
    sum_products = sum(float(row.get("_finite_sum_products", 0.0)) for row in compatible)
    pearson = None
    if finite_count >= 2:
        covariance = sum_products - (sum_a * sum_b / finite_count)
        variance_a = max(0.0, sum_squares_a - (sum_a * sum_a / finite_count))
        variance_b = max(0.0, sum_squares_b - (sum_b * sum_b / finite_count))
        denominator = math.sqrt(variance_a * variance_b)
        if denominator:
            pearson = covariance / denominator
    geometry_component_mismatches = {
        key: sum(not row.get("geometry", {}).get("components_equal", {}).get(key, False) for row in rows)
        for key in ("dimension", "dimensions", "spacing", "origin", "directions", "orientation", "coordinate_space")
    }
    nonfinite_a = {
        key: sum(int(row.get("nonfinite_a", {}).get(key, 0)) for row in compatible)
        for key in ("nan", "positive_inf", "negative_inf")
    }
    nonfinite_b = {
        key: sum(int(row.get("nonfinite_b", {}).get(key, 0)) for row in compatible)
        for key in ("nan", "positive_inf", "negative_inf")
    }
    integer_a = sum(int(row.get("integer_valued_count_a", 0)) for row in compatible)
    integer_b = sum(int(row.get("integer_valued_count_b", 0)) for row in compatible)
    return {
        "total_voxels": total_voxels,
        "unequal_voxels": unequal_voxels,
        "unequal_percentage": 100.0 * unequal_voxels / total_voxels if total_voxels else None,
        "maximum_absolute_difference": max(max_values) if max_values else None,
        "mean_absolute_error": _weighted_metric(compatible, "mean_absolute_error"),
        "rmse": math.sqrt(squared_sum / finite_count) if finite_count else None,
        "mean_signed_difference": _weighted_metric(compatible, "mean_signed_difference"),
        "pearson_correlation": pearson,
        "finite_pairs_for_statistics": finite_count,
        "nonfinite_a": nonfinite_a,
        "nonfinite_b": nonfinite_b,
        "integer_valued_percentage_a": 100.0 * integer_a / total_voxels if total_voxels else None,
        "integer_valued_percentage_b": 100.0 * integer_b / total_voxels if total_voxels else None,
        "all_exact_voxel_equal": bool(rows) and all(row.get("exact_voxel_equal", False) for row in rows),
        "all_byte_equal": bool(rows) and all(row.get("byte_equal", False) for row in rows),
        "all_geometry_equal": bool(rows) and all(row.get("geometry_status") == "IDENTICAL" for row in rows),
        "shape_mismatch_count": sum(not row.get("shape_equal", False) for row in rows),
        "dtype_mismatch_count": sum(not row.get("dtype_equal", False) for row in rows),
        "geometry_mismatch_count": sum(row.get("geometry_status") == "DIFFERENT" for row in rows),
        "geometry_component_mismatch_counts": geometry_component_mismatches,
        "pairs_with_voxel_differences": sum((row.get("unequal_voxels") or 0) > 0 for row in compatible),
        "minimum_pair_unequal_percentage": min(unequal_percentages) if unequal_percentages else None,
        "median_pair_unequal_percentage": float(np.median(unequal_percentages)) if unequal_percentages else None,
        "maximum_pair_unequal_percentage": max(unequal_percentages) if unequal_percentages else None,
        "comparison_error_count": sum(row.get("comparison_result") == "ERROR" for row in rows),
    }


def _aggregate_metric_dicts(metrics: Sequence[dict[str, Any]]) -> dict[str, Any]:
    usable = [item for item in metrics if item and item.get("voxel_count") is not None]
    count = sum(int(item.get("finite_pairs_for_statistics", 0)) for item in usable)
    unequal = sum(int(item.get("unequal_voxels") or 0) for item in usable)
    total_voxels = sum(int(item.get("voxel_count") or 0) for item in usable)
    sum_a = sum(float(item.get("_finite_sum_a", 0.0)) for item in usable)
    sum_b = sum(float(item.get("_finite_sum_b", 0.0)) for item in usable)
    sum_aa = sum(float(item.get("_finite_sum_squares_a", 0.0)) for item in usable)
    sum_bb = sum(float(item.get("_finite_sum_squares_b", 0.0)) for item in usable)
    sum_ab = sum(float(item.get("_finite_sum_products", 0.0)) for item in usable)
    variance_a = max(0.0, sum_aa - sum_a * sum_a / count) if count else 0.0
    variance_b = max(0.0, sum_bb - sum_b * sum_b / count) if count else 0.0
    denominator = math.sqrt(variance_a * variance_b)
    pearson = (sum_ab - sum_a * sum_b / count) / denominator if count and denominator else None
    squared_error = sum_aa + sum_bb - 2.0 * sum_ab
    max_differences = [item["max_absolute_difference"] for item in usable if item.get("max_absolute_difference") is not None]
    return {
        "voxel_count": total_voxels,
        "unequal_voxels": unequal,
        "unequal_percentage": 100.0 * unequal / total_voxels if total_voxels else None,
        "maximum_absolute_difference": max(max_differences) if max_differences else None,
        "mean_absolute_error": (
            sum(float(item["mean_absolute_error"]) * int(item.get("finite_pairs_for_statistics", 0)) for item in usable if item.get("mean_absolute_error") is not None) / count
            if count else None
        ),
        "rmse": math.sqrt(max(0.0, squared_error) / count) if count else None,
        "mean_signed_difference": (sum_a - sum_b) / count if count else None,
        "pearson_correlation": pearson,
        "finite_pairs_for_statistics": count,
        "_finite_sum_a": sum_a,
        "_finite_sum_b": sum_b,
        "_finite_sum_squares_a": sum_aa,
        "_finite_sum_squares_b": sum_bb,
        "_finite_sum_products": sum_ab,
    }


def aggregate_physical_rows(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    physical = [row["physical_comparison"] for row in rows if row.get("physical_comparison")]
    aligned_metrics = [item.get("aligned_intensity_metrics") for item in physical]
    aligned = _aggregate_metric_dicts([item for item in aligned_metrics if item])
    axis_diagnostic = _aggregate_metric_dicts([
        item["axis_reorientation_diagnostic_metrics"]
        for item in physical if item.get("axis_reorientation_diagnostic_metrics")
    ])
    status_counts: dict[str, int] = defaultdict(int)
    method_counts: dict[str, int] = defaultdict(int)
    for item in physical:
        method_counts[item["alignment_method"]] += 1
        for status in item["statuses"]:
            status_counts[status] += 1
    weights = [
        int(row.get("voxel_count") or math.prod(row.get("geometry", {}).get("a", {}).get("dimensions", [1])))
        for row in rows if row.get("physical_comparison")
    ]
    total_weight = sum(weights)
    overlap_values = [item.get("sampling_overlap_percentage_of_a", 0.0) for item in physical]
    support_overlap_a = [item.get("physical_overlap", {}).get("percentage_of_a", 0.0) for item in physical]
    support_overlap_b = [item.get("physical_overlap", {}).get("percentage_of_b", 0.0) for item in physical]
    support_overlap_iou = [item.get("physical_overlap", {}).get("intersection_over_union_percentage", 0.0) for item in physical]
    relation_means = [
        item.get("axis_relation", {}).get("corresponding_sample_mean_discrepancy_mm")
        for item in physical
    ]
    relation_maxima = [
        item.get("axis_relation", {}).get("corresponding_sample_max_discrepancy_mm")
        for item in physical
    ]
    center_differences = [item.get("center_difference_mm") for item in physical if item.get("center_difference_mm") is not None]

    fit = {"available": False}
    n = aligned.get("finite_pairs_for_statistics", 0)
    if n and aligned.get("_finite_sum_squares_b", 0.0) - aligned.get("_finite_sum_b", 0.0) ** 2 / n > 0:
        sum_a = aligned["_finite_sum_a"]
        sum_b = aligned["_finite_sum_b"]
        sum_aa = aligned["_finite_sum_squares_a"]
        sum_bb = aligned["_finite_sum_squares_b"]
        sum_ab = aligned["_finite_sum_products"]
        variance_b = sum_bb - sum_b * sum_b / n
        scale = (sum_ab - sum_a * sum_b / n) / variance_b
        offset = (sum_a - scale * sum_b) / n
        residual_sse = (
            sum_aa + scale * scale * sum_bb + n * offset * offset
            - 2.0 * scale * sum_ab - 2.0 * offset * sum_a + 2.0 * scale * offset * sum_b
        )
        pair_fits = [item.get("diagnostic_affine_intensity_fit", {}) for item in physical]
        pair_fit_weight = sum(
            weight for pair_fit, weight in zip(pair_fits, weights) if pair_fit.get("available")
        )
        pair_residual_mae = sum(
            float(pair_fit.get("residual_mae", 0.0)) * weight
            for pair_fit, weight in zip(pair_fits, weights)
            if pair_fit.get("available")
        ) / pair_fit_weight if pair_fit_weight else None
        fit = {
            "available": True,
            "model": "native = scale * converted + offset",
            "scale": float(scale),
            "offset": float(offset),
            "residual_rmse": float(math.sqrt(max(0.0, residual_sse) / n)),
                "mean_pairwise_fit_residual_mae": float(pair_residual_mae) if pair_residual_mae is not None else None,
        }
    integer_a = sum(item.get("intensity_distribution_a", {}).get("integer_valued_count", 0) for item in physical)
    integer_b = sum(item.get("intensity_distribution_b", {}).get("integer_valued_count", 0) for item in physical)
    spatial_exact_count = sum(item.get("spatial_status") == "SPATIAL_PASS_EXACT" for item in physical)
    spatial_tolerance_count = sum(item.get("spatial_status") == "SPATIAL_PASS_TOLERANCE" for item in physical)
    spatial_pass_count = spatial_exact_count + spatial_tolerance_count
    pixel_exact_count = sum(item.get("pixel_status") == "PIXEL_PASS_EXACT" for item in physical)
    pixel_similar_count = sum(item.get("pixel_status") == "PIXEL_SIMILAR" for item in physical)
    orientation_relationship_counts: dict[str, int] = defaultdict(int)
    for item in physical:
        relation = item.get("axis_relation", {})
        if relation.get("detected"):
            key = (
                f"permutation={relation.get('source_axis_to_target_axis')};"
                f"signs={relation.get('source_axis_signs')}"
            )
        else:
            key = "no signed-axis relation"
        orientation_relationship_counts[key] += 1
    return {
        "pair_count": len(physical),
        "status_counts": dict(sorted(status_counts.items())),
        "alignment_method_counts": dict(sorted(method_counts.items())),
        "orientation_relationship_counts": dict(sorted(orientation_relationship_counts.items())),
        "lossless_reorientation_feasible_count": sum(item.get("lossless_reorientation_feasible", False) for item in physical),
        "spatial_exact_count": spatial_exact_count,
        "spatial_tolerance_count": spatial_tolerance_count,
        "spatial_corresponding_count": spatial_pass_count,
        "spatial_corresponding_percentage": 100.0 * spatial_pass_count / len(physical) if physical else None,
        "spatial_fail_grid_count": sum("SPATIAL_FAIL_GRID" in item.get("spatial_statuses", []) for item in physical),
        "spatial_fail_coverage_count": sum("SPATIAL_FAIL_COVERAGE" in item.get("spatial_statuses", []) for item in physical),
        "spatial_indeterminate_count": sum(item.get("spatial_status") == "SPATIAL_INDETERMINATE" for item in physical),
        "pixel_exact_count": pixel_exact_count,
        "pixel_similar_count": pixel_similar_count,
        "pixel_different_count": sum(item.get("pixel_status") == "PIXEL_DIFFERENT" for item in physical),
        "pixel_incomparable_count": sum(item.get("pixel_status") == "PIXEL_INCOMPARABLE" for item in physical),
        "mean_corresponding_coordinate_discrepancy_mm": (
            sum(float(value) * weight for value, weight in zip(relation_means, weights) if value is not None) /
            sum(weight for value, weight in zip(relation_means, weights) if value is not None)
            if any(value is not None for value in relation_means) else None
        ),
        "maximum_corresponding_coordinate_discrepancy_mm": max(
            (float(value) for value in relation_maxima if value is not None), default=None
        ),
        "mean_center_difference_mm": float(np.mean(center_differences)) if center_differences else None,
        "maximum_center_difference_mm": max(center_differences, default=None),
        "sampling_overlap_percentage_of_a": (
            sum(value * weight for value, weight in zip(overlap_values, weights)) / total_weight
            if total_weight else None
        ),
        "physical_support_overlap_percentage_of_a": (
            sum(value * weight for value, weight in zip(support_overlap_a, weights)) / total_weight
            if total_weight else None
        ),
        "physical_support_overlap_percentage_of_b": (
            sum(value * weight for value, weight in zip(support_overlap_b, weights)) / total_weight
            if total_weight else None
        ),
        "physical_support_intersection_over_union_percentage": (
            sum(value * weight for value, weight in zip(support_overlap_iou, weights)) / total_weight
            if total_weight else None
        ),
        "axis_reorientation_diagnostic_metrics": axis_diagnostic,
        "aligned_intensity_metrics": aligned,
        "diagnostic_affine_intensity_fit": fit,
        "aligned_integer_valued_percentage_a": 100.0 * integer_a / n if n else None,
        "aligned_integer_valued_percentage_b": 100.0 * integer_b / n if n else None,
        "all_pairs_same_physical_grid": bool(physical) and all(
            item.get("spatial_status") in ("SPATIAL_PASS_EXACT", "SPATIAL_PASS_TOLERANCE")
            for item in physical
        ),
    }


def _csv_row(row: dict[str, Any]) -> dict[str, Any]:
    fields = (
        "image_identifier", "file_a", "file_b", "shape_a", "shape_b", "dtype_a", "dtype_b",
        "shape_equal", "dtype_equal", "voxel_count", "unequal_voxels", "unequal_percentage",
        "max_absolute_difference", "mean_absolute_error", "rmse", "mean_signed_difference",
        "pearson_correlation", "byte_equal", "exact_voxel_equal", "geometry_status",
        "integer_valued_count_a", "integer_valued_count_b",
        "geometry_components_equal", "geometry_a", "geometry_b",
        "spatial_status", "spatial_statuses", "pixel_status", "pixel_comparison_basis",
        "spatial_tolerance_mm", "spatial_exact_discrete_grid_equivalence",
        "spatial_tolerance_equivalence", "physical_statuses", "physical_alignment_method", "orientation_a", "orientation_b",
        "origin_a", "origin_b", "spacing_a", "spacing_b", "direction_matrix_a", "direction_matrix_b",
        "physical_bbox_a", "physical_bbox_b", "physical_overlap_percentage_a",
        "sampling_overlap_percentage_a", "lossless_reorientation_feasible",
        "axis_permutation_source_to_target", "reflected_source_axes",
        "mean_coordinate_discrepancy_mm", "max_coordinate_discrepancy_mm",
        "aligned_pearson_correlation", "aligned_mae", "aligned_rmse",
        "aligned_max_absolute_difference", "aligned_mean_signed_difference",
        "aligned_unequal_percentage", "intensity_fit_scale", "intensity_fit_offset",
        "intensity_fit_residual_mae", "intensity_fit_residual_rmse",
        "axis_reorientation_diagnostic_correlation", "axis_reorientation_diagnostic_mae",
        "axis_reorientation_diagnostic_rmse",
        "physical_comparison", "comparison_result", "error",
    )
    output = {}
    for key in fields:
        value = row.get(key)
        if key.startswith("geometry_") and key != "geometry_status":
            geometry_key = key.removeprefix("geometry_")
            value = row.get("geometry", {}).get(geometry_key)
        physical = row.get("physical_comparison", {})
        nested_values = {
            "physical_alignment_method": physical.get("alignment_method"),
            "pixel_comparison_basis": physical.get("pixel_comparison_basis"),
            "spatial_tolerance_mm": physical.get("spatial_tolerance_mm"),
            "spatial_exact_discrete_grid_equivalence": physical.get("spatial_exact_discrete_grid_equivalence"),
            "spatial_tolerance_equivalence": physical.get("spatial_tolerance_equivalence"),
            "orientation_a": physical.get("orientation_a"),
            "orientation_b": physical.get("orientation_b"),
            "origin_a": row.get("geometry", {}).get("a", {}).get("origin"),
            "origin_b": row.get("geometry", {}).get("b", {}).get("origin"),
            "spacing_a": row.get("geometry", {}).get("a", {}).get("spacing"),
            "spacing_b": row.get("geometry", {}).get("b", {}).get("spacing"),
            "direction_matrix_a": row.get("geometry", {}).get("a", {}).get("direction_matrix"),
            "direction_matrix_b": row.get("geometry", {}).get("b", {}).get("direction_matrix"),
            "physical_bbox_a": [
                row.get("geometry", {}).get("a", {}).get("physical_bounding_box_min"),
                row.get("geometry", {}).get("a", {}).get("physical_bounding_box_max"),
            ],
            "physical_bbox_b": [
                row.get("geometry", {}).get("b", {}).get("physical_bounding_box_min"),
                row.get("geometry", {}).get("b", {}).get("physical_bounding_box_max"),
            ],
            "physical_overlap_percentage_a": physical.get("physical_overlap", {}).get("percentage_of_a"),
            "sampling_overlap_percentage_a": physical.get("sampling_overlap_percentage_of_a"),
            "lossless_reorientation_feasible": physical.get("lossless_reorientation_feasible"),
            "axis_permutation_source_to_target": physical.get("axis_relation", {}).get("source_axis_to_target_axis"),
            "reflected_source_axes": physical.get("axis_relation", {}).get("reflected_source_axes"),
            "mean_coordinate_discrepancy_mm": physical.get("axis_relation", {}).get("corresponding_sample_mean_discrepancy_mm"),
            "max_coordinate_discrepancy_mm": physical.get("axis_relation", {}).get("corresponding_sample_max_discrepancy_mm"),
            "aligned_pearson_correlation": physical.get("aligned_intensity_metrics", {}).get("pearson_correlation") if physical.get("aligned_intensity_metrics") else None,
            "aligned_mae": physical.get("aligned_intensity_metrics", {}).get("mean_absolute_error") if physical.get("aligned_intensity_metrics") else None,
            "aligned_rmse": physical.get("aligned_intensity_metrics", {}).get("rmse") if physical.get("aligned_intensity_metrics") else None,
            "aligned_max_absolute_difference": physical.get("aligned_intensity_metrics", {}).get("max_absolute_difference") if physical.get("aligned_intensity_metrics") else None,
            "aligned_mean_signed_difference": physical.get("aligned_intensity_metrics", {}).get("mean_signed_difference") if physical.get("aligned_intensity_metrics") else None,
            "aligned_unequal_percentage": physical.get("aligned_intensity_metrics", {}).get("unequal_percentage") if physical.get("aligned_intensity_metrics") else None,
            "intensity_fit_scale": physical.get("diagnostic_affine_intensity_fit", {}).get("scale"),
            "intensity_fit_offset": physical.get("diagnostic_affine_intensity_fit", {}).get("offset"),
            "intensity_fit_residual_mae": physical.get("diagnostic_affine_intensity_fit", {}).get("residual_mae"),
            "intensity_fit_residual_rmse": physical.get("diagnostic_affine_intensity_fit", {}).get("residual_rmse"),
            "axis_reorientation_diagnostic_correlation": physical.get("axis_reorientation_diagnostic_metrics", {}).get("pearson_correlation") if physical.get("axis_reorientation_diagnostic_metrics") else None,
            "axis_reorientation_diagnostic_mae": physical.get("axis_reorientation_diagnostic_metrics", {}).get("mean_absolute_error") if physical.get("axis_reorientation_diagnostic_metrics") else None,
            "axis_reorientation_diagnostic_rmse": physical.get("axis_reorientation_diagnostic_metrics", {}).get("rmse") if physical.get("axis_reorientation_diagnostic_metrics") else None,
        }
        if key in nested_values:
            value = nested_values[key]
        output[key] = json.dumps(value) if isinstance(value, (list, dict)) else _json_value(value)
    return output


def write_reports(output_dir: Path, summary: dict[str, Any], rows: Sequence[dict[str, Any]]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2, default=_json_value, allow_nan=False) + "\n", encoding="utf-8"
    )
    csv_rows = [_csv_row(row) for row in rows]
    fieldnames = list(csv_rows[0]) if csv_rows else list(_csv_row({}))
    with (output_dir / "image_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)

    counts = summary["image_counts"]
    aggregate = summary["aggregate"]
    physical = summary["physical_space_aggregate"]
    aligned = physical["aligned_intensity_metrics"]
    axis_diagnostic = physical["axis_reorientation_diagnostic_metrics"]
    intensity_fit = physical["diagnostic_affine_intensity_fit"]
    coverage = summary["coverage"]
    issues = summary["result_categories"]
    configuration = summary["validation_configuration"]
    lines = [
        "# Motion-control image comparison", "",
        f"- **Selected validation result ({configuration['validation_target']}):** `{summary['selected_validation_result']}`",
        f"- **Coverage:** `{summary['coverage_status']}`",
        f"- **Spatial correspondence:** `{summary['spatial_overall_status']}`",
        f"- **Pixel correspondence:** `{summary['pixel_overall_status']}`",
        f"- **A:** {summary['label_a']} — `{summary['dir_a']}`",
        f"- **B:** {summary['label_b']} — `{summary['dir_b']}`", "",
        "## A. Image discovery and matching", "",
        f"- Primary images discovered: {counts['a']} in A; {counts['b']} in B",
        f"- Matched unambiguous pairs: {counts['matched_pairs']}",
        f"- Missing from A: {len(coverage['missing_from_a'])}",
        f"- Missing from B: {len(coverage['missing_from_b'])}",
        f"- Ambiguous identifiers: {len(coverage['ambiguous_identifiers'])}",
        f"- Excluded supported image files: {counts['excluded_a']} in A; {counts['excluded_b']} in B",
        "- Matching method: volume and slice indices parsed from FIRE-server filenames; no positional pairing.", "",
        "## B. Spatial correspondence", "",
        f"- Overall status: `{summary['spatial_overall_status']}`",
        f"- Acceptance tolerance: {_format_number(configuration['spatial_tolerance_mm'])} mm",
        f"- Corresponding pairs: {physical['spatial_corresponding_count']} of {physical['pair_count']} ({_format_number(physical['spatial_corresponding_percentage'])}%)",
        f"- Exact discrete grids: {physical['spatial_exact_count']}; tolerance-based grids: {physical['spatial_tolerance_count']}",
        f"- Grid failures: {physical['spatial_fail_grid_count']}; coverage failures: {physical['spatial_fail_coverage_count']}; indeterminate: {physical['spatial_indeterminate_count']}",
        f"- Mean / maximum corresponding-coordinate discrepancy: {_format_number(physical['mean_corresponding_coordinate_discrepancy_mm'])} / {_format_number(physical['maximum_corresponding_coordinate_discrepancy_mm'])} mm",
        f"- Physical support overlap A / B / IoU: {_format_number(physical['physical_support_overlap_percentage_of_a'])}% / {_format_number(physical['physical_support_overlap_percentage_of_b'])}% / {_format_number(physical['physical_support_intersection_over_union_percentage'])}%",
        f"- Orientation relationships: `{json.dumps(physical['orientation_relationship_counts'], sort_keys=True)}`",
        f"- Lossless reorientation feasible: {physical['lossless_reorientation_feasible_count']} of {physical['pair_count']} pairs",
        f"- Alignment methods used only for pixel diagnostics: `{json.dumps(physical['alignment_method_counts'], sort_keys=True)}`",
        "- Coordinate convention: SimpleITK LPS physical coordinates in millimeters.",
        "- Absolute-frame caveat: equivalence of the scanner reference origin cannot be independently proven from NRRD alone.", "",
        "## C. Pixel-value correspondence", "",
        f"- Overall status: `{summary['pixel_overall_status']}`",
        f"- Similarity threshold: Pearson correlation ≥ {_format_number(configuration['pixel_correlation_threshold'])}",
        f"- Pair statuses: exact={physical['pixel_exact_count']}, similar={physical['pixel_similar_count']}, different={physical['pixel_different_count']}, incomparable={physical['pixel_incomparable_count']}", "",
        "### Direct stored-array baseline", "",
        f"- Total voxels compared: {aggregate['total_voxels']:,}",
        f"- Unequal voxels: {aggregate['unequal_voxels']:,} ({_format_number(aggregate['unequal_percentage'])}%)",
        f"- Maximum absolute difference: {_format_number(aggregate['maximum_absolute_difference'])}",
        f"- Mean absolute error: {_format_number(aggregate['mean_absolute_error'])}",
        f"- RMSE: {_format_number(aggregate['rmse'])}",
        f"- Mean signed difference (A - B): {_format_number(aggregate['mean_signed_difference'])}",
        f"- Aggregate Pearson correlation: {_format_number(aggregate['pearson_correlation'])}",
        f"- Integer-valued samples A / B: {_format_number(aggregate['integer_valued_percentage_a'])}% / {_format_number(aggregate['integer_valued_percentage_b'])}%",
        f"- Pairs with voxel differences: {aggregate['pairs_with_voxel_differences']} of {counts['matched_pairs']}",
        f"- Per-pair unequal percentage (min/median/max): {_format_number(aggregate['minimum_pair_unequal_percentage'])}% / {_format_number(aggregate['median_pair_unequal_percentage'])}% / {_format_number(aggregate['maximum_pair_unequal_percentage'])}%",
        f"- All matched arrays exactly equal (including dtype): {aggregate['all_exact_voxel_equal']}",
        f"- All matched array byte representations equal: {aggregate['all_byte_equal']}", "",
        f"- Difference distribution: {summary['difference_characterization']}", "",
        "### Signed-axis diagnostic", "",
        f"- Pearson correlation: {_format_number(axis_diagnostic['pearson_correlation'])}",
        f"- MAE / RMSE: {_format_number(axis_diagnostic['mean_absolute_error'])} / {_format_number(axis_diagnostic['rmse'])}",
        "- For tolerance-accepted grids this is the lossless aligned comparison; for failed grids it remains a nonphysical diagnostic.", "",
        "### Physically aligned comparison", "",
        f"- Methods: `{json.dumps(physical['alignment_method_counts'], sort_keys=True)}`",
        f"- Sampling overlap on A grid: {_format_number(physical['sampling_overlap_percentage_of_a'])}%",
        f"- Overlap voxels: {aligned['voxel_count']:,}",
        f"- Pearson correlation: {_format_number(aligned['pearson_correlation'])}",
        f"- Unequal voxels: {aligned['unequal_voxels']:,} ({_format_number(aligned['unequal_percentage'])}%)",
        f"- Maximum absolute difference: {_format_number(aligned['maximum_absolute_difference'])}",
        f"- MAE / RMSE: {_format_number(aligned['mean_absolute_error'])} / {_format_number(aligned['rmse'])}",
        f"- Mean signed difference (A - aligned B): {_format_number(aligned['mean_signed_difference'])}",
        f"- Integer-valued samples A / aligned B: {_format_number(physical['aligned_integer_valued_percentage_a'])}% / {_format_number(physical['aligned_integer_valued_percentage_b'])}%",
        f"- Diagnostic affine fit: native = {_format_number(intensity_fit.get('scale'))} × converted + {_format_number(intensity_fit.get('offset'))}",
        f"- Affine-fit residual RMSE: {_format_number(intensity_fit.get('residual_rmse'))}",
        f"- Mean per-pair affine-fit residual MAE: {_format_number(intensity_fit.get('mean_pairwise_fit_residual_mae'))}", "",
        "## D. Interpretation", "",
        f"1. **Spatial:** {summary['spatial_conclusion']}",
        f"2. **Pixels:** {summary['pixel_conclusion']}",
        f"3. **Retrospective motion characterization:** {summary['retrospective_motion_interpretation']}", "",
        f"Legacy exact-array categories retained for compatibility: {', '.join(f'`{item}`' for item in issues)}.", "",
    ]
    if coverage["missing_from_a"] or coverage["missing_from_b"] or coverage["ambiguous_identifiers"]:
        lines.extend(["## Coverage issues", "", "```json", json.dumps(coverage, indent=2), "```", ""])
    discrepant = [row for row in rows if row.get("comparison_result") != "PASS_EXACT"][:20]
    if discrepant:
        lines.extend(["## First significant discrepancies", ""])
        for row in discrepant:
            lines.append(
                f"- `{row['image_identifier']}`: `{row['comparison_result']}`; "
                f"spatial=`{row.get('spatial_status')}`, pixel=`{row.get('pixel_status')}`, "
                f"unequal={row.get('unequal_voxels')}"
            )
        lines.append("")
    if summary.get("visualizations"):
        lines.extend(["## Diagnostic visualizations", ""])
        for path in summary["visualizations"]:
            lines.append(f"- `{path}`")
        lines.append("")
    (output_dir / "comparison_report.md").write_text("\n".join(lines), encoding="utf-8")


def _format_number(value: Any) -> str:
    if value is None:
        return "N/A"
    return f"{value:.12g}" if isinstance(value, float) else str(value)


def _display_slice(array: np.ndarray) -> np.ndarray:
    view = np.asarray(array)
    while view.ndim > 2:
        view = view[view.shape[0] // 2]
    return view


def _symmetric_normalized_difference(native: np.ndarray, converted: np.ndarray) -> np.ndarray:
    """Return 200 * (native - converted) / (|native| + |converted|) in percent."""
    native_float = np.asarray(native, dtype=np.float64)
    converted_float = np.asarray(converted, dtype=np.float64)
    scale = np.maximum(np.abs(native_float), np.abs(converted_float))
    scaled_native = np.zeros_like(native_float)
    scaled_converted = np.zeros_like(converted_float)
    nonzero = np.isfinite(native_float) & np.isfinite(converted_float) & (scale != 0.0)
    np.divide(native_float, scale, out=scaled_native, where=nonzero)
    np.divide(converted_float, scale, out=scaled_converted, where=nonzero)
    denominator = np.abs(scaled_native) + np.abs(scaled_converted)
    result = np.zeros_like(native_float)
    np.divide(
        200.0 * (scaled_native - scaled_converted),
        denominator,
        out=result,
        where=denominator != 0.0,
    )
    return result


def _format_intensity_fit_annotation(fit: dict[str, Any], correlation: float | None) -> str:
    correlation_text = "N/A" if correlation is None else f"{correlation:.6f}"
    if not fit.get("available"):
        return f"Least-squares fit unavailable\nPearson r = {correlation_text}"
    offset = float(fit["offset"])
    offset_operator = "+" if offset >= 0.0 else "−"
    return (
        f"Native = {fit['scale']:.6f} × Converted {offset_operator} {abs(offset):.6f}\n"
        f"Pearson r = {correlation_text}"
    )


def _finite_symmetric_limit(values: np.ndarray, fallback: float = 1.0) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if not finite.size:
        return fallback
    limit = float(np.max(np.abs(finite)))
    return limit if limit > 0.0 else fallback


def _shared_intensity_limits(*images: np.ndarray) -> tuple[float, float]:
    finite_parts = [np.asarray(image, dtype=np.float64)[np.isfinite(image)] for image in images]
    finite_parts = [part for part in finite_parts if part.size]
    if not finite_parts:
        return 0.0, 1.0
    minimum = min(float(np.min(part)) for part in finite_parts)
    maximum = max(float(np.max(part)) for part in finite_parts)
    if minimum == maximum:
        padding = max(abs(minimum) * 0.01, 1.0)
        return minimum - padding, maximum + padding
    return minimum, maximum


def _create_diagnostic_figure(
    plt: Any,
    row: dict[str, Any],
    physical: dict[str, Any],
    native_slice: np.ndarray,
    aligned_slice: np.ndarray,
    mask_slice: np.ndarray,
    scatter_converted: np.ndarray,
    scatter_native: np.ndarray,
) -> tuple[Any, np.ndarray]:
    """Create the diagnostic 2x3 figure without changing comparison results."""
    native_float = np.asarray(native_slice, dtype=np.float64)
    aligned_float = np.asarray(aligned_slice, dtype=np.float64)
    mask_slice = np.asarray(mask_slice, dtype=bool)
    absolute_difference = np.where(mask_slice, native_float - aligned_float, np.nan)
    normalized_difference = np.where(
        mask_slice,
        _symmetric_normalized_difference(native_float, aligned_float),
        np.nan,
    )
    fit = physical["diagnostic_affine_intensity_fit"]
    correlation = physical["aligned_intensity_metrics"].get("pearson_correlation")

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    shared_minimum, shared_maximum = _shared_intensity_limits(native_float, aligned_float)
    extent = (-0.5, native_float.shape[1] - 0.5, -0.5, native_float.shape[0] - 0.5)

    image_panels = (
        (axes[0, 0], native_float, "Native (A physical grid)"),
        (axes[0, 1], aligned_float, "Converted aligned (B)"),
    )
    for axis, image, title in image_panels:
        shown = axis.imshow(
            image,
            cmap="gray",
            origin="lower",
            extent=extent,
            vmin=shared_minimum,
            vmax=shared_maximum,
        )
        axis.set_title(title)
        axis.set_xlabel("SimpleITK x index")
        axis.set_ylabel("SimpleITK y index")
        fig.colorbar(shown, ax=axis, shrink=0.72, pad=0.02, label="Intensity")

    axes[0, 2].axis("off")
    axes[0, 2].set_title("Numerical summary")
    method_lines = textwrap.wrap(
        physical["alignment_method"].replace("_", " "),
        width=32,
        subsequent_indent="        ",
    )
    summary_lines = [
        f"Identifier: {row['image_identifier']}",
        f"Method: {method_lines[0]}",
        *method_lines[1:],
        f"Interpolated: {physical['interpolated']}",
        f"Overlap: {physical['sampling_overlap_percentage_of_a']:.3f}%",
        f"Pearson r: {_format_number(correlation)}",
        f"Fit slope: {_format_number(fit.get('scale'))}",
        f"Fit intercept: {_format_number(fit.get('offset'))}",
    ]
    axes[0, 2].text(
        0.0,
        0.98,
        "\n".join(summary_lines),
        va="top",
        family="monospace",
        fontsize=9,
    )

    absolute_limit = _finite_symmetric_limit(absolute_difference)
    absolute_plot = axes[1, 0].imshow(
        absolute_difference,
        cmap="coolwarm",
        origin="lower",
        extent=extent,
        vmin=-absolute_limit,
        vmax=absolute_limit,
    )
    axes[1, 0].set_title("Absolute intensity difference: A − B")
    axes[1, 0].set_xlabel("SimpleITK x index")
    axes[1, 0].set_ylabel("SimpleITK y index")
    fig.colorbar(
        absolute_plot,
        ax=axes[1, 0],
        shrink=0.72,
        pad=0.02,
        label="A − aligned B intensity",
    )

    normalized_plot = axes[1, 1].imshow(
        normalized_difference,
        cmap="coolwarm",
        origin="lower",
        extent=extent,
        vmin=-100.0,
        vmax=100.0,
    )
    axes[1, 1].set_title("Normalized difference (%)\n200(A−B)/(|A|+|B|), display saturated at ±100%")
    axes[1, 1].set_xlabel("SimpleITK x index")
    axes[1, 1].set_ylabel("SimpleITK y index")
    finite_normalized = normalized_difference[np.isfinite(normalized_difference)]
    below_limit = bool(finite_normalized.size and np.min(finite_normalized) < -100.0)
    above_limit = bool(finite_normalized.size and np.max(finite_normalized) > 100.0)
    colorbar_extend = "both" if below_limit and above_limit else (
        "min" if below_limit else "max" if above_limit else "neither"
    )
    fig.colorbar(
        normalized_plot,
        ax=axes[1, 1],
        shrink=0.72,
        pad=0.02,
        label="Symmetric normalized difference (%)",
        extend=colorbar_extend,
    )

    scatter_axis = axes[1, 2]
    scatter_axis.scatter(scatter_converted, scatter_native, s=2, alpha=0.15)
    if fit.get("available") and scatter_converted.size:
        endpoints = np.asarray((np.min(scatter_converted), np.max(scatter_converted)))
        scatter_axis.plot(endpoints, fit["scale"] * endpoints + fit["offset"], color="red")
    scatter_axis.set_title("Native vs aligned converted intensity")
    scatter_axis.set_xlabel("Aligned converted MRD intensity")
    scatter_axis.set_ylabel("Native MRD intensity")
    if scatter_native.size:
        native_minimum = float(np.min(scatter_native))
        native_maximum = float(np.max(scatter_native))
        native_span = native_maximum - native_minimum
        annotation_band = max(native_span * 0.28, max(abs(native_minimum), 1.0) * 0.02)
        scatter_axis.set_ylim(native_minimum, native_maximum + annotation_band)
    scatter_axis.text(
        0.02,
        0.98,
        _format_intensity_fit_annotation(fit, correlation),
        transform=scatter_axis.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        family="monospace",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "0.7", "alpha": 0.9},
    )
    return fig, axes


def generate_diagnostic_visualizations(
    output_dir: Path,
    pair_records: dict[str, tuple[ImageRecord, ImageRecord]],
    rows: Sequence[dict[str, Any]],
    spatial_tolerance_mm: float = DEFAULT_SPATIAL_TOLERANCE_MM,
    pixel_correlation_threshold: float = DEFAULT_PIXEL_CORRELATION_THRESHOLD,
) -> list[str]:
    """Generate at most two physical-alignment figures."""
    candidates = [
        row for row in rows
        if row.get("physical_comparison", {}).get("aligned_intensity_metrics")
        and row["image_identifier"] in pair_records
    ]
    if not candidates:
        return []
    representative = next((row for row in candidates if row["image_identifier"] == "volume:0000"), candidates[0])
    worst = max(
        candidates,
        key=lambda row: row["physical_comparison"]["aligned_intensity_metrics"].get("rmse") or -1.0,
    )
    selected = [representative]
    if worst["image_identifier"] != representative["image_identifier"]:
        selected.append(worst)

    cache_dir = Path("/tmp/motion-control-stack-matplotlib")
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    try:
        import matplotlib
        matplotlib.use("Agg")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Unable to import Axes3D.*", category=UserWarning)
            import matplotlib.pyplot as plt
    except ImportError:
        return []

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    generated = []
    for row in selected:
        record_a, record_b = pair_records[row["image_identifier"]]
        native = load_image(record_a.path)
        converted = load_image(record_b.path)
        physical = compare_physical_space(
            native,
            converted,
            keep_arrays=True,
            spatial_tolerance_mm=spatial_tolerance_mm,
            pixel_correlation_threshold=pixel_correlation_threshold,
        )
        aligned = physical.pop("_aligned_array", None)
        mask = physical.pop("_overlap_mask", None)
        if aligned is None or mask is None:
            continue
        native_slice = _display_slice(native.array)
        aligned_slice = _display_slice(aligned)
        mask_slice = _display_slice(mask)
        valid = mask & np.isfinite(native.array) & np.isfinite(aligned)
        x = np.asarray(aligned[valid], dtype=float)
        y = np.asarray(native.array[valid], dtype=float)
        if x.size > 10000:
            sample_indices = np.linspace(0, x.size - 1, 10000, dtype=int)
            x, y = x[sample_indices], y[sample_indices]
        fig, _ = _create_diagnostic_figure(
            plt,
            row,
            physical,
            native_slice,
            aligned_slice,
            mask_slice,
            x,
            y,
        )
        safe_identifier = row["image_identifier"].replace(":", "-").replace("/", "_")
        destination = figure_dir / f"physical_comparison_{safe_identifier}.png"
        fig.savefig(destination, dpi=140)
        plt.close(fig)
        generated.append(destination.relative_to(output_dir).as_posix())
    return generated


def compare_directories(
    dir_a: Path,
    dir_b: Path,
    output_dir: Path,
    label_a: str | None = None,
    label_b: str | None = None,
    image_pattern: str | None = None,
    verbose: bool = False,
    spatial_tolerance_mm: float = DEFAULT_SPATIAL_TOLERANCE_MM,
    pixel_correlation_threshold: float = DEFAULT_PIXEL_CORRELATION_THRESHOLD,
    validation_target: str = "both",
) -> tuple[dict[str, Any], list[dict[str, Any]], int]:
    dir_a, dir_b, output_dir = dir_a.resolve(), dir_b.resolve(), output_dir.resolve()
    for directory in (dir_a, dir_b):
        if not directory.is_dir():
            raise ValueError(f"Input directory does not exist: {directory}")
    for input_dir in (dir_a, dir_b):
        if output_dir == input_dir or input_dir in output_dir.parents:
            raise ValueError(f"Output directory must not be an input directory or inside one: {output_dir}")
    if spatial_tolerance_mm <= 0:
        raise ValueError("--spatial-tolerance-mm must be greater than zero")
    if not 0.0 <= pixel_correlation_threshold <= 1.0:
        raise ValueError("--pixel-correlation-threshold must be between 0 and 1")
    if validation_target not in ("spatial", "pixels", "both"):
        raise ValueError("validation_target must be spatial, pixels, or both")

    images_a, excluded_a = discover_images(dir_a, image_pattern)
    images_b, excluded_b = discover_images(dir_b, image_pattern)
    matching = match_images(images_a, images_b)
    rows = []
    for index, (identifier, record_a, record_b) in enumerate(matching["pairs"], start=1):
        if verbose and (index == 1 or index % 250 == 0):
            print(f"Comparing {index}/{len(matching['pairs'])}: {identifier}", file=sys.stderr)
        rows.append(compare_pair(
            identifier,
            record_a,
            record_b,
            spatial_tolerance_mm,
            pixel_correlation_threshold,
        ))
    aggregate = aggregate_rows(rows)
    physical_aggregate = aggregate_physical_rows(rows)
    coverage_failure = bool(
        matching["missing_from_a"] or matching["missing_from_b"]
        or matching["ambiguous_identifiers"] or not rows
    )
    voxel_failure = not aggregate["all_exact_voxel_equal"]
    geometry_failure = not aggregate["all_geometry_equal"]
    categories = []
    if not coverage_failure and not voxel_failure and not geometry_failure:
        categories.append("PASS_EXACT")
    if not coverage_failure and not voxel_failure and geometry_failure:
        categories.append("PASS_VOXELS_ONLY")
    if voxel_failure:
        categories.append("FAIL_VOXELS")
    if coverage_failure:
        categories.append("FAIL_COVERAGE")
    if geometry_failure:
        categories.append("FAIL_GEOMETRY")
    legacy_overall = "PASS_EXACT" if categories == ["PASS_EXACT"] else (
        "PASS_VOXELS_ONLY" if not coverage_failure and not voxel_failure else "FAIL"
    )
    pair_count = physical_aggregate["pair_count"]
    if pair_count and physical_aggregate["spatial_exact_count"] == pair_count:
        spatial_overall_status = "SPATIAL_PASS_EXACT"
    elif pair_count and physical_aggregate["spatial_corresponding_count"] == pair_count:
        spatial_overall_status = "SPATIAL_PASS_TOLERANCE"
    elif physical_aggregate["spatial_fail_coverage_count"]:
        spatial_overall_status = "SPATIAL_FAIL_COVERAGE"
    elif physical_aggregate["spatial_fail_grid_count"]:
        spatial_overall_status = "SPATIAL_FAIL_GRID"
    else:
        spatial_overall_status = "SPATIAL_INDETERMINATE"
    if pair_count and physical_aggregate["pixel_exact_count"] == pair_count:
        pixel_overall_status = "PIXEL_PASS_EXACT"
    elif (
        physical_aggregate["aligned_intensity_metrics"].get("pearson_correlation") is not None
        and physical_aggregate["aligned_intensity_metrics"]["pearson_correlation"]
        >= pixel_correlation_threshold
    ):
        pixel_overall_status = "PIXEL_SIMILAR"
    elif physical_aggregate["pixel_different_count"]:
        pixel_overall_status = "PIXEL_DIFFERENT"
    else:
        pixel_overall_status = "PIXEL_INCOMPARABLE"
    spatial_pass = spatial_overall_status in ("SPATIAL_PASS_EXACT", "SPATIAL_PASS_TOLERANCE")
    pixel_pass = pixel_overall_status in ("PIXEL_PASS_EXACT", "PIXEL_SIMILAR")
    selected_failure = coverage_failure or (
        validation_target in ("spatial", "both") and not spatial_pass
    ) or (
        validation_target in ("pixels", "both") and not pixel_pass
    )
    overall = "FAIL" if selected_failure else "PASS"
    physical_statuses = [spatial_overall_status, pixel_overall_status]
    different_pairs = aggregate["pairs_with_voxel_differences"]
    if different_pairs == 0:
        difference_characterization = "No voxel differences were detected."
    elif different_pairs == len(rows):
        difference_characterization = (
            "Differences are distributed across every matched image; per-pair unequal fractions range from "
            f"{_format_number(aggregate['minimum_pair_unequal_percentage'])}% to "
            f"{_format_number(aggregate['maximum_pair_unequal_percentage'])}%."
        )
    else:
        difference_characterization = (
            f"Differences occur in {different_pairs} of {len(rows)} matched images; inspect the CSV to determine localization."
        )
    if spatial_overall_status == "SPATIAL_PASS_EXACT":
        physical_conclusion = "All matched pairs have exact discrete physical-grid correspondence."
    elif spatial_overall_status == "SPATIAL_PASS_TOLERANCE":
        physical_conclusion = (
            f"All matched pairs correspond spatially within the configured {spatial_tolerance_mm:g} mm tolerance; "
            "actual residuals remain reported independently."
        )
    elif spatial_overall_status in ("SPATIAL_FAIL_GRID", "SPATIAL_FAIL_COVERAGE"):
        physical_conclusion = "Spatial correspondence fails because grid and/or physical-coverage residuals exceed the configured tolerance."
    else:
        physical_conclusion = "Spatial correspondence is indeterminate."
    if pixel_overall_status == "PIXEL_PASS_EXACT":
        pixel_conclusion = "Corresponding pixel values are exactly equal in a valid direct or lossless comparison."
    elif pixel_overall_status == "PIXEL_SIMILAR":
        pixel_conclusion = (
            f"Pixel values are nonidentical but similar at the configured correlation threshold "
            f"({pixel_correlation_threshold:g})."
        )
    elif pixel_overall_status == "PIXEL_DIFFERENT":
        pixel_conclusion = "Pixel values differ beyond the configured similarity criterion."
    else:
        pixel_conclusion = "Pixel correspondence could not be evaluated."
    motion_suitable = bool(not coverage_failure and spatial_pass)
    if motion_suitable:
        motion_conclusion = (
            "The metadata support comparative retrospective motion characterization in a consistent physical frame, "
            "but do not guarantee identical registration outputs or prove absolute scanner-space accuracy."
        )
    else:
        motion_conclusion = (
            "A consistent physical frame for comparative retrospective motion characterization is not established."
        )
    conclusion = f"{physical_conclusion} {pixel_conclusion} {motion_conclusion}"
    summary = {
        "comparison_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "dir_a": str(dir_a), "dir_b": str(dir_b),
        "label_a": label_a or dir_a.name, "label_b": label_b or dir_b.name,
        "image_pattern": image_pattern,
        "validation_configuration": {
            "validation_target": validation_target,
            "spatial_tolerance_mm": spatial_tolerance_mm,
            "pixel_correlation_threshold": pixel_correlation_threshold,
        },
        "matching_method": "FIRE filename volume/slice acquisition identifiers; never lexicographic position",
        "image_counts": {
            "a": sum(len(value) for value in images_a.values()),
            "b": sum(len(value) for value in images_b.values()),
            "unique_identifiers_a": len(images_a), "unique_identifiers_b": len(images_b),
            "matched_pairs": len(rows), "excluded_a": len(excluded_a), "excluded_b": len(excluded_b),
        },
        "coverage": {key: matching[key] for key in (
            "missing_from_a", "missing_from_b", "duplicates_a", "duplicates_b", "ambiguous_identifiers"
        )},
        "excluded_files": {"a": excluded_a, "b": excluded_b},
        "aggregate": aggregate,
        "physical_space_aggregate": physical_aggregate,
        "physical_result_statuses": physical_statuses,
        "coverage_status": "COVERAGE_FAIL" if coverage_failure else "COVERAGE_PASS",
        "spatial_overall_status": spatial_overall_status,
        "pixel_overall_status": pixel_overall_status,
        "selected_validation_result": overall,
        "difference_characterization": difference_characterization,
        "representative_geometry_differences": [
            {
                "image_identifier": row["image_identifier"],
                "components_equal": row["geometry"]["components_equal"],
                "a": row["geometry"]["a"],
                "b": row["geometry"]["b"],
            }
            for row in rows if row.get("geometry_status") == "DIFFERENT"
        ][:10],
        "result_categories": categories,
        "overall_result": legacy_overall,
        "validation": {
            "coverage_complete": not coverage_failure,
            "exact_voxel_equality": not voxel_failure,
            "spatial_geometry_equality": not geometry_failure,
            "physical_voxel_grid_agreement": physical_aggregate["all_pairs_same_physical_grid"],
            "spatial_correspondence_pass": spatial_pass,
            "pixel_correspondence_pass": pixel_pass,
            "retrospective_motion_comparison_suitable": motion_suitable,
        },
        "spatial_conclusion": physical_conclusion,
        "pixel_conclusion": pixel_conclusion,
        "retrospective_motion_interpretation": motion_conclusion,
        "conclusion": conclusion,
    }
    pair_records = {identifier: (record_a, record_b) for identifier, record_a, record_b in matching["pairs"]}
    summary["visualizations"] = generate_diagnostic_visualizations(
        output_dir,
        pair_records,
        rows,
        spatial_tolerance_mm,
        pixel_correlation_threshold,
    )
    write_reports(output_dir, summary, rows)
    return summary, rows, 1 if selected_failure else 0


def default_output_dir(dir_a: Path, dir_b: Path) -> Path:
    return Path.cwd() / f"comparison_results_{dir_a.name}_vs_{dir_b.name}"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir-a", required=True, type=Path)
    parser.add_argument("--dir-b", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--label-a")
    parser.add_argument("--label-b")
    parser.add_argument("--image-pattern", help="Glob applied to relative paths and basenames")
    parser.add_argument(
        "--spatial-tolerance-mm",
        type=float,
        default=DEFAULT_SPATIAL_TOLERANCE_MM,
        help=f"Maximum corresponding-coordinate residual for spatial acceptance (default: {DEFAULT_SPATIAL_TOLERANCE_MM})",
    )
    parser.add_argument(
        "--pixel-correlation-threshold",
        type=float,
        default=DEFAULT_PIXEL_CORRELATION_THRESHOLD,
        help=f"Pearson threshold for PIXEL_SIMILAR (default: {DEFAULT_PIXEL_CORRELATION_THRESHOLD})",
    )
    parser.add_argument(
        "--validation-target",
        choices=("spatial", "pixels", "both"),
        default="both",
        help="Assessment(s) that determine the process exit code; all diagnostics are always calculated",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output_dir = args.output_dir or default_output_dir(args.dir_a, args.dir_b)
    try:
        summary, _, exit_code = compare_directories(
            args.dir_a, args.dir_b, output_dir, args.label_a, args.label_b,
            args.image_pattern, args.verbose, args.spatial_tolerance_mm,
            args.pixel_correlation_threshold, args.validation_target,
        )
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(f"Selected validation result: {summary['selected_validation_result']}")
    print(f"Coverage: {summary['coverage_status']}")
    print(f"Spatial: {summary['spatial_overall_status']}")
    print(f"Pixels: {summary['pixel_overall_status']}")
    print(f"Matched pairs: {summary['image_counts']['matched_pairs']}")
    print(f"Unequal voxels: {summary['aggregate']['unequal_voxels']}")
    print(f"Reports: {output_dir.resolve()}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())

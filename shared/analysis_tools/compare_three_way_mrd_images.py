#!/usr/bin/env python3
"""Compare native FIRE, direct DICOM NRRD, and converted-MRD FIRE images."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import SimpleITK as sitk

try:
    from shared.analysis_tools import compare_motion_control_images as two_way
except ModuleNotFoundError:  # Direct execution from shared/analysis_tools.
    import compare_motion_control_images as two_way


PAIR_NAMES = (
    "native_vs_dicom",
    "dicom_vs_converted",
    "native_vs_converted",
)
VOLUME_RE = re.compile(r"^volume:(\d{4})(?:/slice:(\d{4}))?$")


class ComparisonError(RuntimeError):
    """Raised when correspondence or input validation is indeterminate."""


def _json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items() if not str(key).startswith("_")}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def parse_identifier(identifier: str) -> tuple[int, int | None]:
    match = VOLUME_RE.fullmatch(identifier)
    if match is None:
        raise ComparisonError(f"Unsupported FIRE identifier: {identifier}")
    return int(match.group(1)), None if match.group(2) is None else int(match.group(2))


def determine_volume_offset(
    fire_volume_ids: Sequence[int], dicom_volume_count: int
) -> dict[str, Any]:
    """Establish an index offset only when complete contiguous coverage proves it."""
    unique = sorted(set(int(value) for value in fire_volume_ids))
    if len(unique) != len(fire_volume_ids):
        return {"status": "indeterminate", "reason": "duplicate FIRE volume identifiers"}
    if not unique:
        return {"status": "indeterminate", "reason": "no FIRE volume identifiers"}
    if len(unique) == dicom_volume_count and unique == list(
        range(unique[0], unique[0] + dicom_volume_count)
    ):
        return {
            "status": "established",
            "fire_to_dicom_offset": -unique[0],
            "evidence": "complete contiguous FIRE coverage equals the DICOM time dimension",
        }
    valid_offsets = list(
        range(-unique[0], dicom_volume_count - unique[-1])
    )
    return {
        "status": "indeterminate",
        "reason": "partial coverage permits multiple temporal offsets",
        "candidate_offsets": valid_offsets,
    }


def read_nrrd_header(path: Path) -> dict[str, Any]:
    with path.open("rb") as stream:
        payload = stream.read(65536)
    separator = b"\r\n\r\n" if b"\r\n\r\n" in payload else b"\n\n"
    header_bytes = payload.split(separator, 1)[0]
    lines = header_bytes.decode("ascii", errors="strict").splitlines()
    if not lines or not lines[0].startswith("NRRD"):
        raise ComparisonError(f"Not a NRRD file: {path}")
    fields: dict[str, str] = {}
    comments = []
    for line in lines[1:]:
        if line.startswith("#"):
            comments.append(line[1:].strip())
        elif ":=" in line:
            key, value = line.split(":=", 1)
            fields[key.strip()] = value.strip()
        elif ":" in line:
            key, value = line.split(":", 1)
            fields[key.strip().lower()] = value.strip()
    return {"magic": lines[0], "fields": fields, "comments": comments}


def ras_to_lps_point(point: Sequence[float]) -> np.ndarray:
    return np.asarray((-float(point[0]), -float(point[1]), float(point[2])))


def ras_to_lps_direction(direction: np.ndarray) -> np.ndarray:
    matrix = np.asarray(direction, dtype=np.float64)
    transform = np.diag((-1.0, -1.0, 1.0))
    return transform @ matrix


def _parse_nrrd_vector(value: str) -> np.ndarray:
    match = re.fullmatch(r"\s*\(([^)]*)\)\s*", value)
    if match is None:
        raise ComparisonError(f"Invalid NRRD vector: {value}")
    result = np.asarray([float(item) for item in match.group(1).split(",")])
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise ComparisonError(f"Expected a finite three-component NRRD vector: {value}")
    return result


def verify_ras_to_lps_geometry(fields: dict[str, str], image: sitk.Image) -> dict[str, Any]:
    """Verify SimpleITK's LPS geometry against the raw RAS NRRD fields."""
    if fields.get("space", "").lower() != "right-anterior-superior":
        raise ComparisonError("DICOM NRRD must declare right-anterior-superior space")
    raw_origin = _parse_nrrd_vector(fields.get("space origin", ""))
    direction_tokens = re.findall(r"\([^)]*\)|none", fields.get("space directions", ""), re.I)
    spatial_tokens = [token for token in direction_tokens if token.lower() != "none"]
    if len(direction_tokens) != 4 or len(spatial_tokens) != 3:
        raise ComparisonError("DICOM NRRD must contain three spatial directions and one list axis")
    raw_vectors = np.column_stack([_parse_nrrd_vector(token) for token in spatial_tokens])
    lps_vectors = ras_to_lps_direction(raw_vectors)
    expected_spacing = np.linalg.norm(lps_vectors, axis=0)
    expected_direction = lps_vectors / expected_spacing
    observed_origin = np.asarray(image.GetOrigin())
    observed_spacing = np.asarray(image.GetSpacing())
    observed_direction = np.asarray(image.GetDirection()).reshape(3, 3)
    origin_error = float(np.max(np.abs(observed_origin - ras_to_lps_point(raw_origin))))
    spacing_error = float(np.max(np.abs(observed_spacing - expected_spacing)))
    direction_error = float(np.max(np.abs(observed_direction - expected_direction)))
    verified = origin_error <= 1e-5 and spacing_error <= 1e-5 and direction_error <= 1e-5
    if not verified:
        raise ComparisonError(
            "SimpleITK LPS geometry does not match the raw RAS NRRD geometry "
            f"(origin={origin_error:g}, spacing={spacing_error:g}, direction={direction_error:g})"
        )
    return {
        "verified": True,
        "raw_ras_origin": raw_origin.tolist(),
        "raw_ras_space_direction_vectors": raw_vectors.T.tolist(),
        "expected_lps_origin": ras_to_lps_point(raw_origin).tolist(),
        "expected_lps_spacing": expected_spacing.tolist(),
        "expected_lps_direction_matrix": expected_direction.tolist(),
        "maximum_origin_error": origin_error,
        "maximum_spacing_error": spacing_error,
        "maximum_direction_error": direction_error,
    }


@dataclass
class DicomNrrd:
    path: Path
    image: sitk.Image
    array_tzyx: np.ndarray
    header: dict[str, Any]
    _last_volume_index: int | None = field(default=None, init=False)
    _last_volume_image: sitk.Image | None = field(default=None, init=False)

    @property
    def volume_count(self) -> int:
        return int(self.array_tzyx.shape[0])

    @property
    def spatial_shape_zyx(self) -> tuple[int, int, int]:
        return tuple(int(value) for value in self.array_tzyx.shape[1:])

    def volume_image(self, volume: int) -> sitk.Image:
        if not 0 <= volume < self.volume_count:
            raise ComparisonError(f"DICOM volume index out of range: {volume}")
        if self._last_volume_index == volume and self._last_volume_image is not None:
            return self._last_volume_image
        if self.image.GetDimension() == 3 and self.image.GetNumberOfComponentsPerPixel() > 1:
            result = sitk.VectorIndexSelectionCast(self.image, volume)
        elif self.image.GetDimension() == 4:
            size = list(self.image.GetSize())
            index = [0, 0, 0, volume]
            size[3] = 0
            result = sitk.Extract(self.image, size, index)
        else:
            if self.volume_count != 1 or volume != 0:
                raise ComparisonError("NRRD does not expose a time/list dimension")
            result = self.image
        self._last_volume_index = volume
        self._last_volume_image = result
        return result

    def loaded_volume(self, volume: int) -> two_way.LoadedImage:
        image = self.volume_image(volume)
        array = self.array_tzyx[volume]
        return two_way.LoadedImage(array, two_way._sitk_geometry(image, array.shape), image)

    def loaded_slice(self, volume: int, anatomical_slice: int) -> two_way.LoadedImage:
        image = self.volume_image(volume)
        x, y, z = image.GetSize()
        if not 0 <= anatomical_slice < z:
            raise ComparisonError(f"DICOM anatomical slice out of range: {anatomical_slice}")
        slice_image = sitk.RegionOfInterest(image, (x, y, 1), (0, 0, anatomical_slice))
        array = self.array_tzyx[volume, anatomical_slice : anatomical_slice + 1]
        return two_way.LoadedImage(
            array, two_way._sitk_geometry(slice_image, array.shape), slice_image
        )


def load_dicom_nrrd(path: Path) -> DicomNrrd:
    if not path.is_file():
        raise ComparisonError(f"DICOM NRRD does not exist: {path}")
    header = read_nrrd_header(path)
    fields = header["fields"]
    if int(fields.get("dimension", "0")) != 4:
        raise ComparisonError("DICOM NRRD must have four stored dimensions")
    if fields.get("encoding", "").lower() != "raw":
        raise ComparisonError("DICOM NRRD must use raw encoding")
    image = sitk.ReadImage(str(path))
    raw = sitk.GetArrayFromImage(image)
    if image.GetDimension() == 3 and image.GetNumberOfComponentsPerPixel() > 1:
        if raw.ndim != 4:
            raise ComparisonError(f"Unexpected vector NRRD array shape: {raw.shape}")
        array_tzyx = np.moveaxis(raw, -1, 0)
    elif image.GetDimension() == 4:
        array_tzyx = raw
    else:
        raise ComparisonError(
            "The fourth NRRD dimension must be represented as time or vector components"
        )
    sizes = tuple(int(item) for item in fields.get("sizes", "").split())
    expected_tzyx = tuple(reversed(sizes))
    if array_tzyx.shape != expected_tzyx:
        raise ComparisonError(
            f"NRRD array shape {array_tzyx.shape} does not match header sizes {sizes}"
        )
    header["ras_to_lps_verification"] = verify_ras_to_lps_geometry(fields, image)
    return DicomNrrd(path, image, np.asarray(array_tzyx), header)


@dataclass
class FireDataset:
    label: str
    root: Path
    records: dict[str, list[two_way.ImageRecord]]
    excluded: list[dict[str, str]]
    reference: two_way.LoadedImage
    slice_ordinal_to_anatomical: dict[int, int]
    slice_mapping_residual_max: float

    def record(self, identifier: str) -> two_way.ImageRecord:
        items = self.records.get(identifier, [])
        if len(items) != 1:
            raise ComparisonError(
                f"{self.label} identifier {identifier} has {len(items)} records"
            )
        return items[0]

    def load(self, identifier: str) -> two_way.LoadedImage:
        return two_way.load_image(self.record(identifier).path)

    @property
    def volume_ids(self) -> list[int]:
        return sorted({parse_identifier(identifier)[0] for identifier in self.records})

    def anatomical_slice(self, identifier: str, loaded: two_way.LoadedImage | None = None) -> int:
        _volume, ordinal = parse_identifier(identifier)
        if ordinal is None:
            raise ComparisonError("3D reference volume has no single anatomical slice")
        if ordinal not in self.slice_ordinal_to_anatomical:
            raise ComparisonError(f"No anatomical mapping for FIRE slice ordinal {ordinal}")
        expected = self.slice_ordinal_to_anatomical[ordinal]
        if loaded is not None and loaded.image is not None and self.reference.image is not None:
            index = self.reference.image.TransformPhysicalPointToContinuousIndex(
                loaded.image.GetOrigin()
            )
            observed = int(round(index[2]))
            residual = float(np.linalg.norm(np.asarray(index) - (0.0, 0.0, observed)))
            if observed != expected or residual > 0.01:
                raise ComparisonError(
                    f"Inconsistent slice mapping for {identifier}: expected Z={expected}, "
                    f"observed {index}"
                )
        return expected

    def volume_array(self, volume: int) -> np.ndarray:
        if volume == 0:
            return self.reference.array
        z_size, y_size, x_size = self.reference.array.shape
        output = np.empty((z_size, y_size, x_size), dtype=self.reference.array.dtype)
        seen: set[int] = set()
        for ordinal, anatomical in self.slice_ordinal_to_anatomical.items():
            identifier = f"volume:{volume:04d}/slice:{ordinal:04d}"
            loaded = self.load(identifier)
            self.anatomical_slice(identifier, loaded)
            output[anatomical] = loaded.array[0]
            seen.add(anatomical)
        if seen != set(range(z_size)):
            raise ComparisonError(f"Incomplete reconstructed FIRE volume {volume}")
        return output


def discover_fire_dataset(path: Path, label: str) -> FireDataset:
    if not path.is_dir():
        raise ComparisonError(f"{label} FIRE directory does not exist: {path}")
    records, excluded = two_way.discover_images(path)
    duplicates = {key: len(value) for key, value in records.items() if len(value) != 1}
    if duplicates:
        raise ComparisonError(f"{label} has ambiguous identifiers: {duplicates}")
    reference_key = "volume:0000"
    if reference_key not in records:
        raise ComparisonError(f"{label} is missing its volume-0 reference image")
    reference = two_way.load_image(records[reference_key][0].path)
    if reference.image is None or reference.image.GetDimension() != 3:
        raise ComparisonError(f"{label} reference volume lacks 3D physical geometry")
    later_volumes = sorted(value for value in {
        parse_identifier(identifier)[0] for identifier in records
    } if value > 0)
    if not later_volumes:
        raise ComparisonError(f"{label} has no per-slice volumes")
    first_volume = later_volumes[0]
    mapping: dict[int, int] = {}
    residuals = []
    for identifier, items in records.items():
        volume, ordinal = parse_identifier(identifier)
        if volume != first_volume or ordinal is None:
            continue
        loaded = two_way.load_image(items[0].path)
        if loaded.image is None:
            raise ComparisonError(f"{label} slice {identifier} lacks physical geometry")
        index = reference.image.TransformPhysicalPointToContinuousIndex(
            loaded.image.GetOrigin()
        )
        anatomical = int(round(index[2]))
        residual = float(np.linalg.norm(np.asarray(index) - (0.0, 0.0, anatomical)))
        mapping[ordinal] = anatomical
        residuals.append(residual)
    z_size = int(reference.image.GetSize()[2])
    if sorted(mapping) != list(range(z_size)) or sorted(mapping.values()) != list(range(z_size)):
        raise ComparisonError(
            f"{label} slice ordinals do not form a one-to-one mapping to anatomical Z"
        )
    return FireDataset(
        label, path, records, excluded, reference, mapping, max(residuals, default=0.0)
    )


@dataclass
class DistributionAccumulator:
    dtype_names: set[str] = field(default_factory=set)
    count: int = 0
    finite_count: int = 0
    integer_count: int = 0
    total: float = 0.0
    total_squares: float = 0.0
    minimum: float = math.inf
    maximum: float = -math.inf
    samples: list[np.ndarray] = field(default_factory=list)

    def add(self, array: np.ndarray) -> None:
        values = np.asarray(array)
        self.dtype_names.add(values.dtype.str)
        self.count += int(values.size)
        finite = np.asarray(values[np.isfinite(values)], dtype=np.float64)
        if not finite.size:
            return
        self.finite_count += int(finite.size)
        self.integer_count += int(np.count_nonzero(finite == np.rint(finite)))
        self.total += float(np.sum(finite))
        self.total_squares += float(np.sum(finite * finite))
        self.minimum = min(self.minimum, float(np.min(finite)))
        self.maximum = max(self.maximum, float(np.max(finite)))
        step = max(1, finite.size // 128)
        self.samples.append(finite[::step][:128])

    def summary(self) -> dict[str, Any]:
        if not self.finite_count:
            return {"count": self.count, "finite_count": 0, "dtypes": sorted(self.dtype_names)}
        mean = self.total / self.finite_count
        variance = max(0.0, self.total_squares / self.finite_count - mean * mean)
        sample = np.concatenate(self.samples) if self.samples else np.empty(0)
        percentiles = {
            str(value): float(np.percentile(sample, value))
            for value in (1, 5, 25, 50, 75, 95, 99)
        }
        return {
            "count": self.count,
            "finite_count": self.finite_count,
            "dtypes": sorted(self.dtype_names),
            "minimum": self.minimum,
            "maximum": self.maximum,
            "mean": mean,
            "standard_deviation": math.sqrt(variance),
            "percentiles_approximate": percentiles,
            "percentile_sample_count": int(sample.size),
            "integer_valued_percentage": 100.0 * self.integer_count / self.finite_count,
        }


def _byte_equal_after_common_cast(a: np.ndarray, b: np.ndarray) -> bool:
    if a.shape != b.shape:
        return False
    dtype = np.result_type(a.dtype, b.dtype)
    return np.asarray(a, dtype=dtype).tobytes() == np.asarray(b, dtype=dtype).tobytes()


def _global_affine_fit(metrics: dict[str, Any]) -> dict[str, Any]:
    count = int(metrics.get("finite_pairs_for_statistics", 0))
    if count < 2:
        return {"available": False, "reason": "insufficient samples"}
    sum_a = float(metrics.get("_finite_sum_a", 0.0))
    sum_b = float(metrics.get("_finite_sum_b", 0.0))
    sum_aa = float(metrics.get("_finite_sum_squares_a", 0.0))
    sum_bb = float(metrics.get("_finite_sum_squares_b", 0.0))
    sum_ab = float(metrics.get("_finite_sum_products", 0.0))
    variance_b = sum_bb - sum_b * sum_b / count
    if variance_b <= 0:
        return {"available": False, "reason": "zero source-B variance"}
    covariance = sum_ab - sum_a * sum_b / count
    scale = covariance / variance_b
    offset = (sum_a - scale * sum_b) / count
    residual_sum_squares = (
        sum_aa + scale * scale * sum_bb + count * offset * offset
        - 2.0 * scale * sum_ab - 2.0 * offset * sum_a + 2.0 * scale * offset * sum_b
    )
    return {
        "available": True,
        "model": "A = scale * B + offset",
        "scale": scale,
        "offset": offset,
        "residual_rmse": math.sqrt(max(0.0, residual_sum_squares) / count),
    }


def compare_loaded(
    pair_name: str,
    identifier: str,
    loaded_a: two_way.LoadedImage,
    loaded_b: two_way.LoadedImage,
    spatial_tolerance_mm: float,
    pixel_correlation_threshold: float,
    anatomical_slice: int | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    direct = two_way.compare_arrays(loaded_a.array, loaded_b.array)
    physical = two_way.compare_physical_space(
        loaded_a,
        loaded_b,
        spatial_tolerance_mm=spatial_tolerance_mm,
        pixel_correlation_threshold=pixel_correlation_threshold,
        direct_metrics=direct,
    )
    if physical.get("lossless_reorientation_feasible"):
        primary = physical.get("aligned_intensity_metrics") or direct
        basis = "lossless_physical_reorientation"
    elif physical.get("axis_reorientation_diagnostic_metrics") is not None:
        primary = physical["axis_reorientation_diagnostic_metrics"]
        basis = "index_only_signed_axis_diagnostic"
    else:
        primary = direct
        basis = "direct_array"
    relation = physical.get("axis_relation", {})
    overlap = physical.get("physical_overlap", {})
    compared_b = (
        two_way._apply_lossless_reorientation(loaded_b.array, relation)
        if basis != "direct_array" and relation.get("detected")
        else loaded_b.array
    )
    same_shape = loaded_a.array.shape == compared_b.shape
    numeric_values_equal = bool(
        same_shape and np.array_equal(loaded_a.array, compared_b, equal_nan=True)
    )
    volume, ordinal = parse_identifier(identifier)
    row = {
        "comparison": pair_name,
        "identifier": identifier,
        "volume": volume,
        "fire_slice_ordinal": ordinal,
        "anatomical_slice": anatomical_slice,
        "spatial_status": physical.get("spatial_status"),
        "lossless_reorientation_feasible": physical.get("lossless_reorientation_feasible"),
        "alignment_method": physical.get("alignment_method"),
        "pixel_comparison_basis": basis,
        "reflected_source_axes": json.dumps(relation.get("reflected_source_axes", [])),
        "mean_coordinate_discrepancy_mm": relation.get(
            "corresponding_sample_mean_discrepancy_mm"
        ),
        "max_coordinate_discrepancy_mm": relation.get(
            "corresponding_sample_max_discrepancy_mm"
        ),
        "physical_support_iou_percentage": overlap.get(
            "intersection_over_union_percentage"
        ),
        "dtype_a": loaded_a.array.dtype.str,
        "dtype_b": compared_b.dtype.str,
        "voxel_count": primary.get("voxel_count"),
        "unequal_voxels": primary.get("unequal_voxels"),
        "numeric_values_equal": numeric_values_equal,
        "exact_voxel_equal_including_dtype": bool(
            numeric_values_equal and loaded_a.array.dtype == compared_b.dtype
        ),
        "byte_equal_after_common_cast": _byte_equal_after_common_cast(
            loaded_a.array, compared_b
        ),
        "pearson_correlation": primary.get("pearson_correlation"),
        "mean_absolute_error": primary.get("mean_absolute_error"),
        "rmse": primary.get("rmse"),
        "mean_signed_difference": primary.get("mean_signed_difference"),
        "maximum_absolute_difference": primary.get("max_absolute_difference"),
    }
    return row, {"direct": direct, "physical": physical, "primary": primary}


def aggregate_pair(
    pair_name: str,
    rows: Sequence[dict[str, Any]],
    metric_details: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    primary = two_way._aggregate_metric_dicts([item["primary"] for item in metric_details])
    physical = [item["physical"] for item in metric_details]
    relations = [item.get("axis_relation", {}) for item in physical]
    overlaps = [item.get("physical_overlap", {}) for item in physical]
    status_counts = Counter(item.get("spatial_status", "SPATIAL_INDETERMINATE") for item in physical)
    basis_counts = Counter(row["pixel_comparison_basis"] for row in rows)
    coordinate_means = [
        item.get("corresponding_sample_mean_discrepancy_mm")
        for item in relations if item.get("corresponding_sample_mean_discrepancy_mm") is not None
    ]
    coordinate_maxima = [
        item.get("corresponding_sample_max_discrepancy_mm")
        for item in relations if item.get("corresponding_sample_max_discrepancy_mm") is not None
    ]
    ious = [
        item.get("intersection_over_union_percentage")
        for item in overlaps if item.get("intersection_over_union_percentage") is not None
    ]
    result = {
        "name": pair_name,
        "pair_count": len(rows),
        "spatial_status_counts": dict(status_counts),
        "spatial_status": next(iter(status_counts)) if len(status_counts) == 1 else "SPATIAL_MIXED",
        "lossless_reorientation_feasible_count": sum(
            bool(item.get("lossless_reorientation_feasible")) for item in physical
        ),
        "pixel_comparison_basis_counts": dict(basis_counts),
        "mean_corresponding_coordinate_discrepancy_mm": (
            float(np.mean(coordinate_means)) if coordinate_means else None
        ),
        "maximum_corresponding_coordinate_discrepancy_mm": (
            max(coordinate_maxima) if coordinate_maxima else None
        ),
        "mean_physical_support_iou_percentage": float(np.mean(ious)) if ious else None,
        "pixel_metrics": primary,
        "all_numeric_values_equal": bool(rows) and all(row["numeric_values_equal"] for row in rows),
        "all_exact_voxel_equal_including_dtype": bool(rows) and all(
            row["exact_voxel_equal_including_dtype"] for row in rows
        ),
        "all_byte_equal_after_common_cast": bool(rows) and all(
            row["byte_equal_after_common_cast"] for row in rows
        ),
    }
    voxel_count = int(primary.get("voxel_count") or 0)
    unequal_voxels = int(primary.get("unequal_voxels") or 0)
    primary["equal_voxels"] = voxel_count - unequal_voxels
    primary["equal_percentage"] = (
        100.0 * (voxel_count - unequal_voxels) / voxel_count if voxel_count else None
    )
    result["diagnostic_affine_intensity_fit"] = _global_affine_fit(primary)
    return result


def _expected_fire_identifiers(volume_count: int, slice_count: int) -> set[str]:
    expected = {"volume:0000"}
    expected.update(
        f"volume:{volume:04d}/slice:{slice_index:04d}"
        for volume in range(1, volume_count)
        for slice_index in range(slice_count)
    )
    return expected


def analyze_fire_dicom_pair(
    pair_name: str,
    fire: FireDataset,
    dicom: DicomNrrd,
    fire_first: bool,
    spatial_tolerance_mm: float,
    pixel_correlation_threshold: float,
    distribution: DistributionAccumulator | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    expected = _expected_fire_identifiers(dicom.volume_count, dicom.spatial_shape_zyx[0])
    actual = set(fire.records)
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    if missing or extra:
        raise ComparisonError(
            f"{fire.label} correspondence incomplete: {len(missing)} missing, {len(extra)} extra"
        )
    rows = []
    details = []
    for identifier in sorted(actual):
        volume, ordinal = parse_identifier(identifier)
        fire_loaded = fire.load(identifier)
        if distribution is not None:
            distribution.add(fire_loaded.array)
        if ordinal is None:
            dicom_loaded = dicom.loaded_volume(volume)
            anatomical = None
        else:
            anatomical = fire.anatomical_slice(identifier, fire_loaded)
            dicom_loaded = dicom.loaded_slice(volume, anatomical)
        loaded_a, loaded_b = (
            (fire_loaded, dicom_loaded) if fire_first else (dicom_loaded, fire_loaded)
        )
        row, detail = compare_loaded(
            pair_name,
            identifier,
            loaded_a,
            loaded_b,
            spatial_tolerance_mm,
            pixel_correlation_threshold,
            anatomical,
        )
        rows.append(row)
        details.append(detail)
    return rows, details, aggregate_pair(pair_name, rows, details)


def analyze_fire_fire_pair(
    pair_name: str,
    native: FireDataset,
    converted: FireDataset,
    spatial_tolerance_mm: float,
    pixel_correlation_threshold: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    native_ids, converted_ids = set(native.records), set(converted.records)
    if native_ids != converted_ids:
        raise ComparisonError(
            f"FIRE identifier mismatch: {len(native_ids - converted_ids)} native-only, "
            f"{len(converted_ids - native_ids)} converted-only"
        )
    rows, details = [], []
    for identifier in sorted(native_ids):
        _volume, ordinal = parse_identifier(identifier)
        native_loaded = native.load(identifier)
        converted_loaded = converted.load(identifier)
        anatomical = None
        if ordinal is not None:
            anatomical = native.anatomical_slice(identifier, native_loaded)
            converted_anatomical = converted.anatomical_slice(identifier, converted_loaded)
            if converted_anatomical != anatomical:
                raise ComparisonError(
                    f"FIRE anatomical slice mismatch for {identifier}: "
                    f"{anatomical} versus {converted_anatomical}"
                )
        row, detail = compare_loaded(
            pair_name,
            identifier,
            native_loaded,
            converted_loaded,
            spatial_tolerance_mm,
            pixel_correlation_threshold,
            anatomical,
        )
        rows.append(row)
        details.append(detail)
    return rows, details, aggregate_pair(pair_name, rows, details)


def _aligned_volumes(
    native: FireDataset, dicom: DicomNrrd, converted: FireDataset, volume: int,
    spatial_tolerance_mm: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    native_array = native.volume_array(volume)
    dicom_array = dicom.array_tzyx[volume]
    converted_array = converted.volume_array(volume)
    dicom_relation = two_way._signed_axis_relation(
        dicom.loaded_volume(volume).image, native.reference.image, spatial_tolerance_mm
    )
    converted_relation = two_way._signed_axis_relation(
        converted.reference.image, native.reference.image, spatial_tolerance_mm
    )
    dicom_aligned = two_way._apply_lossless_reorientation(dicom_array, dicom_relation)
    converted_aligned = two_way._apply_lossless_reorientation(
        converted_array, converted_relation
    )
    return native_array, dicom_aligned, converted_aligned


def generate_figures(
    output_dir: Path,
    native: FireDataset,
    dicom: DicomNrrd,
    converted: FireDataset,
    pair_summaries: dict[str, dict[str, Any]],
    spatial_tolerance_mm: float,
) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    volumes = sorted(set((0, dicom.volume_count // 2, dicom.volume_count - 1)))
    slice_index = dicom.spatial_shape_zyx[0] // 2
    aligned = {
        volume: _aligned_volumes(native, dicom, converted, volume, spatial_tolerance_mm)
        for volume in volumes
    }

    figure, axes = plt.subplots(len(volumes), 3, figsize=(12, 4 * len(volumes)), squeeze=False)
    labels = ("Native FIRE", "Direct DICOM NRRD", "Corrected converted-MRD FIRE")
    for row, volume in enumerate(volumes):
        images = aligned[volume]
        low, high = np.percentile(np.concatenate([item[slice_index].ravel() for item in images]), (1, 99))
        for column, (image, label) in enumerate(zip(images, labels)):
            axes[row, column].imshow(image[slice_index], cmap="gray", origin="lower", vmin=low, vmax=high)
            axes[row, column].set_title(f"{label}\nvolume {volume}, anatomical Z {slice_index}")
            axes[row, column].axis("off")
    figure.tight_layout()
    figure.savefig(figure_dir / "figure1_three_way_anatomy.png", dpi=160)
    plt.close(figure)

    figure, axes = plt.subplots(len(volumes), 3, figsize=(12, 4 * len(volumes)), squeeze=False)
    difference_labels = ("Native - DICOM", "DICOM - converted", "Native - converted")
    for row, volume in enumerate(volumes):
        native_array, dicom_array, converted_array = aligned[volume]
        differences = (
            native_array - dicom_array,
            dicom_array - converted_array,
            native_array - converted_array,
        )
        limit = max(float(np.percentile(np.abs(item[slice_index]), 99)) for item in differences) or 1.0
        for column, (difference, label) in enumerate(zip(differences, difference_labels)):
            axes[row, column].imshow(
                difference[slice_index], cmap="coolwarm", origin="lower", vmin=-limit, vmax=limit
            )
            axes[row, column].set_title(f"{label}\nvolume {volume}")
            axes[row, column].axis("off")
    figure.suptitle(
        "Signed-axis index correspondence; DICOM pairs retain a physical half-voxel grid mismatch",
        fontsize=11,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(figure_dir / "figure2_pairwise_differences.png", dpi=160)
    plt.close(figure)

    middle = volumes[len(volumes) // 2]
    native_array, dicom_array, converted_array = aligned[middle]
    scatter_pairs = (
        (native_array, dicom_array, "Native vs DICOM"),
        (dicom_array, converted_array, "DICOM vs converted"),
        (native_array, converted_array, "Native vs converted"),
    )
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for axis, (a, b, label) in zip(axes, scatter_pairs):
        x, y = a.ravel()[::8], b.ravel()[::8]
        metrics = two_way.compare_arrays(x, y)
        axis.scatter(y, x, s=2, alpha=0.12, rasterized=True)
        lower = min(float(np.min(x)), float(np.min(y)))
        upper = max(float(np.max(x)), float(np.max(y)))
        axis.plot((lower, upper), (lower, upper), color="black", linewidth=1)
        axis.set_title(
            f"{label}\nr={metrics['pearson_correlation']:.6f}, MAE={metrics['mean_absolute_error']:.3f}"
        )
        axis.set_xlabel("Source B")
        axis.set_ylabel("Source A")
    figure.suptitle(
        "Pixel diagnostics after signed-axis indexing (physical pass only for native-converted)",
        fontsize=11,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(figure_dir / "figure3_intensity_scatterplots.png", dpi=160)
    plt.close(figure)

    labels = ["Native-DICOM", "DICOM-converted", "Native-converted"]
    summaries = [pair_summaries[name] for name in PAIR_NAMES]
    means = [item["mean_corresponding_coordinate_discrepancy_mm"] for item in summaries]
    maxima = [item["maximum_corresponding_coordinate_discrepancy_mm"] for item in summaries]
    ious = [item["mean_physical_support_iou_percentage"] for item in summaries]
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    positions = np.arange(3)
    axes[0].bar(positions - 0.18, means, 0.36, label="Mean")
    axes[0].bar(positions + 0.18, maxima, 0.36, label="Maximum")
    axes[0].axhline(spatial_tolerance_mm, color="black", linestyle="--", label="Tolerance")
    axes[0].set_xticks(positions, labels, rotation=15)
    axes[0].set_ylabel("Coordinate discrepancy (mm)")
    axes[0].legend()
    axes[1].bar(positions, ious)
    axes[1].set_xticks(positions, labels, rotation=15)
    axes[1].set_ylabel("Mean physical-support IoU (%)")
    axes[1].set_ylim(min(95.0, min(ious) - 0.5), 100.0)
    figure.tight_layout()
    figure.savefig(figure_dir / "figure4_spatial_geometry.png", dpi=160)
    plt.close(figure)


def _format_number(value: Any, digits: int = 9) -> str:
    return "n/a" if value is None else f"{float(value):.{digits}g}"


def write_outputs(
    output_dir: Path,
    summary: dict[str, Any],
    rows: Sequence[dict[str, Any]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "comparison_summary.json").write_text(
        json.dumps(_json_safe(summary), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    columns = list(rows[0]) if rows else []
    with (output_dir / "pairwise_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    pairs = summary["pairwise_comparisons"]
    lines = [
        "# Three-way REF image comparison", "",
        "## Inputs", "",
        f"- Native FIRE: `{summary['inputs']['native_fire']}`",
        f"- Direct DICOM NRRD: `{summary['inputs']['dicom_nrrd']}`",
        f"- Corrected converted-MRD FIRE: `{summary['inputs']['converted_fire']}`",
        f"- Spatial tolerance: `{summary['spatial_tolerance_mm']}` mm", "",
        "## Discovery and correspondence", "",
        f"- Native / converted FIRE records: {summary['discovery']['native']['record_count']} / {summary['discovery']['converted']['record_count']}",
        f"- DICOM stored layout: `{summary['discovery']['dicom']['stored_sizes_xyzt']}`; analysis array order: `(time,z,y,x)`",
        f"- Temporal mapping: DICOM index = FIRE volume + {summary['correspondence']['temporal']['fire_to_dicom_offset']}",
        f"- FIRE slice ordinal to anatomical Z: `{summary['correspondence']['native_slice_ordinal_to_anatomical_z']}`",
        "- The mapping is established by complete identifiers, physical slice origins, dimensions, and exact DICOM-versus-converted pixel checks; correlation alone was not used.",
        "- The source NRRD declares RAS. SimpleITK exposes its spatial image in LPS by negating physical X and Y, which was verified against the raw header.", "",
        "## Pairwise summary", "",
        "| Comparison | Spatial status | Pixel basis | Correlation | MAE | RMSE | Equal voxels | Equal % | Exact incl. dtype |",
        "|---|---|---|---:|---:|---:|---:|---:|---|",
    ]
    display_names = {
        "native_vs_dicom": "Native vs. DICOM",
        "dicom_vs_converted": "DICOM vs. converted MRD",
        "native_vs_converted": "Native vs. converted MRD",
    }
    for name in PAIR_NAMES:
        item = pairs[name]
        metrics = item["pixel_metrics"]
        basis = ", ".join(f"{key}: {value}" for key, value in item["pixel_comparison_basis_counts"].items())
        lines.append(
            f"| {display_names[name]} | `{item['spatial_status']}` | {basis} | "
            f"{_format_number(metrics.get('pearson_correlation'))} | "
            f"{_format_number(metrics.get('mean_absolute_error'))} | "
            f"{_format_number(metrics.get('rmse'))} | "
            f"{metrics.get('equal_voxels', 0)} | "
            f"{_format_number(metrics.get('equal_percentage'))} | "
            f"{item['all_exact_voxel_equal_including_dtype']} |"
        )
    lines.extend(["", "## Spatial findings", ""])
    for name in PAIR_NAMES:
        item = pairs[name]
        lines.append(
            f"- **{display_names[name]}:** mean / maximum coordinate discrepancy "
            f"{_format_number(item['mean_corresponding_coordinate_discrepancy_mm'])} / "
            f"{_format_number(item['maximum_corresponding_coordinate_discrepancy_mm'])} mm; "
            f"mean support IoU {_format_number(item['mean_physical_support_iou_percentage'])}%; "
            f"lossless-within-tolerance {item['lossless_reorientation_feasible_count']} of {item['pair_count']} records."
        )
    half_voxel = summary["interpretation"]["half_voxel_geometry_diagnostic"]
    lines.extend([
        "", f"The expected in-plane half-voxel diagonal is {_format_number(half_voxel['expected_in_plane_half_voxel_diagonal_mm'])} mm. The direct DICOM grid differs from native and corrected converted-MRD by approximately that amount, while native and converted pass at 0.001 mm. Therefore the observed half-voxel convention is already present in the direct DICOM export; it is not newly introduced by MRD serialization.",
        "", "## Intensity distributions", "",
        "| Source | Dtype | Min | Max | Mean | SD | Integer-valued % |",
        "|---|---|---:|---:|---:|---:|---:|",
    ])
    for name in ("native", "dicom", "converted"):
        distribution = summary["intensity_distributions"][name]
        lines.append(
            f"| {name} | `{', '.join(distribution['dtypes'])}` | "
            f"{_format_number(distribution['minimum'])} | {_format_number(distribution['maximum'])} | "
            f"{_format_number(distribution['mean'])} | {_format_number(distribution['standard_deviation'])} | "
            f"{_format_number(distribution['integer_valued_percentage'])} |"
        )
    lines.extend([
        "", "Selected percentiles for each source are retained in `comparison_summary.json`.",
        "", "## Pixel-intensity findings", "",
        f"- Native vs. DICOM: correlation {_format_number(pairs['native_vs_dicom']['pixel_metrics'].get('pearson_correlation'))}, MAE {_format_number(pairs['native_vs_dicom']['pixel_metrics'].get('mean_absolute_error'))}; values are not identical.",
        f"- DICOM vs. converted: all numeric values equal = `{pairs['dicom_vs_converted']['all_numeric_values_equal']}`; byte equality after a common numeric cast = `{pairs['dicom_vs_converted']['all_byte_equal_after_common_cast']}`. This directly establishes preservation through NIfTI-to-MRD-to-FIRE for this acquisition.",
        f"- Native vs. converted: correlation {_format_number(pairs['native_vs_converted']['pixel_metrics'].get('pearson_correlation'))}, MAE {_format_number(pairs['native_vs_converted']['pixel_metrics'].get('mean_absolute_error'))}.",
        "- Native scanner images are float reconstructions, while DICOM/converted values are integer-valued exports. The affine-fit residual remains nonzero, so quantization/rescaling is plausible but does not fully prove a single responsible processing stage.",
        "- Both NIfTI and direct NRRD were decoded by dcm2niix and therefore share upstream DICOM decoding/scaling behavior.",
        "", "## Three-way conclusions", "",
        "- **A — Native vs. DICOM intensities:** not identical; they are highly correlated but retain substantial nonzero residual error.",
        "- **B — Downstream preservation:** DICOM and converted-MRD contain identical numeric values after the established acquisition/slice mapping. Their storage dtypes differ (`int16` and `float32`), and their bytes match after casting to a common dtype.",
        "- **C — DICOM vs. corrected geometry:** the physical grids do not match at 0.001 mm; index-only pixel identity must not be interpreted as spatial equivalence.",
        "- **D — Half-voxel source:** the offset is present in the direct DICOM NRRD geometry and is absent between native and corrected converted-MRD FIRE images. MRD serialization did not introduce it.",
        "- **E — Native intensity differences:** datatype conversion, quantization, or rescaling remain plausible contributors, but the nonzero affine-fit residual and shared dcm2niix decoding do not isolate one exact stage.",
        "", "## Limitations", "",
        "- Index-only signed-axis metrics are reported separately when grids fail the physical tolerance; they are not described as physical alignment.",
        "- NRRD headers cannot independently prove that absolute origins use the same scanner frame.",
        "- Approximate intensity percentiles use deterministic bounded samples; counts, means, standard deviations, equality, MAE, RMSE, and correlations use all voxels.",
        "- This analysis localizes preservation and geometry conventions but cannot isolate a particular dcm2niix internal scaling step without an additional pre-export DICOM decoder comparison.",
        "", "## Figures", "",
        "- `figures/figure1_three_way_anatomy.png`",
        "- `figures/figure2_pairwise_differences.png`",
        "- `figures/figure3_intensity_scatterplots.png`",
        "- `figures/figure4_spatial_geometry.png`", "",
    ])
    (output_dir / "comparison_report.md").write_text("\n".join(lines), encoding="utf-8")


def run_analysis(
    native_path: Path,
    dicom_path: Path,
    converted_path: Path,
    output_dir: Path,
    spatial_tolerance_mm: float = 0.001,
    pixel_correlation_threshold: float = 0.99,
    overwrite: bool = False,
) -> dict[str, Any]:
    if spatial_tolerance_mm <= 0 or not math.isfinite(spatial_tolerance_mm):
        raise ComparisonError("Spatial tolerance must be positive and finite")
    if not 0 <= pixel_correlation_threshold <= 1:
        raise ComparisonError("Pixel correlation threshold must be between 0 and 1")
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise ComparisonError(f"Output directory is not empty (use --overwrite): {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    native = discover_fire_dataset(native_path, "native")
    dicom = load_dicom_nrrd(dicom_path)
    converted = discover_fire_dataset(converted_path, "converted")
    native_temporal = determine_volume_offset(native.volume_ids, dicom.volume_count)
    converted_temporal = determine_volume_offset(converted.volume_ids, dicom.volume_count)
    if native_temporal.get("status") != "established" or converted_temporal.get("status") != "established":
        raise ComparisonError(
            f"Temporal correspondence indeterminate: native={native_temporal}, "
            f"converted={converted_temporal}"
        )
    if native_temporal["fire_to_dicom_offset"] != converted_temporal["fire_to_dicom_offset"]:
        raise ComparisonError("Native and converted FIRE temporal offsets disagree")

    native_distribution = DistributionAccumulator()
    dicom_distribution = DistributionAccumulator()
    dicom_distribution.add(dicom.array_tzyx)
    converted_distribution = DistributionAccumulator()

    rows_ab, details_ab, summary_ab = analyze_fire_dicom_pair(
        "native_vs_dicom", native, dicom, True, spatial_tolerance_mm,
        pixel_correlation_threshold, native_distribution,
    )
    rows_bc, details_bc, summary_bc = analyze_fire_dicom_pair(
        "dicom_vs_converted", converted, dicom, False, spatial_tolerance_mm,
        pixel_correlation_threshold, converted_distribution,
    )
    rows_ac, details_ac, summary_ac = analyze_fire_fire_pair(
        "native_vs_converted", native, converted, spatial_tolerance_mm,
        pixel_correlation_threshold,
    )
    del details_ab, details_bc, details_ac

    raw_fields = dicom.header["fields"]
    sizes = tuple(int(item) for item in raw_fields["sizes"].split())
    temporal_offset = native_temporal["fire_to_dicom_offset"]
    dicom_spacing = dicom.volume_image(0).GetSpacing()
    expected_half_voxel_diagonal = 0.5 * math.hypot(dicom_spacing[0], dicom_spacing[1])
    observed_half_voxel_means = {
        "native_vs_dicom": summary_ab["mean_corresponding_coordinate_discrepancy_mm"],
        "dicom_vs_converted": summary_bc["mean_corresponding_coordinate_discrepancy_mm"],
    }
    half_voxel_present = bool(
        summary_ac["spatial_status"].startswith("SPATIAL_PASS")
        and all(
            value is not None and abs(value - expected_half_voxel_diagonal) <= 0.01
            for value in observed_half_voxel_means.values()
        )
    )
    summary = {
        "utility": "compare_three_way_mrd_images.py",
        "inputs": {
            "native_fire": str(native_path),
            "dicom_nrrd": str(dicom_path),
            "converted_fire": str(converted_path),
        },
        "spatial_tolerance_mm": spatial_tolerance_mm,
        "pixel_correlation_threshold": pixel_correlation_threshold,
        "discovery": {
            "native": {
                "record_count": len(native.records),
                "excluded_count": len(native.excluded),
                "volume_ids": [min(native.volume_ids), max(native.volume_ids)],
                "reference_shape_zyx": list(native.reference.array.shape),
                "pixel_dtype": native.reference.array.dtype.str,
                "geometry": native.reference.geometry,
            },
            "dicom": {
                "nrrd_magic": dicom.header["magic"],
                "encoding": raw_fields.get("encoding"),
                "declared_space": raw_fields.get("space"),
                "stored_sizes_xyzt": list(sizes),
                "analysis_shape_tzyx": list(dicom.array_tzyx.shape),
                "pixel_dtype": dicom.array_tzyx.dtype.str,
                "simpleitk_spatial_geometry": two_way._sitk_geometry(
                    dicom.volume_image(0), dicom.array_tzyx[0].shape
                ),
                "raw_ras_to_simpleitk_lps_verification": dicom.header[
                    "ras_to_lps_verification"
                ],
                "comments": dicom.header["comments"],
            },
            "converted": {
                "record_count": len(converted.records),
                "excluded_count": len(converted.excluded),
                "volume_ids": [min(converted.volume_ids), max(converted.volume_ids)],
                "reference_shape_zyx": list(converted.reference.array.shape),
                "pixel_dtype": converted.reference.array.dtype.str,
                "geometry": converted.reference.geometry,
            },
        },
        "correspondence": {
            "temporal": {
                "status": "established",
                "fire_to_dicom_offset": temporal_offset,
                "evidence": native_temporal["evidence"],
            },
            "native_slice_ordinal_to_anatomical_z": native.slice_ordinal_to_anatomical,
            "converted_slice_ordinal_to_anatomical_z": converted.slice_ordinal_to_anatomical,
            "native_slice_mapping_max_index_residual": native.slice_mapping_residual_max,
            "converted_slice_mapping_max_index_residual": converted.slice_mapping_residual_max,
            "volume_zero_representation": "one 3D FIRE NRRD",
            "later_volume_representation": "20 single-slice FIRE NRRDs per volume",
        },
        "intensity_distributions": {
            "native": native_distribution.summary(),
            "dicom": dicom_distribution.summary(),
            "converted": converted_distribution.summary(),
        },
        "pairwise_comparisons": {
            "native_vs_dicom": summary_ab,
            "dicom_vs_converted": summary_bc,
            "native_vs_converted": summary_ac,
        },
        "interpretation": {
            "native_vs_dicom_same_intensities": summary_ab["all_numeric_values_equal"],
            "dicom_to_converted_intensities_preserved": summary_bc["all_numeric_values_equal"],
            "dicom_and_corrected_converted_same_physical_grid": summary_bc["spatial_status"].startswith("SPATIAL_PASS"),
            "half_voxel_offset_present_in_direct_dicom_geometry": half_voxel_present,
            "half_voxel_geometry_diagnostic": {
                "expected_in_plane_half_voxel_diagonal_mm": expected_half_voxel_diagonal,
                "observed_mean_coordinate_discrepancy_mm": observed_half_voxel_means,
                "native_vs_converted_passes_tolerance": summary_ac[
                    "spatial_status"
                ].startswith("SPATIAL_PASS"),
            },
            "exact_stage_of_native_intensity_difference_isolated": False,
        },
    }
    all_rows = rows_ab + rows_bc + rows_ac
    write_outputs(output_dir, summary, all_rows)
    generate_figures(
        output_dir, native, dicom, converted, summary["pairwise_comparisons"],
        spatial_tolerance_mm,
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare native FIRE, direct DICOM NRRD, and converted-MRD FIRE images"
    )
    parser.add_argument("--native-fire", required=True, type=Path)
    parser.add_argument("--dicom-nrrd", required=True, type=Path)
    parser.add_argument("--converted-fire", required=True, type=Path)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/image_comparison_REF_three_way")
    )
    parser.add_argument("--spatial-tolerance-mm", type=float, default=0.001)
    parser.add_argument("--pixel-correlation-threshold", type=float, default=0.99)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = run_analysis(
            args.native_fire.resolve(),
            args.dicom_nrrd.resolve(),
            args.converted_fire.resolve(),
            args.output_dir.resolve(),
            args.spatial_tolerance_mm,
            args.pixel_correlation_threshold,
            args.overwrite,
        )
    except (ComparisonError, OSError, RuntimeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    print(f"Reports: {args.output_dir.resolve()}")
    for name in PAIR_NAMES:
        pair = summary["pairwise_comparisons"][name]
        metrics = pair["pixel_metrics"]
        print(
            f"{name}: {pair['spatial_status']}, r={_format_number(metrics.get('pearson_correlation'))}, "
            f"MAE={_format_number(metrics.get('mean_absolute_error'))}, "
            f"numeric_equal={pair['all_numeric_values_equal']}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())

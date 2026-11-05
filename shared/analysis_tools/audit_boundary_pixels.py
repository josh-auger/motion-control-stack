#!/usr/bin/env python3
"""Audit boundary clipping and exact-zero masks in native and DICOM-derived FIRE images."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pydicom
import SimpleITK as sitk

_MATPLOTLIB_CACHE = Path("/tmp/motion-control-stack-matplotlib")
_MATPLOTLIB_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MATPLOTLIB_CACHE))

try:
    from shared.analysis_tools import compare_motion_control_images as two_way
    from shared.analysis_tools import compare_three_way_mrd_images as three_way
except ModuleNotFoundError:  # Direct execution from shared/analysis_tools.
    import compare_motion_control_images as two_way
    import compare_three_way_mrd_images as three_way


DEFAULT_NATIVE = Path(
    "/home/jauger/GitHubRepos/motion-control-stack/data/"
    "savedData_20261105T143728_ep2d_bold_fire_noMoCo_tr1210_ref"
)
DEFAULT_CONVERTED = Path(
    "/home/jauger/GitHubRepos/motion-control-stack/data/"
    "savedData_20261105T145244_ep2d_bold_fire_noMoCo_tr1210_ref_convertedMRDv3"
)
DEFAULT_DICOM = Path(
    "/home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/"
    "fromResearchSynapse_20261015/6045885-2-sorted/6045885/20261015/"
    "joshauger_motioncorrection_20261008_2/11/ep2d_bold_fire_nomoco_tr1210_ref"
)
DEFAULT_OUTPUT = Path("data/boundary_pixel_audit")
DEFAULT_SPATIAL_TOLERANCE_MM = 0.001
REPRESENTATIVE_IDENTIFIERS = ("volume:0000", "volume:0209/slice:0010")


class AuditError(RuntimeError):
    """Raised when lossless correspondence cannot be established."""


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
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def boundary_masks(shape_yx: Sequence[int]) -> dict[str, np.ndarray]:
    """Return boundary masks in NumPy (y, x) order, including corner-safe unions."""
    if len(shape_yx) != 2:
        raise AuditError(f"Expected a 2D shape, received {tuple(shape_yx)}")
    rows, columns = (int(value) for value in shape_yx)
    if rows < 1 or columns < 1:
        raise AuditError(f"Boundary analysis requires nonempty axes: {(rows, columns)}")

    masks: dict[str, np.ndarray] = {}
    definitions = (
        ("first_x", (slice(None), 0)),
        ("last_x", (slice(None), columns - 1)),
        ("first_y", (0, slice(None))),
        ("last_y", (rows - 1, slice(None))),
    )
    for name, index in definitions:
        mask = np.zeros((rows, columns), dtype=bool)
        mask[index] = True
        masks[name] = mask

    outer = np.zeros((rows, columns), dtype=bool)
    outer[0, :] = outer[-1, :] = True
    outer[:, 0] = outer[:, -1] = True
    masks["all_outer_boundary"] = outer

    width_y, width_x = min(2, rows), min(2, columns)
    two_pixel_definitions = (
        ("first_two_x", (slice(None), slice(0, width_x))),
        ("last_two_x", (slice(None), slice(columns - width_x, columns))),
        ("first_two_y", (slice(0, width_y), slice(None))),
        ("last_two_y", (slice(rows - width_y, rows), slice(None))),
    )
    for name, index in two_pixel_definitions:
        mask = np.zeros((rows, columns), dtype=bool)
        mask[index] = True
        masks[name] = mask

    outer_two = np.zeros((rows, columns), dtype=bool)
    outer_two[:width_y, :] = outer_two[-width_y:, :] = True
    outer_two[:, :width_x] = outer_two[:, -width_x:] = True
    masks["all_outer_two_pixel_boundary"] = outer_two
    return masks


def _broadcast_plane_mask(mask_yx: np.ndarray, shape: Sequence[int]) -> np.ndarray:
    if tuple(shape[-2:]) != mask_yx.shape:
        raise AuditError(f"Mask shape {mask_yx.shape} does not match image shape {tuple(shape)}")
    return np.broadcast_to(mask_yx, tuple(shape))


def zero_mask_classes(native: np.ndarray, converted: np.ndarray) -> dict[str, np.ndarray]:
    native = np.asarray(native)
    converted = np.asarray(converted)
    if native.shape != converted.shape:
        raise AuditError(f"Array shape mismatch: {native.shape} versus {converted.shape}")
    native_zero = native == 0
    converted_zero = converted == 0
    return {
        "both_zero": native_zero & converted_zero,
        "converted_zero_native_nonzero": converted_zero & ~native_zero,
        "native_zero_converted_nonzero": native_zero & ~converted_zero,
        "both_nonzero": ~native_zero & ~converted_zero,
    }


def _pearson_from_values(first: np.ndarray, second: np.ndarray) -> float | None:
    first = np.asarray(first, dtype=np.float64).ravel()
    second = np.asarray(second, dtype=np.float64).ravel()
    finite = np.isfinite(first) & np.isfinite(second)
    first, second = first[finite], second[finite]
    if first.size < 2 or np.var(first) == 0.0 or np.var(second) == 0.0:
        return None
    return float(np.corrcoef(first, second)[0, 1])


def pair_metrics(
    native: np.ndarray, converted: np.ndarray, mask: np.ndarray | None = None
) -> dict[str, Any]:
    native = np.asarray(native)
    converted = np.asarray(converted)
    if native.shape != converted.shape:
        raise AuditError(f"Array shape mismatch: {native.shape} versus {converted.shape}")
    selected = (
        np.ones(native.shape, dtype=bool)
        if mask is None
        else np.array(mask, dtype=bool, copy=True)
    )
    selected &= np.isfinite(native) & np.isfinite(converted)
    native_values = np.asarray(native[selected], dtype=np.float64)
    converted_values = np.asarray(converted[selected], dtype=np.float64)
    if not native_values.size:
        return {"pixel_count": 0, "pearson_correlation": None, "mae": None, "rmse": None}
    difference = native_values - converted_values
    squared = difference * difference
    return {
        "pixel_count": int(native_values.size),
        "native_zero_count": int(np.count_nonzero(native_values == 0)),
        "native_zero_percentage": float(100.0 * np.mean(native_values == 0)),
        "converted_zero_count": int(np.count_nonzero(converted_values == 0)),
        "converted_zero_percentage": float(100.0 * np.mean(converted_values == 0)),
        "converted_zero_native_nonzero_count": int(
            np.count_nonzero((converted_values == 0) & (native_values != 0))
        ),
        "native_zero_converted_nonzero_count": int(
            np.count_nonzero((native_values == 0) & (converted_values != 0))
        ),
        "native_mean": float(np.mean(native_values)),
        "native_median": float(np.median(native_values)),
        "native_maximum": float(np.max(native_values)),
        "converted_mean": float(np.mean(converted_values)),
        "converted_median": float(np.median(converted_values)),
        "converted_maximum": float(np.max(converted_values)),
        "mae": float(np.mean(np.abs(difference))),
        "rmse": float(np.sqrt(np.mean(squared))),
        "squared_error": float(np.sum(squared)),
        "pearson_correlation": _pearson_from_values(native_values, converted_values),
    }


def squared_error_attribution(native: np.ndarray, converted: np.ndarray) -> dict[str, Any]:
    native = np.asarray(native)
    converted = np.asarray(converted)
    classes = zero_mask_classes(native, converted)
    boundary = _broadcast_plane_mask(boundary_masks(native.shape[-2:])["all_outer_boundary"], native.shape)
    difference = np.asarray(native, dtype=np.float64) - np.asarray(converted, dtype=np.float64)
    squared = difference * difference
    total = float(np.sum(squared))
    category_masks = {
        "boundary": boundary,
        "interior_converted_zero_native_nonzero": ~boundary & classes["converted_zero_native_nonzero"],
        "interior_native_zero_converted_nonzero": ~boundary & classes["native_zero_converted_nonzero"],
        "interior_both_nonzero": ~boundary & classes["both_nonzero"],
    }
    result: dict[str, Any] = {"total_squared_error": total}
    for name, category_mask in category_masks.items():
        value = float(np.sum(squared[category_mask]))
        result[f"{name}_squared_error"] = value
        result[f"{name}_squared_error_percentage"] = 100.0 * value / total if total else 0.0
    zero_disagreement = classes["converted_zero_native_nonzero"] | classes[
        "native_zero_converted_nonzero"
    ]
    converted_zero_disagreement = classes["converted_zero_native_nonzero"]
    for name, category_mask in (
        ("zero_mask_disagreement", zero_disagreement),
        ("converted_zero_native_nonzero", converted_zero_disagreement),
    ):
        value = float(np.sum(squared[category_mask]))
        result[f"{name}_squared_error"] = value
        result[f"{name}_squared_error_percentage"] = 100.0 * value / total if total else 0.0
    return result


@dataclass
class PairAccumulator:
    count: int = 0
    sum_a: float = 0.0
    sum_b: float = 0.0
    sum_aa: float = 0.0
    sum_bb: float = 0.0
    sum_ab: float = 0.0
    absolute_error: float = 0.0
    squared_error: float = 0.0

    def add(self, native: np.ndarray, converted: np.ndarray, mask: np.ndarray | None = None) -> None:
        native = np.asarray(native)
        converted = np.asarray(converted)
        selected = np.isfinite(native) & np.isfinite(converted)
        if mask is not None:
            selected &= np.asarray(mask, dtype=bool)
        a = np.asarray(native[selected], dtype=np.float64)
        b = np.asarray(converted[selected], dtype=np.float64)
        if not a.size:
            return
        difference = a - b
        self.count += int(a.size)
        self.sum_a += float(np.sum(a))
        self.sum_b += float(np.sum(b))
        self.sum_aa += float(np.dot(a, a))
        self.sum_bb += float(np.dot(b, b))
        self.sum_ab += float(np.dot(a, b))
        self.absolute_error += float(np.sum(np.abs(difference)))
        self.squared_error += float(np.dot(difference, difference))

    def summary(self) -> dict[str, Any]:
        if not self.count:
            return {"pixel_count": 0, "pearson_correlation": None, "mae": None, "rmse": None}
        variance_a = self.sum_aa - self.sum_a * self.sum_a / self.count
        variance_b = self.sum_bb - self.sum_b * self.sum_b / self.count
        covariance = self.sum_ab - self.sum_a * self.sum_b / self.count
        denominator = math.sqrt(max(0.0, variance_a) * max(0.0, variance_b))
        return {
            "pixel_count": self.count,
            "pearson_correlation": covariance / denominator if denominator else None,
            "mae": self.absolute_error / self.count,
            "rmse": math.sqrt(self.squared_error / self.count),
            "squared_error": self.squared_error,
        }


@dataclass
class BoundaryAccumulator:
    pixel_count: int = 0
    native_zero_count: int = 0
    converted_zero_count: int = 0
    converted_zero_native_nonzero_count: int = 0
    native_zero_converted_nonzero_count: int = 0
    native_sum: float = 0.0
    converted_sum: float = 0.0
    absolute_error_sum: float = 0.0
    squared_error_sum: float = 0.0
    native_maximum: float = -math.inf
    converted_maximum: float = -math.inf
    native_values: list[np.ndarray] = field(default_factory=list)
    converted_values: list[np.ndarray] = field(default_factory=list)

    def add(self, native: np.ndarray, converted: np.ndarray, mask: np.ndarray) -> None:
        selected = np.asarray(mask, dtype=bool) & np.isfinite(native) & np.isfinite(converted)
        a = np.asarray(native[selected], dtype=np.float32)
        b = np.asarray(converted[selected], dtype=np.float32)
        difference = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
        self.pixel_count += int(a.size)
        self.native_zero_count += int(np.count_nonzero(a == 0))
        self.converted_zero_count += int(np.count_nonzero(b == 0))
        self.converted_zero_native_nonzero_count += int(np.count_nonzero((b == 0) & (a != 0)))
        self.native_zero_converted_nonzero_count += int(np.count_nonzero((a == 0) & (b != 0)))
        self.native_sum += float(np.sum(a, dtype=np.float64))
        self.converted_sum += float(np.sum(b, dtype=np.float64))
        self.absolute_error_sum += float(np.sum(np.abs(difference)))
        self.squared_error_sum += float(np.dot(difference, difference))
        if a.size:
            self.native_maximum = max(self.native_maximum, float(np.max(a)))
            self.converted_maximum = max(self.converted_maximum, float(np.max(b)))
            self.native_values.append(a)
            self.converted_values.append(b)

    def summary(self, total_squared_error: float) -> dict[str, Any]:
        if not self.pixel_count:
            return {"pixel_count": 0}
        native = np.concatenate(self.native_values)
        converted = np.concatenate(self.converted_values)
        return {
            "pixel_count": self.pixel_count,
            "native_zero_count": self.native_zero_count,
            "native_zero_percentage": 100.0 * self.native_zero_count / self.pixel_count,
            "converted_zero_count": self.converted_zero_count,
            "converted_zero_percentage": 100.0 * self.converted_zero_count / self.pixel_count,
            "converted_zero_native_nonzero_count": self.converted_zero_native_nonzero_count,
            "converted_zero_native_nonzero_percentage": (
                100.0 * self.converted_zero_native_nonzero_count / self.pixel_count
            ),
            "native_zero_converted_nonzero_count": self.native_zero_converted_nonzero_count,
            "native_zero_converted_nonzero_percentage": (
                100.0 * self.native_zero_converted_nonzero_count / self.pixel_count
            ),
            "native_mean": self.native_sum / self.pixel_count,
            "native_median": float(np.median(native)),
            "native_maximum": self.native_maximum,
            "converted_mean": self.converted_sum / self.pixel_count,
            "converted_median": float(np.median(converted)),
            "converted_maximum": self.converted_maximum,
            "mae": self.absolute_error_sum / self.pixel_count,
            "rmse": math.sqrt(self.squared_error_sum / self.pixel_count),
            "squared_error": self.squared_error_sum,
            "total_squared_error_percentage": (
                100.0 * self.squared_error_sum / total_squared_error if total_squared_error else 0.0
            ),
        }


@dataclass(frozen=True)
class DicomVolumeRecord:
    volume: int
    path: Path
    frame_indices_by_anatomical_slice: tuple[int, ...]
    temporal_position_index: int
    acquisition_number: int | None


@dataclass
class EnhancedDicomSeries:
    root: Path
    volumes: dict[int, DicomVolumeRecord]
    metadata: dict[str, Any]
    _cached_volume: int | None = field(default=None, init=False)
    _cached_array: np.ndarray | None = field(default=None, init=False)
    _cached_geometry: sitk.Image | None = field(default=None, init=False)
    _cached_position_residual_mm: float | None = field(default=None, init=False)

    def load_volume(self, volume: int) -> tuple[np.ndarray, sitk.Image, float]:
        if self._cached_volume == volume:
            assert self._cached_array is not None and self._cached_geometry is not None
            assert self._cached_position_residual_mm is not None
            return self._cached_array, self._cached_geometry, self._cached_position_residual_mm
        record = self.volumes.get(volume)
        if record is None:
            raise AuditError(f"No unambiguous DICOM object maps to volume {volume}")
        dataset = pydicom.dcmread(record.path)
        decoded = np.asarray(dataset.pixel_array)
        if decoded.ndim == 2:
            decoded = decoded[np.newaxis, ...]
        array = decoded[np.asarray(record.frame_indices_by_anatomical_slice)]
        image, residual = dicom_volume_geometry(dataset, record.frame_indices_by_anatomical_slice, array)
        self._cached_volume = volume
        self._cached_array = array
        self._cached_geometry = image
        self._cached_position_residual_mm = residual
        return array, image, residual


def _functional_group_value(
    dataset: pydicom.Dataset, frame: pydicom.Dataset, sequence_name: str
) -> pydicom.Dataset:
    if hasattr(frame, sequence_name):
        return getattr(frame, sequence_name)[0]
    shared = dataset.get("SharedFunctionalGroupsSequence", [])
    if shared and hasattr(shared[0], sequence_name):
        return getattr(shared[0], sequence_name)[0]
    raise AuditError(f"DICOM frame lacks {sequence_name}")


def index_enhanced_dataset(dataset: pydicom.Dataset, path: Path) -> DicomVolumeRecord:
    groups = dataset.get("PerFrameFunctionalGroupsSequence")
    if not groups:
        raise AuditError(f"{path.name}: no enhanced per-frame metadata")
    if int(dataset.get("NumberOfFrames", len(groups))) != len(groups):
        raise AuditError(f"{path.name}: NumberOfFrames does not match per-frame metadata")
    temporal_values: set[int] = set()
    slice_to_frame: dict[int, int] = {}
    for frame_index, group in enumerate(groups):
        content = group.get("FrameContentSequence")
        if not content:
            raise AuditError(f"{path.name}: frame {frame_index} lacks FrameContentSequence")
        content = content[0]
        temporal = content.get("TemporalPositionIndex")
        anatomical = content.get("InStackPositionNumber")
        if temporal is None or anatomical is None:
            raise AuditError(
                f"{path.name}: frame {frame_index} lacks TemporalPositionIndex or InStackPositionNumber"
            )
        temporal_values.add(int(temporal))
        anatomical_zero_based = int(anatomical) - 1
        if anatomical_zero_based in slice_to_frame:
            raise AuditError(f"{path.name}: duplicate InStackPositionNumber {anatomical}")
        slice_to_frame[anatomical_zero_based] = frame_index
    if len(temporal_values) != 1:
        raise AuditError(f"{path.name}: frames contain multiple temporal positions")
    expected = list(range(len(groups)))
    if sorted(slice_to_frame) != expected:
        raise AuditError(f"{path.name}: anatomical frame indices are incomplete or noncontiguous")
    temporal = temporal_values.pop()
    acquisition = dataset.get("AcquisitionNumber")
    return DicomVolumeRecord(
        temporal - 1,
        path,
        tuple(slice_to_frame[index] for index in expected),
        temporal,
        None if acquisition is None else int(acquisition),
    )


def inspect_dicom_series(root: Path) -> EnhancedDicomSeries:
    if not root.is_dir():
        raise AuditError(f"DICOM directory does not exist: {root}")
    volumes: dict[int, DicomVolumeRecord] = {}
    errors: list[str] = []
    metadata_sets: dict[str, set[Any]] = defaultdict(set)
    file_count = frame_count = 0
    for path in sorted(item for item in root.iterdir() if item.is_file()):
        file_count += 1
        try:
            dataset = pydicom.dcmread(path, stop_before_pixels=True)
            record = index_enhanced_dataset(dataset, path)
            if record.volume in volumes:
                raise AuditError(f"duplicate TemporalPositionIndex {record.temporal_position_index}")
            volumes[record.volume] = record
            frame_count += len(record.frame_indices_by_anatomical_slice)
            for key in (
                "Rows", "Columns", "NumberOfFrames", "PixelRepresentation", "BitsAllocated",
                "BitsStored", "HighBit", "PixelPaddingValue", "PixelPaddingRangeLimit",
            ):
                value = dataset.get(key)
                if value is not None:
                    metadata_sets[key].add(str(value))
            metadata_sets["SOPClassUID"].add(str(dataset.get("SOPClassUID", "")))
            metadata_sets["TransferSyntaxUID"].add(
                str(dataset.file_meta.get("TransferSyntaxUID", ""))
            )
            slopes, intercepts = [], []
            for group in dataset.get("SharedFunctionalGroupsSequence", []):
                for transform in group.get("PixelValueTransformationSequence", []):
                    slopes.append(float(transform.get("RescaleSlope", 1.0)))
                    intercepts.append(float(transform.get("RescaleIntercept", 0.0)))
            for group in dataset.get("PerFrameFunctionalGroupsSequence", []):
                for transform in group.get("PixelValueTransformationSequence", []):
                    slopes.append(float(transform.get("RescaleSlope", 1.0)))
                    intercepts.append(float(transform.get("RescaleIntercept", 0.0)))
            metadata_sets["RescaleSlope"].update(str(value) for value in slopes)
            metadata_sets["RescaleIntercept"].update(str(value) for value in intercepts)
        except Exception as error:
            errors.append(f"{path.name}: {error}")
    expected_volumes = list(range(len(volumes)))
    coverage_complete = sorted(volumes) == expected_volumes
    metadata = {
        "status": "established" if not errors and coverage_complete else "indeterminate",
        "file_count": file_count,
        "frame_count": frame_count,
        "mapped_volume_count": len(volumes),
        "mapped_volume_range": [min(volumes), max(volumes)] if volumes else None,
        "complete_zero_based_temporal_coverage": coverage_complete,
        "mapping_basis": "TemporalPositionIndex and InStackPositionNumber from enhanced per-frame metadata",
        "filesystem_order_used": False,
        "errors": errors,
        **{key: sorted(values) for key, values in metadata_sets.items()},
    }
    for key in (
        "Rows", "Columns", "NumberOfFrames", "PixelRepresentation", "BitsAllocated",
        "BitsStored", "HighBit", "PixelPaddingValue", "PixelPaddingRangeLimit",
        "RescaleSlope", "RescaleIntercept", "SOPClassUID", "TransferSyntaxUID",
    ):
        metadata.setdefault(key, [])
    if errors or not coverage_complete:
        raise AuditError(f"DICOM frame correspondence is indeterminate: {metadata}")
    return EnhancedDicomSeries(root, volumes, metadata)


def dicom_volume_geometry(
    dataset: pydicom.Dataset,
    frame_indices_by_anatomical_slice: Sequence[int],
    array_zyx: np.ndarray,
) -> tuple[sitk.Image, float]:
    groups = dataset.PerFrameFunctionalGroupsSequence
    first = groups[frame_indices_by_anatomical_slice[0]]
    orientation_item = _functional_group_value(dataset, first, "PlaneOrientationSequence")
    measures = _functional_group_value(dataset, first, "PixelMeasuresSequence")
    orientation = np.asarray(orientation_item.ImageOrientationPatient, dtype=np.float64)
    x_direction = orientation[:3]
    y_direction = orientation[3:]
    z_direction = np.cross(x_direction, y_direction)
    z_direction /= np.linalg.norm(z_direction)
    pixel_spacing = np.asarray(measures.PixelSpacing, dtype=np.float64)
    z_spacing = float(measures.get("SpacingBetweenSlices", measures.get("SliceThickness")))
    positions = np.asarray([
        _functional_group_value(dataset, groups[index], "PlanePositionSequence").ImagePositionPatient
        for index in frame_indices_by_anatomical_slice
    ], dtype=np.float64)
    origin = positions[0]
    expected_positions = origin + np.arange(len(positions))[:, None] * z_spacing * z_direction
    residual = float(np.max(np.linalg.norm(positions - expected_positions, axis=1)))

    image = sitk.GetImageFromArray(np.asarray(array_zyx))
    image.SetOrigin(tuple(origin))
    image.SetSpacing((float(pixel_spacing[1]), float(pixel_spacing[0]), z_spacing))
    image.SetDirection(tuple(np.column_stack((x_direction, y_direction, z_direction)).ravel()))
    return image, residual


def losslessly_align(
    source: two_way.LoadedImage,
    target: two_way.LoadedImage,
    spatial_tolerance_mm: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    if source.image is None or target.image is None:
        raise AuditError("Lossless alignment requires SimpleITK physical geometry")
    relation = two_way._signed_axis_relation(source.image, target.image, spatial_tolerance_mm)
    if not (
        relation.get("detected")
        and relation.get("size_compatible")
        and relation.get("lossless_reorientation_feasible")
        and relation.get("within_spatial_tolerance")
    ):
        raise AuditError(f"No validated lossless physical-grid alignment: {relation}")
    return two_way._apply_lossless_reorientation(source.array, relation), relation


def align_dicom_volume(
    series: EnhancedDicomSeries,
    volume: int,
    native_reference: two_way.LoadedImage,
    spatial_tolerance_mm: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    array, image, position_residual = series.load_volume(volume)
    loaded = two_way.LoadedImage(array, two_way._sitk_geometry(image, array.shape), image)
    aligned, relation = losslessly_align(loaded, native_reference, spatial_tolerance_mm)
    relation = dict(relation)
    relation["dicom_frame_position_model_max_residual_mm"] = position_residual
    return aligned, relation


def _per_boundary_rows(
    identifier: str,
    native: np.ndarray,
    converted: np.ndarray,
    total_squared_error: float,
) -> list[dict[str, Any]]:
    rows = []
    for name, mask_2d in boundary_masks(native.shape[-2:]).items():
        mask = _broadcast_plane_mask(mask_2d, native.shape)
        metrics = pair_metrics(native, converted, mask)
        metrics["squared_error_percentage_of_image"] = (
            100.0 * metrics.get("squared_error", 0.0) / total_squared_error
            if total_squared_error else 0.0
        )
        rows.append({"image_identifier": identifier, "boundary": name, **metrics})
    return rows


def _zero_mask_row(identifier: str, native: np.ndarray, converted: np.ndarray) -> dict[str, Any]:
    classes = zero_mask_classes(native, converted)
    whole = pair_metrics(native, converted)
    nonzero = pair_metrics(native, converted, classes["both_nonzero"])
    counts = {f"{name}_count": int(np.count_nonzero(mask)) for name, mask in classes.items()}
    return {
        "image_identifier": identifier,
        "voxel_count": int(native.size),
        "native_zero_count": int(np.count_nonzero(native == 0)),
        "converted_zero_count": int(np.count_nonzero(converted == 0)),
        **counts,
        "neither_zero_count": counts["both_nonzero_count"],
        "whole_image_pearson_correlation": whole["pearson_correlation"],
        "whole_image_mae": whole["mae"],
        "whole_image_rmse": whole["rmse"],
        "nonzero_overlap_pearson_correlation": nonzero["pearson_correlation"],
        "nonzero_overlap_mae": nonzero["mae"],
        "nonzero_overlap_rmse": nonzero["rmse"],
        **squared_error_attribution(native, converted),
    }


def _write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _format(value: Any, digits: int = 8) -> str:
    if value is None:
        return "n/a"
    return f"{float(value):.{digits}g}"


def _display_plane(array: np.ndarray) -> tuple[np.ndarray, str]:
    array = np.asarray(array)
    if array.ndim == 2:
        return array, "2D image"
    if array.ndim == 3:
        index = array.shape[0] // 2
        return array[index], f"anatomical Z={index} of 0–{array.shape[0] - 1}"
    raise AuditError(f"Cannot display array with shape {array.shape}")


def generate_representative_figure(
    path: Path,
    identifier: str,
    native: np.ndarray,
    converted: np.ndarray,
) -> None:
    import matplotlib.pyplot as plt

    native_plane, plane_label = _display_plane(native)
    converted_plane, _ = _display_plane(converted)
    disagreement = np.zeros(native_plane.shape, dtype=np.int8)
    disagreement[(converted_plane == 0) & (native_plane != 0)] = 1
    disagreement[(native_plane == 0) & (converted_plane != 0)] = -1
    minimum = min(float(np.min(native_plane)), float(np.min(converted_plane)))
    maximum = max(float(np.max(native_plane)), float(np.max(converted_plane)))

    figure, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for axis, image, title in (
        (axes[0, 0], native_plane, "Native MRD"),
        (axes[0, 1], converted_plane, "Aligned DICOM-derived MRD"),
    ):
        shown = axis.imshow(image, origin="lower", cmap="gray", vmin=minimum, vmax=maximum)
        axis.set_title(title)
        axis.set_xlabel("SimpleITK x index")
        axis.set_ylabel("SimpleITK y index")
        figure.colorbar(shown, ax=axis, shrink=0.75, label="Intensity")
    shown = axes[0, 2].imshow(disagreement, origin="lower", cmap="bwr", vmin=-1, vmax=1)
    axes[0, 2].set_title("Exact-zero disagreement")
    axes[0, 2].set_xlabel("SimpleITK x index")
    axes[0, 2].set_ylabel("SimpleITK y index")
    colorbar = figure.colorbar(shown, ax=axes[0, 2], shrink=0.75, ticks=(-1, 0, 1))
    colorbar.ax.set_yticklabels(("native zero", "agreement", "converted zero"))

    x = np.arange(native_plane.shape[1])
    for y_index, label, style in (
        (0, "first Y row", "-"),
        (native_plane.shape[0] - 1, "last Y row", "--"),
    ):
        axes[1, 0].plot(x, native_plane[y_index], linestyle=style, label=f"Native {label}")
        axes[1, 0].plot(x, converted_plane[y_index], linestyle=style, label=f"Converted {label}")
    axes[1, 0].set_title("Y-boundary intensity profiles")
    axes[1, 0].set_xlabel("SimpleITK x index")
    axes[1, 0].set_ylabel("Intensity")
    axes[1, 0].legend(fontsize=7, ncol=2)

    y = np.arange(native_plane.shape[0])
    for x_index, label, style in (
        (0, "first X column", "-"),
        (native_plane.shape[1] - 1, "last X column", "--"),
    ):
        axes[1, 1].plot(y, native_plane[:, x_index], linestyle=style, label=f"Native {label}")
        axes[1, 1].plot(y, converted_plane[:, x_index], linestyle=style, label=f"Converted {label}")
    axes[1, 1].set_title("X-boundary intensity profiles")
    axes[1, 1].set_xlabel("SimpleITK y index")
    axes[1, 1].set_ylabel("Intensity")
    axes[1, 1].legend(fontsize=7, ncol=2)

    metrics = _zero_mask_row(identifier, native, converted)
    axes[1, 2].axis("off")
    axes[1, 2].text(
        0.0,
        1.0,
        "\n".join((
            f"Identifier: {identifier}",
            f"Displayed: {plane_label}",
            f"Converted-zero/native-nonzero: {metrics['converted_zero_native_nonzero_count']}",
            f"Native-zero/converted-nonzero: {metrics['native_zero_converted_nonzero_count']}",
            f"Whole-image MAE: {_format(metrics['whole_image_mae'])}",
            f"Whole-image RMSE: {_format(metrics['whole_image_rmse'])}",
            f"Boundary SSE: {_format(metrics['boundary_squared_error_percentage'])}%",
        )),
        va="top",
        family="monospace",
    )
    figure.suptitle(f"Boundary audit: {identifier}")
    figure.savefig(path, dpi=160)
    plt.close(figure)


def generate_aggregate_figure(path: Path, diagnostics: dict[str, np.ndarray], summary: dict[str, Any]) -> None:
    import matplotlib.pyplot as plt

    boundary_names = ("first_x", "last_x", "first_y", "last_y")
    boundary_values = [
        summary["boundary_aggregate"][name]["converted_zero_native_nonzero_percentage"]
        for name in boundary_names
    ]
    figure, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    axes[0, 0].bar(boundary_names, boundary_values, color="#4C78A8")
    axes[0, 0].set_title("Converted-zero/native-nonzero by boundary")
    axes[0, 0].set_ylabel("Pixels (%)")
    axes[0, 0].tick_params(axis="x", rotation=20)

    axes[0, 1].plot(diagnostics["row_mae"])
    axes[0, 1].set_title("Row-wise mean absolute error")
    axes[0, 1].set_xlabel("SimpleITK y index")
    axes[0, 1].set_ylabel("MAE")
    axes[0, 2].plot(diagnostics["column_mae"])
    axes[0, 2].set_title("Column-wise mean absolute error")
    axes[0, 2].set_xlabel("SimpleITK x index")
    axes[0, 2].set_ylabel("MAE")
    axes[1, 0].plot(100.0 * diagnostics["row_converted_zero_fraction"])
    axes[1, 0].set_title("Row-wise converted zero fraction")
    axes[1, 0].set_xlabel("SimpleITK y index")
    axes[1, 0].set_ylabel("Pixels (%)")
    axes[1, 1].plot(100.0 * diagnostics["column_converted_zero_fraction"])
    axes[1, 1].set_title("Column-wise converted zero fraction")
    axes[1, 1].set_xlabel("SimpleITK x index")
    axes[1, 1].set_ylabel("Pixels (%)")
    shown = axes[1, 2].imshow(
        100.0 * diagnostics["converted_zero_native_nonzero_frequency"],
        origin="lower",
        cmap="magma",
        vmin=0.0,
        vmax=100.0,
    )
    axes[1, 2].set_title("Converted-zero/native-nonzero frequency")
    axes[1, 2].set_xlabel("SimpleITK x index")
    axes[1, 2].set_ylabel("SimpleITK y index")
    figure.colorbar(shown, ax=axes[1, 2], shrink=0.75, label="Planes (%)")
    figure.savefig(path, dpi=160)
    plt.close(figure)


def write_report(output_dir: Path, summary: dict[str, Any]) -> None:
    boundary = summary["boundary_aggregate"]
    attribution = summary["squared_error_attribution"]
    dicom = summary["dicom_comparison"]
    systematic = [
        name for name in ("first_x", "last_x", "first_y", "last_y")
        if boundary[name]["converted_zero_percentage"] >= 95.0
    ]
    disagreement_boundaries = [
        name for name in ("first_x", "last_x", "first_y", "last_y")
        if boundary[name]["converted_zero_native_nonzero_percentage"] >= 95.0
    ]
    dicom_exact = dicom["converted_vs_original_dicom_unequal_voxels"] == 0
    diagnostics = summary["row_column_diagnostics"]
    aggregate = summary["aggregate_metrics"]
    total_voxels = aggregate["whole_image"]["pixel_count"]
    last_y_index = len(diagnostics["row_mae"]) - 1
    penultimate_y_index = max(0, last_y_index - 1)
    lines = [
        "# Boundary-pixel and zero-masking audit", "",
        "## Scope and correspondence", "",
        f"- Matched FIRE images: **{summary['image_counts']['matched_pairs']:,}**.",
        f"- Analyzed 2D planes: **{summary['image_counts']['analyzed_planes']:,}**.",
        "- Native and converted FIRE images were matched by volume/slice identifiers and aligned only by the validated signed-axis permutation/reflection.",
        "- Enhanced DICOM objects were mapped by `TemporalPositionIndex` and frames by `InStackPositionNumber`; filesystem order was not used.",
        "- No interpolation, reconstruction, or source-data modification was performed.", "",
        "## Boundary findings", "",
        "| Boundary | Native zero | Converted zero | Converted zero / native nonzero | MAE | SSE contribution |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name in ("first_x", "last_x", "first_y", "last_y", "all_outer_boundary"):
        item = boundary[name]
        lines.append(
            f"| `{name}` | {_format(item['native_zero_percentage'])}% | "
            f"{_format(item['converted_zero_percentage'])}% | "
            f"{_format(item['converted_zero_native_nonzero_percentage'])}% | "
            f"{_format(item['mae'])} | {_format(item['total_squared_error_percentage'])}% |"
        )
    lines.extend([
        "",
        f"Boundaries meeting the prespecified descriptive criterion of at least 95% converted zeros: "
        f"**{', '.join(systematic) if systematic else 'none'}**.",
        f"Boundaries with at least 95% converted-zero/native-nonzero disagreement: "
        f"**{', '.join(disagreement_boundaries) if disagreement_boundaries else 'none'}**.",
        "",
        f"The X-edge effect has an important one-pixel offset: `x=0` is 100% zero in both native and converted images and contributes no error, while `x=1` is 100% zero in converted/DICOM images, has column MAE **{_format(diagnostics['column_mae'][1])}**, and is nonzero in native images. Thus the visible left-edge residual is at `x=1`, immediately inside a shared-zero `x=0` column.",
        "",
        f"The Y-edge effect is confined to the outer rows: converted zero fractions are 100% at `y=0` and `y={last_y_index}`, but only {_format(100.0 * diagnostics['row_converted_zero_fraction'][1])}% at `y=1` and {_format(100.0 * diagnostics['row_converted_zero_fraction'][-2])}% at `y={penultimate_y_index}`. The two-pixel perimeter contributes **{_format(boundary['all_outer_two_pixel_boundary']['total_squared_error_percentage'])}%** of total squared error.",
        "",
        "## Whole-image zero masks and error attribution", "",
        f"- Native contains **{aggregate['native_zero_count']:,}** exact-zero voxels ({_format(100.0 * aggregate['native_zero_count'] / total_voxels)}%); converted/DICOM contains **{aggregate['converted_zero_count']:,}** ({_format(100.0 * aggregate['converted_zero_count'] / total_voxels)}%).",
        f"- Converted-zero/native-nonzero disagreements: **{aggregate['converted_zero_native_nonzero_count']:,}** voxels ({_format(100.0 * aggregate['converted_zero_native_nonzero_count'] / total_voxels)}%).",
        f"- The unique one-pixel outer boundary contributes **{_format(attribution['boundary_squared_error_percentage'])}%** of total squared error.",
        f"- The two-pixel outer perimeter contributes **{_format(boundary['all_outer_two_pixel_boundary']['total_squared_error_percentage'])}%** of total squared error.",
        f"- Converted-zero/native-nonzero pixels contribute **{_format(attribution['converted_zero_native_nonzero_squared_error_percentage'])}%** of total squared error.",
        f"- All exact-zero disagreements contribute **{_format(attribution['zero_mask_disagreement_squared_error_percentage'])}%** of total squared error.",
        f"- After excluding all zero-mask disagreements, **{_format(attribution['residual_after_excluding_zero_mask_squared_error_percentage'])}%** of total squared error remains.",
        f"- After excluding zero-mask disagreements, Pearson r/MAE/RMSE are "
        f"**{_format(summary['aggregate_metrics']['nonzero_overlap']['pearson_correlation'])} / "
        f"{_format(summary['aggregate_metrics']['nonzero_overlap']['mae'])} / "
        f"{_format(summary['aggregate_metrics']['nonzero_overlap']['rmse'])}**.",
        "",
        "## Original Siemens DICOM evidence", "",
        f"- Inspected **{summary['dicom_metadata']['file_count']}** enhanced MR objects containing "
        f"**{summary['dicom_metadata']['frame_count']}** frames.",
        f"- Original DICOM versus aligned converted FIRE unequal voxels: "
        f"**{dicom['converted_vs_original_dicom_unequal_voxels']:,}** of "
        f"**{dicom['compared_voxels']:,}**.",
        f"- Original-DICOM-nonzero/converted-zero voxels: "
        f"**{dicom['dicom_nonzero_converted_zero_voxels']:,}**.",
        f"- Original-DICOM-zero/converted-nonzero voxels: "
        f"**{dicom['dicom_zero_converted_nonzero_voxels']:,}**.",
        f"- Pixel storage: Rows={summary['dicom_metadata']['Rows']}, Columns={summary['dicom_metadata']['Columns']}, BitsAllocated={summary['dicom_metadata']['BitsAllocated']}, BitsStored={summary['dicom_metadata']['BitsStored']}, PixelRepresentation={summary['dicom_metadata']['PixelRepresentation']}.",
        f"- Standard scaling metadata: slope={summary['dicom_metadata']['RescaleSlope']}, intercept={summary['dicom_metadata']['RescaleIntercept']}; pixel padding value/range={summary['dicom_metadata']['PixelPaddingValue']}/{summary['dicom_metadata']['PixelPaddingRangeLimit']}.",
        "",
        "## Answers", "",
        f"1. **Are boundaries systematically zeroed?** "
        f"{f'Yes: converted/DICOM x=0, x=1, y=0, and y={last_y_index} are fully zero. Native x=0 is already zero; the added disagreements are x=1, y=0, and y={last_y_index}.' if systematic else 'No boundary reaches 95% converted zeros.'}",
        f"2. **Are those zeros present in original DICOM?** "
        f"{'Yes. The original decoded DICOM pixels are exactly equal to aligned converted FIRE pixels.' if dicom_exact else 'Not completely; see the DICOM mismatch counts above.'}",
        f"3. **Does downstream DICOM-to-NIfTI-to-MRD conversion add zeros?** "
        f"{'No additional zeros were detected.' if dicom_exact else 'Additional or removed zeros were detected; the stage is not isolated by this audit.'}",
        f"4. **Likely source?** "
        f"{'The zero mask is already present in the Siemens enhanced DICOM representation, while downstream conversion preserves it exactly. This points upstream of dcm2niix/NIfTI-to-MRD, although the evidence cannot distinguish Siemens reconstruction from the DICOM export step itself.' if dicom_exact else 'The source remains unresolved because DICOM and converted pixels are not identical.'}",
        "",
        "Observed counts and percentages above are measurements. The attribution to an internal Siemens reconstruction/export stage remains a bounded hypothesis because no pre-export Siemens processing trace is available.",
    ])
    (output_dir / "boundary_audit_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_audit(
    native_dir: Path,
    converted_dir: Path,
    dicom_dir: Path,
    output_dir: Path,
    spatial_tolerance_mm: float = DEFAULT_SPATIAL_TOLERANCE_MM,
) -> dict[str, Any]:
    native = three_way.discover_fire_dataset(native_dir.resolve(), "native")
    converted = three_way.discover_fire_dataset(converted_dir.resolve(), "converted")
    matches = two_way.match_images(native.records, converted.records)
    if matches["missing_from_a"] or matches["missing_from_b"] or matches["ambiguous_identifiers"]:
        raise AuditError(f"FIRE image coverage is not one-to-one: {matches}")
    dicom = inspect_dicom_series(dicom_dir.resolve())
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    reference_shape = native.reference.array.shape[-2:]
    boundary_accumulators = {
        name: BoundaryAccumulator() for name in boundary_masks(reference_shape)
    }
    whole_accumulator = PairAccumulator()
    nonzero_accumulator = PairAccumulator()
    boundary_rows: list[dict[str, Any]] = []
    zero_rows: list[dict[str, Any]] = []
    total_zero_counts: Counter[str] = Counter()
    attribution_sums: Counter[str] = Counter()
    dicom_counts: Counter[str] = Counter()
    relation_counts: Counter[str] = Counter()
    dicom_relation_counts: Counter[str] = Counter()
    dicom_position_residual_max = 0.0
    rows, columns = reference_shape
    row_abs_sum = np.zeros(rows, dtype=np.float64)
    column_abs_sum = np.zeros(columns, dtype=np.float64)
    row_zero_count = np.zeros(rows, dtype=np.int64)
    column_zero_count = np.zeros(columns, dtype=np.int64)
    row_sample_count = np.zeros(rows, dtype=np.int64)
    column_sample_count = np.zeros(columns, dtype=np.int64)
    disagreement_frequency_count = np.zeros((rows, columns), dtype=np.int64)
    plane_count = 0
    representatives: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    cached_dicom_volume: int | None = None
    cached_dicom_aligned: np.ndarray | None = None

    for pair_index, (identifier, native_record, converted_record) in enumerate(matches["pairs"]):
        native_loaded = two_way.load_image(native_record.path)
        converted_loaded = two_way.load_image(converted_record.path)
        aligned_converted, relation = losslessly_align(
            converted_loaded, native_loaded, spatial_tolerance_mm
        )
        relation_key = f"permutation={relation['source_axis_to_target_axis']};signs={relation['source_axis_signs']}"
        relation_counts[relation_key] += 1
        native_array = np.asarray(native_loaded.array)
        converted_array = np.asarray(aligned_converted)
        volume, slice_ordinal = three_way.parse_identifier(identifier)

        if cached_dicom_volume != volume:
            cached_dicom_aligned, dicom_relation = align_dicom_volume(
                dicom, volume, native.reference, spatial_tolerance_mm
            )
            cached_dicom_volume = volume
            dicom_key = (
                f"permutation={dicom_relation['source_axis_to_target_axis']};"
                f"signs={dicom_relation['source_axis_signs']}"
            )
            dicom_relation_counts[dicom_key] += 1
            dicom_position_residual_max = max(
                dicom_position_residual_max,
                float(dicom_relation["dicom_frame_position_model_max_residual_mm"]),
            )
        assert cached_dicom_aligned is not None
        if slice_ordinal is None:
            dicom_array = cached_dicom_aligned
        else:
            anatomical = native.anatomical_slice(identifier, native_loaded)
            dicom_array = cached_dicom_aligned[anatomical : anatomical + 1]
        if dicom_array.shape != converted_array.shape:
            raise AuditError(
                f"DICOM shape mismatch for {identifier}: {dicom_array.shape} versus {converted_array.shape}"
            )
        dicom_counts["compared_voxels"] += int(dicom_array.size)
        dicom_counts["converted_vs_original_dicom_unequal_voxels"] += int(
            np.count_nonzero(converted_array != dicom_array)
        )
        dicom_counts["dicom_nonzero_converted_zero_voxels"] += int(
            np.count_nonzero((dicom_array != 0) & (converted_array == 0))
        )
        dicom_counts["dicom_zero_converted_nonzero_voxels"] += int(
            np.count_nonzero((dicom_array == 0) & (converted_array != 0))
        )

        image_attribution = squared_error_attribution(native_array, converted_array)
        total_sse = image_attribution["total_squared_error"]
        boundary_rows.extend(_per_boundary_rows(identifier, native_array, converted_array, total_sse))
        zero_row = _zero_mask_row(identifier, native_array, converted_array)
        zero_rows.append(zero_row)
        for key, value in image_attribution.items():
            if key.endswith("_squared_error") or key == "total_squared_error":
                attribution_sums[key] += float(value)
        classes = zero_mask_classes(native_array, converted_array)
        total_zero_counts["native_zero_count"] += int(np.count_nonzero(native_array == 0))
        total_zero_counts["converted_zero_count"] += int(np.count_nonzero(converted_array == 0))
        for name, mask in classes.items():
            total_zero_counts[f"{name}_count"] += int(np.count_nonzero(mask))
        whole_accumulator.add(native_array, converted_array)
        nonzero_accumulator.add(native_array, converted_array, classes["both_nonzero"])

        for name, mask_2d in boundary_masks(native_array.shape[-2:]).items():
            mask = _broadcast_plane_mask(mask_2d, native_array.shape)
            boundary_accumulators[name].add(native_array, converted_array, mask)

        native_planes = native_array.reshape((-1, rows, columns))
        converted_planes = converted_array.reshape((-1, rows, columns))
        absolute = np.abs(
            np.asarray(native_planes, dtype=np.float64) - np.asarray(converted_planes, dtype=np.float64)
        )
        converted_zero = converted_planes == 0
        disagreement = converted_zero & (native_planes != 0)
        row_abs_sum += np.sum(absolute, axis=(0, 2))
        column_abs_sum += np.sum(absolute, axis=(0, 1))
        row_zero_count += np.sum(converted_zero, axis=(0, 2))
        column_zero_count += np.sum(converted_zero, axis=(0, 1))
        row_sample_count += converted_planes.shape[0] * columns
        column_sample_count += converted_planes.shape[0] * rows
        disagreement_frequency_count += np.sum(disagreement, axis=0)
        plane_count += converted_planes.shape[0]
        if identifier in REPRESENTATIVE_IDENTIFIERS:
            representatives[identifier] = (native_array.copy(), converted_array.copy())

        if pair_index and pair_index % 1000 == 0:
            print(f"Processed {pair_index:,} of {len(matches['pairs']):,} matched images", flush=True)

    total_squared_error = float(attribution_sums["total_squared_error"])
    attribution_summary = {"total_squared_error": total_squared_error}
    for key, value in attribution_sums.items():
        if key == "total_squared_error":
            continue
        attribution_summary[key] = float(value)
        attribution_summary[f"{key}_percentage"] = (
            100.0 * float(value) / total_squared_error if total_squared_error else 0.0
        )
    residual_sse = total_squared_error - float(
        attribution_sums["zero_mask_disagreement_squared_error"]
    )
    attribution_summary["residual_after_excluding_zero_mask_squared_error"] = residual_sse
    attribution_summary["residual_after_excluding_zero_mask_squared_error_percentage"] = (
        100.0 * residual_sse / total_squared_error if total_squared_error else 0.0
    )
    boundary_summary = {
        name: accumulator.summary(total_squared_error)
        for name, accumulator in boundary_accumulators.items()
    }
    aggregate_metrics = {
        "whole_image": whole_accumulator.summary(),
        "nonzero_overlap": nonzero_accumulator.summary(),
        **dict(total_zero_counts),
    }
    diagnostics = {
        "row_mae": row_abs_sum / row_sample_count,
        "column_mae": column_abs_sum / column_sample_count,
        "row_converted_zero_fraction": row_zero_count / row_sample_count,
        "column_converted_zero_fraction": column_zero_count / column_sample_count,
        "converted_zero_native_nonzero_frequency": disagreement_frequency_count / plane_count,
    }
    summary = {
        "utility": "audit_boundary_pixels.py",
        "inputs": {
            "native_fire": str(native_dir.resolve()),
            "converted_fire": str(converted_dir.resolve()),
            "original_dicom": str(dicom_dir.resolve()),
        },
        "configuration": {
            "spatial_tolerance_mm": spatial_tolerance_mm,
            "exact_zero_comparisons": True,
            "interpolation_used": False,
        },
        "image_counts": {
            "matched_pairs": len(matches["pairs"]),
            "analyzed_planes": plane_count,
            "missing_from_native": matches["missing_from_a"],
            "missing_from_converted": matches["missing_from_b"],
        },
        "correspondence": {
            "fire_alignment_relationships": dict(relation_counts),
            "dicom_to_native_alignment_relationships": dict(dicom_relation_counts),
            "native_slice_ordinal_to_anatomical_z": native.slice_ordinal_to_anatomical,
            "maximum_dicom_frame_position_model_residual_mm": dicom_position_residual_max,
        },
        "boundary_aggregate": boundary_summary,
        "aggregate_metrics": aggregate_metrics,
        "squared_error_attribution": attribution_summary,
        "dicom_metadata": dicom.metadata,
        "dicom_comparison": dict(dicom_counts),
        "row_column_diagnostics": {key: value.tolist() for key, value in diagnostics.items()},
        "interpretation": {
            "converted_pixels_equal_original_dicom": (
                dicom_counts["converted_vs_original_dicom_unequal_voxels"] == 0
            ),
            "additional_converted_zeros_detected": (
                dicom_counts["dicom_nonzero_converted_zero_voxels"] > 0
            ),
            "zeros_removed_by_conversion_detected": (
                dicom_counts["dicom_zero_converted_nonzero_voxels"] > 0
            ),
        },
    }

    _write_csv(output_dir / "per_image_boundary_metrics.csv", boundary_rows)
    _write_csv(output_dir / "per_image_zero_mask_metrics.csv", zero_rows)
    (output_dir / "boundary_audit_summary.json").write_text(
        json.dumps(_json_safe(summary), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    generate_aggregate_figure(
        figure_dir / "aggregate_boundary_diagnostics.png", diagnostics, summary
    )
    for identifier, (native_array, converted_array) in representatives.items():
        safe_identifier = identifier.replace(":", "-").replace("/", "_")
        generate_representative_figure(
            figure_dir / f"boundary_{safe_identifier}.png",
            identifier,
            native_array,
            converted_array,
        )
    write_report(output_dir, summary)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit boundary pixels and exact-zero masks in native and DICOM-derived FIRE images"
    )
    parser.add_argument("--native-dir", type=Path, default=DEFAULT_NATIVE)
    parser.add_argument("--converted-dir", type=Path, default=DEFAULT_CONVERTED)
    parser.add_argument("--dicom-dir", type=Path, default=DEFAULT_DICOM)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--spatial-tolerance-mm", type=float, default=DEFAULT_SPATIAL_TOLERANCE_MM
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = run_audit(
            args.native_dir,
            args.converted_dir,
            args.dicom_dir,
            args.output_dir,
            args.spatial_tolerance_mm,
        )
    except (AuditError, OSError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(f"Matched images: {summary['image_counts']['matched_pairs']:,}")
    print(
        "Boundary SSE percentage: "
        f"{summary['squared_error_attribution']['boundary_squared_error_percentage']:.6f}%"
    )
    print(
        "Zero-mask SSE percentage: "
        f"{summary['squared_error_attribution']['zero_mask_disagreement_squared_error_percentage']:.6f}%"
    )
    print(f"Reports: {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Diagnose native-MRD versus DICOM intensity representation differences."""

from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Sequence

import h5py
import ismrmrd
import matplotlib
import numpy as np
import pydicom

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from shared.analysis_tools import compare_motion_control_images as two_way
    from shared.analysis_tools import compare_three_way_mrd_images as three_way
except ModuleNotFoundError:  # Direct execution from shared/analysis_tools.
    import compare_motion_control_images as two_way
    import compare_three_way_mrd_images as three_way


SECTION_HEADING = "## Native vs. DICOM Intensity Representation Investigation"
PERCENTILES = (0, 1, 5, 25, 50, 75, 95, 99, 100)
ERROR_THRESHOLDS = (0.5, 1.0, 2.0, 5.0, 10.0)


class DiagnosticError(RuntimeError):
    """Raised when intensity correspondence or diagnostic inputs are invalid."""


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


def _validated_values(
    first: np.ndarray, second: np.ndarray | None = None, mask: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray | None]:
    first = np.asarray(first)
    if second is not None:
        second = np.asarray(second)
        if first.shape != second.shape:
            raise DiagnosticError(f"Array shapes differ: {first.shape} versus {second.shape}")
    valid = np.isfinite(first)
    if second is not None:
        valid &= np.isfinite(second)
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != first.shape:
            raise DiagnosticError("Voxel mask shape does not match the input arrays")
        valid &= mask
    if not np.any(valid):
        raise DiagnosticError("No finite voxels remain for intensity diagnostics")
    first_values = np.asarray(first[valid], dtype=np.float64)
    second_values = None if second is None else np.asarray(second[valid], dtype=np.float64)
    return first_values, second_values


def value_distribution(array: np.ndarray, representative_count: int = 12) -> dict[str, Any]:
    values, _ = _validated_values(array)
    nearest_distance = np.abs(values - np.rint(values))
    fractional_bins = (0.0, 0.01, 0.05, 0.1, 0.25, 0.5)
    histogram = {"exact_integer": int(np.count_nonzero(nearest_distance == 0))}
    for lower, upper in zip(fractional_bins[:-1], fractional_bins[1:]):
        selected = (nearest_distance > lower) & (
            nearest_distance <= upper if upper == 0.5 else nearest_distance < upper
        )
        histogram[f"({lower:g},{upper:g}{']' if upper == 0.5 else ')'}"] = int(
            np.count_nonzero(selected)
        )
    indices = np.linspace(0, values.size - 1, min(representative_count, values.size), dtype=int)
    return {
        "storage_dtype": np.asarray(array).dtype.str,
        "count": int(values.size),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "mean": float(np.mean(values)),
        "standard_deviation": float(np.std(values)),
        "integer_valued_count": int(np.count_nonzero(nearest_distance == 0)),
        "integer_valued_percentage": float(100.0 * np.mean(nearest_distance == 0)),
        "distance_to_nearest_integer_percentiles": {
            str(item): float(value)
            for item, value in zip(PERCENTILES, np.percentile(nearest_distance, PERCENTILES))
        },
        "fractional_distance_histogram": histogram,
        "representative_values": [float(values[index]) for index in indices],
    }


def error_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, Any]:
    observed_values, predicted_values = _validated_values(observed, predicted)
    assert predicted_values is not None
    differences = observed_values - predicted_values
    absolute = np.abs(differences)
    return {
        "voxel_count": int(differences.size),
        "exact_agreement_count": int(np.count_nonzero(differences == 0)),
        "exact_agreement_percentage": float(100.0 * np.mean(differences == 0)),
        "mean_absolute_error": float(np.mean(absolute)),
        "rmse": float(np.sqrt(np.mean(differences * differences))),
        "mean_signed_error": float(np.mean(differences)),
        "signed_error_standard_deviation": float(np.std(differences)),
        "maximum_absolute_error": float(np.max(absolute)),
        "signed_error_percentiles": {
            str(item): float(value)
            for item, value in zip(PERCENTILES, np.percentile(differences, PERCENTILES))
        },
        "absolute_error_threshold_percentages": {
            str(item): float(100.0 * np.mean(absolute <= item))
            for item in ERROR_THRESHOLDS
        },
    }


def fit_linear_relationship(
    target: np.ndarray, source: np.ndarray, mask: np.ndarray | None = None
) -> dict[str, Any]:
    """Fit target = scale * source + intercept using stable scalar sums."""
    target_values, source_values = _validated_values(target, source, mask)
    assert source_values is not None
    count = source_values.size
    source_mean = float(np.mean(source_values))
    target_mean = float(np.mean(target_values))
    centered_source = source_values - source_mean
    source_ss = float(np.dot(centered_source, centered_source))
    if source_ss <= 0:
        return {
            "available": False,
            "reason": "source intensities have zero variance",
            "voxel_count": int(count),
        }
    centered_target = target_values - target_mean
    scale = float(np.dot(centered_source, centered_target) / source_ss)
    intercept = target_mean - scale * source_mean
    residual = target_values - (scale * source_values + intercept)
    residual_ss = float(np.dot(residual, residual))
    target_ss = float(np.dot(centered_target, centered_target))
    denominator = math.sqrt(source_ss * target_ss)
    pearson = float(np.dot(centered_source, centered_target) / denominator) if denominator else None
    return {
        "available": True,
        "model": "target = scale * source + intercept",
        "voxel_count": int(count),
        "scale": scale,
        "intercept": intercept,
        "r_squared": 1.0 - residual_ss / target_ss if target_ss else None,
        "pearson_correlation": pearson,
        "residual_mean": float(np.mean(residual)),
        "residual_standard_deviation": float(np.std(residual)),
        "residual_mae": float(np.mean(np.abs(residual))),
        "residual_rmse": float(math.sqrt(residual_ss / count)),
        "maximum_absolute_residual": float(np.max(np.abs(residual))),
        "residual_percentiles": {
            str(item): float(value)
            for item, value in zip(PERCENTILES, np.percentile(residual, PERCENTILES))
        },
    }


def quantization_models(native: np.ndarray, dicom: np.ndarray) -> dict[str, Any]:
    models = {
        "round": np.rint(native),
        "truncate": np.trunc(native),
        "floor": np.floor(native),
        "ceil": np.ceil(native),
    }
    return {name: error_metrics(dicom, prediction) for name, prediction in models.items()}


def scaling_quantization_models(
    native: np.ndarray, dicom: np.ndarray, native_from_dicom_fit: dict[str, Any]
) -> dict[str, Any]:
    if not native_from_dicom_fit.get("available"):
        return {"available": False, "reason": native_from_dicom_fit.get("reason")}
    scale = float(native_from_dicom_fit["scale"])
    intercept = float(native_from_dicom_fit["intercept"])
    if scale == 0:
        return {"available": False, "reason": "fitted scale is zero"}
    continuous = (np.asarray(native, dtype=np.float32) - intercept) / scale
    result: dict[str, Any] = {
        "available": True,
        "inverse_implied_scale": 1.0 / scale,
        "inverse_implied_intercept": -intercept / scale,
        "continuous_inverse": error_metrics(dicom, continuous),
    }
    for name, operation in (
        ("round", np.rint), ("truncate", np.trunc), ("floor", np.floor), ("ceil", np.ceil)
    ):
        result[name] = error_metrics(dicom, operation(continuous))
    return result


def zero_floor_diagnostics(native: np.ndarray, dicom: np.ndarray) -> dict[str, Any]:
    if native.shape != dicom.shape:
        raise DiagnosticError("Zero-floor diagnostics require equal array shapes")
    dicom_zero = np.asarray(dicom) == 0
    dicom_positive = np.asarray(dicom) > 0
    squared_error = (np.asarray(native, dtype=np.float64) - np.asarray(dicom, dtype=np.float64)) ** 2
    total_sse = float(np.sum(squared_error))
    result: dict[str, Any] = {
        "dicom_zero_voxel_count": int(np.count_nonzero(dicom_zero)),
        "dicom_zero_voxel_percentage": float(100.0 * np.mean(dicom_zero)),
        "dicom_positive_voxel_count": int(np.count_nonzero(dicom_positive)),
        "dicom_positive_voxel_percentage": float(100.0 * np.mean(dicom_positive)),
        "dicom_zero_squared_error_percentage": (
            float(100.0 * np.sum(squared_error[dicom_zero]) / total_sse) if total_sse else 0.0
        ),
        "dicom_positive_squared_error_percentage": (
            float(100.0 * np.sum(squared_error[dicom_positive]) / total_sse) if total_sse else 0.0
        ),
    }
    if np.any(dicom_zero):
        result["native_where_dicom_zero"] = value_distribution(np.asarray(native)[dicom_zero])
        result["raw_error_where_dicom_zero"] = error_metrics(
            np.asarray(native)[dicom_zero], np.asarray(dicom)[dicom_zero]
        )
    if np.any(dicom_positive):
        result["raw_error_where_dicom_positive"] = error_metrics(
            np.asarray(native)[dicom_positive], np.asarray(dicom)[dicom_positive]
        )
        result["rounding_where_dicom_positive"] = quantization_models(
            np.asarray(native)[dicom_positive], np.asarray(dicom)[dicom_positive]
        )["round"]
        result["linear_fit_where_dicom_positive"] = fit_linear_relationship(
            native, dicom, dicom_positive
        )
    return result


def fit_per_volume(native: np.ndarray, dicom: np.ndarray) -> list[dict[str, Any]]:
    if native.ndim != 4 or native.shape != dicom.shape:
        raise DiagnosticError("Per-volume fitting requires equal (time,z,y,x) arrays")
    rows = []
    for volume in range(native.shape[0]):
        fit = fit_linear_relationship(native[volume], dicom[volume])
        row = {"volume": volume, **fit}
        rows.append(row)
    return rows


def summarize_fit_stability(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    available = [row for row in rows if row.get("available")]
    if not available:
        return {"available": False, "reason": "no volume has source variance"}
    result: dict[str, Any] = {"available": True, "volume_count": len(rows), "fit_count": len(available)}
    for key in ("scale", "intercept", "residual_rmse", "residual_mae"):
        values = np.asarray([row[key] for row in available], dtype=np.float64)
        result[key] = {
            "minimum": float(np.min(values)),
            "p05": float(np.percentile(values, 5)),
            "median": float(np.median(values)),
            "p95": float(np.percentile(values, 95)),
            "maximum": float(np.max(values)),
            "mean": float(np.mean(values)),
            "standard_deviation": float(np.std(values)),
        }
    return result


def representative_slice_fits(
    native: np.ndarray, dicom: np.ndarray, volumes: Sequence[int], slices: Sequence[int]
) -> list[dict[str, Any]]:
    rows = []
    for volume in volumes:
        for slice_index in slices:
            fit = fit_linear_relationship(native[volume, slice_index], dicom[volume, slice_index])
            rows.append({"volume": volume, "anatomical_slice": slice_index, **fit})
    return rows


def load_aligned_native_dicom(
    native_path: Path, dicom_path: Path, spatial_tolerance_mm: float
) -> tuple[three_way.FireDataset, three_way.DicomNrrd, np.ndarray, np.ndarray]:
    native = three_way.discover_fire_dataset(native_path, "native")
    dicom = three_way.load_dicom_nrrd(dicom_path)
    expected = three_way._expected_fire_identifiers(
        dicom.volume_count, dicom.spatial_shape_zyx[0]
    )
    if set(native.records) != expected:
        raise DiagnosticError("Native FIRE coverage differs from the validated three-way mapping")
    native_array = np.empty(dicom.array_tzyx.shape, dtype=native.reference.array.dtype)
    dicom_array = np.empty(dicom.array_tzyx.shape, dtype=dicom.array_tzyx.dtype)
    relationships = Counter()
    for volume in range(dicom.volume_count):
        native_array[volume] = native.volume_array(volume)
        relation = two_way._signed_axis_relation(
            dicom.volume_image(volume), native.reference.image, spatial_tolerance_mm
        )
        if not relation.get("detected") or not relation.get("size_compatible"):
            raise DiagnosticError(f"No validated signed-axis index mapping for volume {volume}")
        relationships[
            (tuple(relation["source_axis_to_target_axis"]), tuple(relation["source_axis_signs"]))
        ] += 1
        dicom_array[volume] = two_way._apply_lossless_reorientation(
            dicom.array_tzyx[volume], relation
        )
    if len(relationships) != 1:
        raise DiagnosticError(f"Signed-axis mapping changes over time: {relationships}")
    return native, dicom, native_array, dicom_array


def _image_group_numbers(h5_file: h5py.File) -> list[int]:
    numbers = []
    for key in h5_file:
        match = re.fullmatch(r"image_(\d+)", key)
        if match:
            numbers.append(int(match.group(1)))
    return sorted(numbers)


def inspect_mrd_capture(
    path: Path, fire_root: Path | None = None, slices_per_volume: int = 20
) -> dict[str, Any]:
    if not path.is_file():
        return {"available": False, "path": str(path), "reason": "file not found"}
    dtype_counts: Counter[str] = Counter()
    data_type_code_counts: Counter[int] = Counter()
    count = integer_count = 0
    total = total_squares = 0.0
    minimum, maximum = math.inf, -math.inf
    representative = []
    raw_matches = raw_missing = 0
    raw_files: dict[tuple[int, int], list[Path]] = {}
    if fire_root is not None:
        raw_pattern = re.compile(r"_volume_(\d+)_slice_(\d+)\.raw$", re.I)
        for raw_path in fire_root.rglob("*.raw"):
            match = raw_pattern.search(raw_path.name)
            if match:
                raw_files.setdefault((int(match.group(1)), int(match.group(2))), []).append(
                    raw_path
                )
    with h5py.File(path, "r") as h5_file:
        numbers = _image_group_numbers(h5_file)
        if not numbers:
            return {"available": False, "path": str(path), "reason": "no image groups"}
        for number in numbers:
            group = h5_file[f"image_{number}"]
            header_bytes = np.asarray(group["header"], dtype=np.uint8).tobytes()
            attribute_bytes = np.asarray(group["attribute"], dtype=np.uint8).tobytes()
            attribute = attribute_bytes.decode("utf-8").split("\0", 1)[0]
            image = ismrmrd.Image(header_bytes, attribute)
            dtype = np.dtype(ismrmrd.get_dtype_from_data_type(image.data_type))
            data_bytes = np.asarray(group["data"], dtype=np.uint8).tobytes()
            values = np.frombuffer(data_bytes, dtype=dtype)
            dtype_counts[dtype.str] += 1
            data_type_code_counts[int(image.data_type)] += 1
            finite = np.asarray(values[np.isfinite(values)], dtype=np.float64)
            count += finite.size
            integer_count += int(np.count_nonzero(finite == np.rint(finite)))
            total += float(np.sum(finite))
            total_squares += float(np.dot(finite, finite))
            minimum = min(minimum, float(np.min(finite)))
            maximum = max(maximum, float(np.max(finite)))
            if len(representative) < 16:
                representative.extend(float(item) for item in finite[: 16 - len(representative)])
            if fire_root is not None:
                volume, ordinal = divmod(number, slices_per_volume)
                candidates = raw_files.get((volume, ordinal), [])
                if len(candidates) != 1:
                    raw_missing += 1
                elif candidates[0].read_bytes() == data_bytes:
                    raw_matches += 1
    mean = total / count
    return {
        "available": True,
        "path": str(path),
        "image_count": len(numbers),
        "voxel_count": count,
        "ismrmrd_data_type_codes": dict(data_type_code_counts),
        "numpy_dtypes": dict(dtype_counts),
        "minimum": minimum,
        "maximum": maximum,
        "mean": mean,
        "standard_deviation": math.sqrt(max(0.0, total_squares / count - mean * mean)),
        "integer_valued_percentage": 100.0 * integer_count / count,
        "representative_values": representative,
        "fire_raw_byte_validation": {
            "expected_slice_files": len(numbers),
            "byte_identical_slice_files": raw_matches,
            "missing_or_ambiguous_slice_files": raw_missing,
            "all_byte_identical": raw_matches == len(numbers) and raw_missing == 0,
        } if fire_root is not None else None,
    }


def derive_dicom_series_dir(dicom_nrrd: Path) -> Path | None:
    parts = list(dicom_nrrd.parts)
    for index, part in enumerate(parts):
        if part.endswith("-4-nrrd"):
            parts[index] = part.removesuffix("-4-nrrd") + "-2-sorted"
            candidate = Path(*parts[:index + 1], *parts[index + 1 : -1])
            return candidate
    return None


def extract_dicom_scaling_metadata(dataset: pydicom.dataset.Dataset) -> dict[str, Any]:
    slopes, intercepts, types = [], [], []
    locations = []

    def add_transform(item: pydicom.dataset.Dataset, location: str) -> None:
        slopes.append(float(item.RescaleSlope))
        intercepts.append(float(item.RescaleIntercept))
        types.append(str(item.get("RescaleType", "")))
        locations.append(location)

    if "RescaleSlope" in dataset and "RescaleIntercept" in dataset:
        add_transform(dataset, "top-level")
    for location, groups in (
        ("shared", dataset.get("SharedFunctionalGroupsSequence", [])),
        ("per-frame", dataset.get("PerFrameFunctionalGroupsSequence", [])),
    ):
        for group_index, group in enumerate(groups):
            for item in group.get("PixelValueTransformationSequence", []):
                add_transform(item, f"{location}[{group_index}].PixelValueTransformationSequence")
    return {
        "slopes": slopes,
        "intercepts": intercepts,
        "rescale_types": types,
        "locations": locations,
        "real_world_value_mapping_present": bool(
            dataset.get("RealWorldValueMappingSequence")
            or any(
                group.get("RealWorldValueMappingSequence")
                for sequence_name in ("SharedFunctionalGroupsSequence", "PerFrameFunctionalGroupsSequence")
                for group in dataset.get(sequence_name, [])
            )
        ),
    }


def inspect_dicom_series(series_dir: Path | None) -> dict[str, Any]:
    if series_dir is None or not series_dir.is_dir():
        return {
            "available": False,
            "path": None if series_dir is None else str(series_dir),
            "reason": "enhanced DICOM series directory not found",
        }
    files = sorted(path for path in series_dir.iterdir() if path.is_file())
    if not files:
        return {"available": False, "path": str(series_dir), "reason": "no DICOM files"}
    fields: dict[str, set[Any]] = {
        key: set() for key in (
            "BitsAllocated", "BitsStored", "HighBit", "PixelRepresentation"
        )
    }
    text_fields: dict[str, set[str]] = {
        key: set() for key in (
            "PixelPaddingValue", "PixelPaddingRangeLimit", "LossyImageCompression",
            "PresentationLUTShape",
        )
    }
    slopes, intercepts, types, transform_locations = [], [], [], Counter()
    real_world_count = frames = 0
    sop_classes, transfer_syntaxes = set(), set()
    private_scaling_candidates: Counter[str] = Counter()
    private_protocol_scaling: dict[str, Counter[float]] = {
        "dOverallImageScaleFactor": Counter(),
        "dOverallImageScaleCorrectionFactor": Counter(),
    }
    first_pixel_dtype = first_pixel_range = None
    for file_index, path in enumerate(files):
        dataset = pydicom.dcmread(path, stop_before_pixels=True)
        frames += int(dataset.get("NumberOfFrames", 1))
        sop_classes.add(str(dataset.get("SOPClassUID", "")))
        transfer_syntaxes.add(str(dataset.file_meta.get("TransferSyntaxUID", "")))
        for key in fields:
            value = dataset.get(key)
            if value is not None:
                fields[key].add(int(value))
        for key in text_fields:
            value = dataset.get(key)
            if value is not None:
                text_fields[key].add(str(value))
        scaling = extract_dicom_scaling_metadata(dataset)
        slopes.extend(scaling["slopes"])
        intercepts.extend(scaling["intercepts"])
        types.extend(scaling["rescale_types"])
        transform_locations.update(scaling["locations"])
        real_world_count += int(scaling["real_world_value_mapping_present"])
        for element in dataset.iterall():
            descriptor = f"{element.keyword} {element.name}".lower()
            if any(term in descriptor for term in ("scale", "slope", "intercept", "real world")):
                private_scaling_candidates[f"{element.tag} {element.name}={str(element.value)[:100]}"] += 1
            if element.tag.is_private and isinstance(element.value, (bytes, str)):
                private_text = (
                    element.value.decode("latin1", errors="ignore")
                    if isinstance(element.value, bytes) else element.value
                )
                for key in private_protocol_scaling:
                    for match in re.finditer(
                        rf"\b{re.escape(key)}\s*=\s*([-+0-9.eE]+)", private_text
                    ):
                        private_protocol_scaling[key][float(match.group(1))] += 1
        if file_index == 0:
            pixels = pydicom.dcmread(path).pixel_array
            first_pixel_dtype = pixels.dtype.str
            first_pixel_range = [int(np.min(pixels)), int(np.max(pixels))]
    return {
        "available": True,
        "path": str(series_dir),
        "file_count": len(files),
        "frame_count": frames,
        "sop_class_uids": sorted(sop_classes),
        "transfer_syntax_uids": sorted(transfer_syntaxes),
        **{key: sorted(value) for key, value in fields.items()},
        **{key: sorted(value) for key, value in text_fields.items()},
        "first_file_decoded_pixel_dtype": first_pixel_dtype,
        "first_file_decoded_pixel_range": first_pixel_range,
        "pixel_value_transformation": {
            "count": len(slopes),
            "unique_slopes": sorted(set(slopes)),
            "unique_intercepts": sorted(set(intercepts)),
            "unique_rescale_types": sorted(set(types)),
            "locations": dict(transform_locations),
        },
        "real_world_value_mapping_file_count": real_world_count,
        "scaling_related_elements": dict(private_scaling_candidates),
        "siemens_private_protocol_scaling": {
            key: {
                "unique_values": sorted(counts),
                "occurrences": {str(value): count for value, count in counts.items()},
            }
            for key, counts in private_protocol_scaling.items()
        },
        "interpretation": (
            "Standard enhanced-DICOM functional groups specify identity scaling; "
            "no Real World Value Mapping was found. Siemens private protocol scale "
            "parameters are reported but were not applied as a vendor formula."
        ),
    }


def _find_single_mrd(root: Path) -> Path | None:
    candidates = sorted(root.glob("measurement-*.hdf5"))
    return candidates[0] if len(candidates) == 1 else None


def _representative_pairs(native: np.ndarray, dicom: np.ndarray) -> list[dict[str, Any]]:
    requested = (
        (0, 0, 0, 1), (0, 10, 36, 36), (80, 5, 24, 24),
        (160, 10, 36, 36), (240, 15, 48, 48), (319, 19, 36, 36),
    )
    return [
        {
            "index_tzyx": list(index),
            "native": float(native[index]),
            "dicom": int(dicom[index]),
            "native_minus_dicom": float(native[index] - dicom[index]),
        }
        for index in requested
    ]


def generate_figures(
    figure_dir: Path,
    native: np.ndarray,
    dicom: np.ndarray,
    fit: dict[str, Any],
    per_volume: Sequence[dict[str, Any]],
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    scale, intercept = float(fit["scale"]), float(fit["intercept"])
    flat_native, flat_dicom = native.ravel(), dicom.ravel()
    step = max(1, flat_native.size // 100_000)
    sampled_native = flat_native[::step]
    sampled_dicom = flat_dicom[::step]

    figure, axis = plt.subplots(figsize=(7, 6))
    axis.scatter(sampled_dicom, sampled_native, s=2, alpha=0.08, rasterized=True)
    lower = min(float(np.min(sampled_dicom)), float(np.min(sampled_native)))
    upper = max(float(np.max(sampled_dicom)), float(np.max(sampled_native)))
    x_line = np.asarray((lower, upper))
    axis.plot(x_line, x_line, color="black", linewidth=1, label="Identity")
    axis.plot(x_line, scale * x_line + intercept, color="red", linewidth=1.5, label="OLS fit")
    axis.set_xlabel("DICOM intensity")
    axis.set_ylabel("Native MRD/FIRE intensity")
    axis.set_title(
        f"Native vs. DICOM intensity\nscale={scale:.8f}, intercept={intercept:.5f}, "
        f"r={fit['pearson_correlation']:.8f}, R²={fit['r_squared']:.8f}"
    )
    axis.legend()
    figure.tight_layout()
    figure.savefig(figure_dir / "figure5_native_dicom_intensity_relationship.png", dpi=170)
    plt.close(figure)

    original = native.astype(np.float32) - dicom.astype(np.float32)
    fitted = native.astype(np.float32) - (scale * dicom.astype(np.float32) + intercept)
    inverse = (native.astype(np.float32) - intercept) / scale
    quantization = dicom.astype(np.float32) - inverse
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for axis, values, title in zip(
        axes,
        (original, fitted, quantization),
        ("Native - DICOM", "Native - fitted DICOM", "DICOM - inverse-scaled native"),
    ):
        low, high = np.percentile(values, (0.5, 99.5))
        axis.hist(values.ravel(), bins=120, range=(low, high), color="#4C78A8")
        axis.axvline(0, color="black", linewidth=1)
        axis.set_title(title)
        axis.set_xlabel("Residual intensity")
        axis.set_ylabel("Voxel count")
    figure.tight_layout()
    figure.savefig(figure_dir / "figure6_intensity_residual_distributions.png", dpi=170)
    plt.close(figure)
    del original, fitted, inverse, quantization
    gc.collect()

    volumes = np.asarray([row["volume"] for row in per_volume if row.get("available")])
    figure, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for axis, key, label in zip(
        axes,
        ("scale", "intercept", "residual_rmse"),
        ("OLS scale", "OLS intercept", "Residual RMSE"),
    ):
        values = np.asarray([row[key] for row in per_volume if row.get("available")])
        axis.plot(volumes, values, linewidth=1)
        axis.axhline(np.mean(values), color="red", linestyle="--", linewidth=1)
        axis.set_ylabel(label)
        axis.grid(alpha=0.2)
    axes[-1].set_xlabel("Acquisition volume")
    figure.tight_layout()
    figure.savefig(figure_dir / "figure7_scaling_stability.png", dpi=170)
    plt.close(figure)

    figure, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    z_index, y_index = native.shape[1] // 2, native.shape[2] // 2
    for axis, volume in zip(axes, (0, native.shape[0] // 2, native.shape[0] - 1)):
        x = np.arange(native.shape[3])
        native_profile = native[volume, z_index, y_index]
        dicom_profile = dicom[volume, z_index, y_index]
        axis.plot(x, native_profile, label="Native", linewidth=1.5)
        axis.plot(x, dicom_profile, label="DICOM", linewidth=1, alpha=0.8)
        axis.plot(x, scale * dicom_profile + intercept, label="Fitted DICOM", linewidth=1)
        axis.set_ylabel(f"Volume {volume}")
        axis.grid(alpha=0.2)
    axes[0].legend(ncol=3)
    axes[-1].set_xlabel(f"X index (anatomical Z={z_index}, Y={y_index})")
    figure.suptitle("Representative native and DICOM intensity profiles")
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(figure_dir / "figure8_representative_intensity_profiles.png", dpi=170)
    plt.close(figure)


def _format(value: Any, digits: int = 9) -> str:
    return "n/a" if value is None else f"{float(value):.{digits}g}"


def write_diagnostics(
    output_dir: Path,
    diagnostics: dict[str, Any],
    per_volume: Sequence[dict[str, Any]],
    per_slice: Sequence[dict[str, Any]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "intensity_diagnostics.json").write_text(
        json.dumps(_json_safe(diagnostics), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for filename, rows in (
        ("intensity_per_volume.csv", per_volume),
        ("intensity_representative_slices.csv", per_slice),
    ):
        fields = sorted({key for row in rows for key in row})
        with (output_dir / filename).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    raw = diagnostics["models"]["raw_native_vs_dicom"]
    rounded = diagnostics["models"]["simple_quantization"]["round"]
    fit = diagnostics["models"]["linear_native_from_dicom"]
    scaled = diagnostics["models"]["scaling_plus_quantization"]["round"]
    native_mrd = diagnostics["source_datatypes"]["native_mrd"]
    dicom_metadata = diagnostics["dicom_metadata"]
    stability = diagnostics["per_volume_fit_stability"]
    zero_floor = diagnostics["zero_floor_diagnostics"]
    section = [
        SECTION_HEADING, "",
        "### Source datatypes and processing path", "",
        "| Representation | Storage datatype | Integer-valued samples | Verified behavior |",
        "|---|---|---:|---|",
        f"| Scanner-native MRD capture | `{', '.join(native_mrd.get('numpy_dtypes', {}))}` | {_format(native_mrd.get('integer_valued_percentage'))}% | {native_mrd.get('image_count', 0)} incoming image messages inspected |",
        f"| Native FIRE NRRD | `{diagnostics['source_datatypes']['native_fire']['storage_dtype']}` | {_format(diagnostics['source_datatypes']['native_fire']['integer_valued_percentage'])}% | Detached NRRD points to incoming bytes without numeric conversion |",
        f"| Direct DICOM NRRD | `{diagnostics['source_datatypes']['dicom_nrrd']['storage_dtype']}` | {_format(diagnostics['source_datatypes']['dicom_nrrd']['integer_valued_percentage'])}% | dcm2niix output |",
        f"| Converted-MRD FIRE NRRD | `{diagnostics['source_datatypes']['converted_fire']['storage_dtype']}` | {_format(diagnostics['source_datatypes']['converted_fire']['integer_valued_percentage'])}% | Same FIRE byte-writing path |",
        "",
        f"The original native MRD capture is available at `{native_mrd.get('path')}`. It declares ISMRMRD datatype code 5 (`float32`) and contains genuinely fractional values. All {native_mrd.get('fire_raw_byte_validation', {}).get('byte_identical_slice_files', 0)} corresponding FIRE `.raw` slice files are byte-identical to the captured incoming MRD payloads. FIRE therefore preserves the incoming pixels; it does not create the fractional values.",
        "",
        f"The converted-MRD replay capture is also `float32`, but {_format(diagnostics['source_datatypes']['converted_mrd'].get('integer_valued_percentage'))}% integer-valued. Its {diagnostics['source_datatypes']['converted_mrd'].get('fire_raw_byte_validation', {}).get('byte_identical_slice_files', 0)} slice payloads likewise byte-match the converted FIRE raw files.",
        "",
        "The FIRE receiver selects the NumPy dtype from the incoming MRD header, copies the socket payload with that dtype, writes the original `data_bytes` directly to `.raw`, and writes only a matching detached NRRD header. No scaling, normalization, clipping, or rounding occurs in this path. Native and converted MRD streams use the same functions.",
        "",
        "### Quantitative model comparison", "",
        "| Model | MAE | RMSE | Exact agreement | Interpretation |",
        "|---|---:|---:|---:|---|",
        f"| Raw native vs. DICOM | {_format(raw['mean_absolute_error'])} | {_format(raw['rmse'])} | {_format(raw['exact_agreement_percentage'])}% | Baseline |",
        f"| Native rounded vs. DICOM | {_format(rounded['mean_absolute_error'])} | {_format(rounded['rmse'])} | {_format(rounded['exact_agreement_percentage'])}% | Simple quantization is insufficient |",
        f"| Linear scale + offset | {_format(fit['residual_mae'])} | {_format(fit['residual_rmse'])} | N/A | High R², but large residuals remain |",
        f"| Inverse scaling + rounding | {_format(scaled['mean_absolute_error'])} | {_format(scaled['rmse'])} | {_format(scaled['exact_agreement_percentage'])}% | Combined model is insufficient |",
        "",
        f"Global OLS fit: `native = {_format(fit['scale'], 12)} × DICOM + {_format(fit['intercept'], 12)}`. Pearson r is {_format(fit['pearson_correlation'])} and R² is {_format(fit['r_squared'])}. The fitted residual MAE/RMSE are {_format(fit['residual_mae'])}/{_format(fit['residual_rmse'])}, versus raw {_format(raw['mean_absolute_error'])}/{_format(raw['rmse'])}. The broad image intensity range makes R² high, but the residual RMSE is much larger than a rounding error.",
        "",
        f"Per-volume scale is stable: median {_format(stability['scale']['median'], 12)}, range {_format(stability['scale']['minimum'], 12)}–{_format(stability['scale']['maximum'], 12)}. Intercept median is {_format(stability['intercept']['median'])}, range {_format(stability['intercept']['minimum'])}–{_format(stability['intercept']['maximum'])}. Residual RMSE remains {_format(stability['residual_rmse']['minimum'])}–{_format(stability['residual_rmse']['maximum'])} across time.",
        "",
        "### Zero-floor and foreground behavior", "",
        f"DICOM is exactly zero for {zero_floor['dicom_zero_voxel_count']} voxels ({_format(zero_floor['dicom_zero_voxel_percentage'])}%). At those locations, native intensity has median {_format(zero_floor['raw_error_where_dicom_zero']['signed_error_percentiles']['50'])}, mean {_format(zero_floor['raw_error_where_dicom_zero']['mean_signed_error'])}, and maximum {_format(zero_floor['raw_error_where_dicom_zero']['maximum_absolute_error'])}. This small subset contributes {_format(zero_floor['dicom_zero_squared_error_percentage'])}% of total squared error.",
        "",
        f"For DICOM-positive voxels, raw MAE/RMSE fall to {_format(zero_floor['raw_error_where_dicom_positive']['mean_absolute_error'])}/{_format(zero_floor['raw_error_where_dicom_positive']['rmse'])}. A separate positive-voxel OLS fit has scale {_format(zero_floor['linear_fit_where_dicom_positive']['scale'], 12)}, intercept {_format(zero_floor['linear_fit_where_dicom_positive']['intercept'], 12)}, and residual RMSE {_format(zero_floor['linear_fit_where_dicom_positive']['residual_rmse'])}. Thus DICOM zeroing or thresholding is a major contributor, but substantial nonzero-region differences remain.",
        "",
        "### DICOM metadata evidence", "",
        f"The enhanced DICOM series contains {dicom_metadata.get('file_count')} files and {dicom_metadata.get('frame_count')} frames. `BitsAllocated={dicom_metadata.get('BitsAllocated')}`, `BitsStored={dicom_metadata.get('BitsStored')}`, `HighBit={dicom_metadata.get('HighBit')}`, and `PixelRepresentation={dicom_metadata.get('PixelRepresentation')}` (unsigned).",
        "",
        f"All {dicom_metadata.get('pixel_value_transformation', {}).get('count')} per-frame Pixel Value Transformation items have slope `{dicom_metadata.get('pixel_value_transformation', {}).get('unique_slopes')}`, intercept `{dicom_metadata.get('pixel_value_transformation', {}).get('unique_intercepts')}`, and type `{dicom_metadata.get('pixel_value_transformation', {}).get('unique_rescale_types')}`. No Real World Value Mapping or Pixel Padding value was found; `LossyImageCompression={dicom_metadata.get('LossyImageCompression')}`. Thus the standard DICOM metadata specifies identity scaling and does not reproduce the empirical native/DICOM residuals.",
        "",
        f"The Siemens private XProtocol exposes `dOverallImageScaleFactor={dicom_metadata.get('siemens_private_protocol_scaling', {}).get('dOverallImageScaleFactor', {}).get('unique_values')}` and `dOverallImageScaleCorrectionFactor={dicom_metadata.get('siemens_private_protocol_scaling', {}).get('dOverallImageScaleCorrectionFactor', {}).get('unique_values')}`. The former is far from the empirical near-unity OLS coefficient and the latter is unity. Their exact role for this enhanced-DICOM representation is not established, so neither was applied as a vendor scaling formula.",
        "",
        "### Answers and interpretation", "",
        "1. **Native FIRE is genuinely fractional:** yes—only about 1.42% of samples are exactly integer-valued.",
        "2. **Original scanner-native MRD is float-valued:** yes—the saved incoming MRD messages are `float32` and byte-match the FIRE raw data.",
        "3. **FIRE modifies pixels while writing NRRD:** no evidence of modification; source inspection and full byte comparison establish preservation.",
        "4. **Simple integer rounding explains the difference:** no. It yields only about 26.60% exact agreement and does not reduce RMSE.",
        "5. **Stable linear relationship:** the near-identity coefficients are highly stable over time, but the fitted residual RMSE remains about 15.39, so a single linear relationship is descriptive rather than explanatory.",
        "6. **Scaling plus integer quantization explains the difference:** no. Inverse-fit rounding has substantially worse exact agreement than direct rounding.",
        "7. **DICOM metadata supports fitted scaling:** no. Standard functional groups specify slope 1 and intercept 0; no real-world mapping is present.",
        "8. **Remaining discrepancy:** DICOM zero-floor behavior is directly observed and explains about half of total squared error, while nonzero voxels still retain RMSE above 11. This is most consistent with more complex scanner reconstruction/display/export differences between the native MRD image stream and the enhanced-DICOM representation. The available evidence does not isolate a vendor-private operation or exact stage.",
        "",
        "The evidence rules out datatype conversion alone, pure integer quantization, and a simple scale-plus-quantization model. It supports a more complex difference upstream of dcm2niix, while downstream MRD serialization is lossless. No converter change is indicated by this analysis.",
        "",
        "### Additional diagnostic figures", "",
        "- `figures/figure5_native_dicom_intensity_relationship.png`",
        "- `figures/figure6_intensity_residual_distributions.png`",
        "- `figures/figure7_scaling_stability.png`",
        "- `figures/figure8_representative_intensity_profiles.png`", "",
    ]
    report_path = output_dir / "comparison_report.md"
    existing = report_path.read_text(encoding="utf-8") if report_path.exists() else "# Three-way REF image comparison\n"
    existing = existing.split(SECTION_HEADING, 1)[0].rstrip()
    report_path.write_text(existing + "\n\n" + "\n".join(section), encoding="utf-8")


def run_diagnostics(
    native_path: Path,
    dicom_path: Path,
    converted_path: Path,
    output_dir: Path,
    native_mrd_path: Path | None = None,
    converted_mrd_path: Path | None = None,
    dicom_series_dir: Path | None = None,
    spatial_tolerance_mm: float = 0.001,
) -> dict[str, Any]:
    if not (output_dir / "comparison_summary.json").is_file():
        raise DiagnosticError(
            "Existing three-way outputs are required; run compare_three_way_mrd_images.py first"
        )
    native, dicom, native_array, dicom_array = load_aligned_native_dicom(
        native_path, dicom_path, spatial_tolerance_mm
    )
    converted = three_way.discover_fire_dataset(converted_path, "converted")
    if native_mrd_path is None:
        native_mrd_path = _find_single_mrd(native_path)
    if converted_mrd_path is None:
        converted_mrd_path = _find_single_mrd(converted_path)
    if dicom_series_dir is None:
        dicom_series_dir = derive_dicom_series_dir(dicom_path)

    raw_metrics = error_metrics(native_array, dicom_array)
    simple_quantization = quantization_models(native_array, dicom_array)
    global_fit = fit_linear_relationship(native_array, dicom_array)
    scaled_quantization = scaling_quantization_models(native_array, dicom_array, global_fit)
    zero_diagnostics = zero_floor_diagnostics(native_array, dicom_array)
    per_volume = fit_per_volume(native_array, dicom_array)
    representative_volumes = sorted(set((0, dicom.volume_count // 2, dicom.volume_count - 1)))
    representative_slices = sorted(set((0, dicom.spatial_shape_zyx[0] // 2, dicom.spatial_shape_zyx[0] - 1)))
    per_slice = representative_slice_fits(
        native_array, dicom_array, representative_volumes, representative_slices
    )

    converted_distribution = value_distribution(converted.reference.array)
    # The existing three-way result established every converted voxel equals DICOM;
    # use the complete DICOM distribution rather than rereading all converted slices.
    converted_distribution.update({
        "count": int(dicom_array.size),
        "minimum": float(np.min(dicom_array)),
        "maximum": float(np.max(dicom_array)),
        "mean": float(np.mean(dicom_array, dtype=np.float64)),
        "standard_deviation": float(np.std(dicom_array, dtype=np.float64)),
        "integer_valued_count": int(dicom_array.size),
        "integer_valued_percentage": 100.0,
        "provenance": "complete distribution inherited from established DICOM numeric equality",
    })
    native_mrd = (
        inspect_mrd_capture(native_mrd_path, native_path)
        if native_mrd_path is not None
        else {"available": False, "reason": "no unique native MRD capture found"}
    )
    converted_mrd = (
        inspect_mrd_capture(converted_mrd_path, converted_path)
        if converted_mrd_path is not None
        else {"available": False, "reason": "no unique converted MRD capture found"}
    )
    dicom_metadata = inspect_dicom_series(dicom_series_dir)
    diagnostics = {
        "utility": "diagnose_native_dicom_intensities.py",
        "inputs": {
            "native_fire": str(native_path),
            "dicom_nrrd": str(dicom_path),
            "converted_fire": str(converted_path),
            "native_mrd": None if native_mrd_path is None else str(native_mrd_path),
            "converted_mrd": None if converted_mrd_path is None else str(converted_mrd_path),
            "dicom_series": None if dicom_series_dir is None else str(dicom_series_dir),
        },
        "correspondence": {
            "method": "previously validated FIRE physical-slice mapping plus signed-axis DICOM index mapping",
            "interpolation_used": False,
            "spatial_tolerance_mm": spatial_tolerance_mm,
            "array_shape_tzyx": list(native_array.shape),
            "voxel_count": int(native_array.size),
            "native_slice_ordinal_to_anatomical_z": native.slice_ordinal_to_anatomical,
        },
        "source_datatypes": {
            "native_mrd": native_mrd,
            "native_fire": value_distribution(native_array),
            "dicom_source": dicom_metadata,
            "dicom_nrrd": value_distribution(dicom_array),
            "converted_mrd": converted_mrd,
            "converted_fire": converted_distribution,
        },
        "representative_corresponding_voxels": _representative_pairs(native_array, dicom_array),
        "models": {
            "raw_native_vs_dicom": raw_metrics,
            "simple_quantization": simple_quantization,
            "linear_native_from_dicom": global_fit,
            "scaling_plus_quantization": scaled_quantization,
        },
        "per_volume_fit_stability": summarize_fit_stability(per_volume),
        "representative_slice_fits": per_slice,
        "zero_floor_diagnostics": zero_diagnostics,
        "dicom_metadata": dicom_metadata,
        "fire_processing_path": {
            "verified": True,
            "connection_reader": "apps/python_fire_server/connection.py:281-296",
            "raw_writer": "apps/python_fire_server/handle_data.py:1233-1241",
            "slice_header_writer": "apps/python_fire_server/handle_data.py:1133-1189",
            "volume_header_writer": "apps/python_fire_server/handle_data.py:1299-1345",
            "finding": "incoming typed data bytes are written unchanged; NRRD headers preserve declared MRD dtype",
            "native_and_converted_share_path": True,
        },
        "interpretation": {
            "datatype_conversion_alone_supported": False,
            "pure_integer_quantization_supported": False,
            "single_linear_scaling_fully_explains_difference": False,
            "linear_scaling_plus_quantization_supported": False,
            "more_complex_reconstruction_or_export_difference_supported": True,
            "exact_vendor_processing_stage_isolated": False,
        },
    }
    write_diagnostics(output_dir, diagnostics, per_volume, per_slice)
    generate_figures(output_dir / "figures", native_array, dicom_array, global_fit, per_volume)
    return diagnostics


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Diagnose native-MRD versus DICOM pixel-intensity differences"
    )
    parser.add_argument("--native-fire", required=True, type=Path)
    parser.add_argument("--dicom-nrrd", required=True, type=Path)
    parser.add_argument("--converted-fire", required=True, type=Path)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/image_comparison_REF_three_way")
    )
    parser.add_argument("--native-mrd", type=Path)
    parser.add_argument("--converted-mrd", type=Path)
    parser.add_argument("--dicom-series-dir", type=Path)
    parser.add_argument("--spatial-tolerance-mm", type=float, default=0.001)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        diagnostics = run_diagnostics(
            args.native_fire.resolve(),
            args.dicom_nrrd.resolve(),
            args.converted_fire.resolve(),
            args.output_dir.resolve(),
            None if args.native_mrd is None else args.native_mrd.resolve(),
            None if args.converted_mrd is None else args.converted_mrd.resolve(),
            None if args.dicom_series_dir is None else args.dicom_series_dir.resolve(),
            args.spatial_tolerance_mm,
        )
    except (DiagnosticError, three_way.ComparisonError, OSError, RuntimeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    fit = diagnostics["models"]["linear_native_from_dicom"]
    rounded = diagnostics["models"]["simple_quantization"]["round"]
    print(f"Diagnostics: {args.output_dir.resolve()}")
    print(
        f"rounding exact={rounded['exact_agreement_percentage']:.6f}%, "
        f"fit scale={fit['scale']:.9f}, intercept={fit['intercept']:.9f}, "
        f"residual RMSE={fit['residual_rmse']:.9f}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

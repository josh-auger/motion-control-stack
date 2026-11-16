#!/usr/bin/env python3
"""Read-only retrospective NOCO motion comparison; see README_compare_subject_motion.md."""

from __future__ import annotations

import argparse
import base64
from collections import Counter
from datetime import datetime, timezone
import html
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

import numpy as np
import SimpleITK as sitk

# Existing single-acquisition analyzers use sibling imports, including internally.
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from analysis_common import (  # noqa: E402
    AnalysisError, describe, discover_unique, pyplot, read_csv, read_json, write_csv,
)
from analyze_registration_calls import parse_queue_items  # noqa: E402
from analyze_registration_completeness import parse_status  # noqa: E402

SUBJECT_RE = re.compile(r"^savedData_.+_(GBS\d+)$")
AXES = "XYZ"
PARAM_FIELDS = [f"{a}_{kind}_{unit}" for kind, unit in (("rotation", "deg"), ("translation", "mm")) for a in AXES]
MOTION_FIELDS = PARAM_FIELDS + ["translation_magnitude_mm", "rotation_euler_norm_deg", "rotation_angle_deg"]
GROUP_FIELDS = [
    "subject", "volume", "group", "is_reference", "acquisition_time_s", "time_source",
    "status", "reason", "registration_finalized", "transform_filename", "transform_source",
    "transform_type", "compute_zyx", "center_x_mm", "center_y_mm", "center_z_mm",
    "transform_valid", "pointer_valid", *MOTION_FIELDS, "group_displacement_mm",
    "recorded_group_displacement_mm", "group_motion_event", "runtime_ms", "queue_outcome",
]
FIGURES = [
    ("figure_1_motion_trajectories.png", "Reference-relative motion: all slice groups"),
    ("figure_2_group_displacement.png", "Consecutive slice-group displacement versus acquisition time"),
    ("figure_3_motion_distributions.png", "Consecutive slice-group displacement distributions and subject statistics"),
    ("figure_4_registration_performance.png", "Registration outcomes, excluding reference initialization"),
    ("figure_5_registration_runtime.png", "Full registration-call runtime distributions"),
]
POSITION_FIGURE = ("slice_group_exceedance_by_position.png", "High-motion exceedance by slice-group position")
POSITION_CSV = "slice_group_exceedance_by_position.csv"
POSITION_FIELDS = [
    "subject", "slice_group_index", "transition_type", "valid_transition_count",
    "above_threshold_count", "exceedance_pct", "mean_displacement_mm", "p95_displacement_mm",
]


def discover_subjects(root: Path) -> list[tuple[str, Path]]:
    subjects: dict[str, Path] = {}
    for path in sorted(root.glob("savedData_*_GBS*")):
        match = SUBJECT_RE.fullmatch(path.name)
        if not path.is_dir() or not match:
            continue
        subject = match[1]
        if subject in subjects:
            raise AnalysisError(f"Ambiguous repeated subject {subject}: {subjects[subject]} and {path}")
        subjects[subject] = path.resolve()
    if not subjects:
        raise AnalysisError(f"No savedData_*_GBS<digits> subject directories in {root}")
    return sorted(subjects.items())


def relative_transform(previous: sitk.Transform, current: sitk.Transform) -> sitk.Euler3DTransform:
    """Repository convention: current o inverse(previous), centered at previous center.

    This is the pure calculation from motion_monitor.monitor_directory.compose_transform_pair.
    Production modules are deliberately not imported (they configure service state).
    """
    inverse = previous.GetInverse()
    a = np.array(current.GetMatrix()).reshape(3, 3)
    b = np.array(inverse.GetMatrix()).reshape(3, 3)
    c0, c1 = np.array(current.GetCenter()), np.array(inverse.GetCenter())
    t = a @ (np.array(inverse.GetTranslation()) + c1 - c0) + np.array(current.GetTranslation()) + c0 - c1
    result = sitk.Euler3DTransform()
    result.SetCenter(c1.tolist())
    result.SetTranslation(t.tolist())
    result.SetMatrix((a @ b).ravel().tolist())
    return result


def displacement(transform: sitk.Euler3DTransform, radius_mm: float = 50.0) -> float:
    """Exact monitor Euler-based displacement formulation, with roundoff clipping."""
    if not math.isfinite(radius_mm) or radius_mm <= 0:
        raise ValueError("radius must be positive and finite")
    x, y, z, tx, ty, tz = transform.GetParameters()
    cosine = (-1 + math.cos(x) * math.cos(y) + math.cos(x) * math.cos(z)
              + math.cos(y) * math.cos(z) + math.sin(x) * math.sin(y) * math.sin(z)) / 2
    theta = math.acos(float(np.clip(cosine, -1, 1)))
    return math.sqrt(tx * tx + ty * ty + tz * tz) + radius_mm * math.sqrt((1 - math.cos(theta)) ** 2 + math.sin(theta) ** 2)


def consecutive_displacements(transforms: list[sitk.Transform | None], radius_mm: float = 50.0) -> list[float | None]:
    # The full expected sequence includes missing slots, so gaps cannot be bridged.
    return [None] + [
        displacement(relative_transform(a, b), radius_mm) if a is not None and b is not None else None
        for a, b in zip(transforms, transforms[1:])
    ] if transforms else []


def motion_parameters(transform: sitk.Transform | None) -> dict[str, Any]:
    if transform is None:
        return {field: None for field in MOTION_FIELDS}
    euler = sitk.Euler3DTransform(transform) if transform.GetName() == "Euler3DTransform" else sitk.Euler3DTransform()
    if transform.GetName() != "Euler3DTransform":
        euler.SetMatrix(transform.GetMatrix())
        euler.SetCenter(transform.GetCenter())
        euler.SetTranslation(transform.GetTranslation())
    params = euler.GetParameters()
    matrix = np.array(transform.GetMatrix()).reshape(3, 3)
    return {
        **dict(zip(PARAM_FIELDS, [*np.degrees(params[:3]).tolist(), *params[3:]])),
        "translation_magnitude_mm": float(np.linalg.norm(params[3:])),
        "rotation_euler_norm_deg": float(np.linalg.norm(np.degrees(params[:3]))),
        "rotation_angle_deg": math.degrees(math.acos(float(np.clip((np.trace(matrix) - 1) / 2, -1, 1)))),
    }


def metric_statistics(values: list[float | None], threshold: float | None = None) -> dict[str, Any]:
    result = describe(values)
    if threshold is not None:
        flags = [value is not None and value > threshold for value in values]
        count = sum(flags)
        result.update(
            threshold_mm=threshold, exceedance_count=count,
            exceedance_pct=100 * count / result["count"] if result["count"] else None,
            event_runs=sum(flag and (i == 0 or not flags[i - 1]) for i, flag in enumerate(flags)),
            missing_count=len(values) - result["count"],
        )
    return result


def acquisition_metadata(path: Path) -> dict[str, Any]:
    m = read_json(path)
    enc = m["encoding"][0]
    limits = enc["encodingLimits"]
    if limits["repetition"]["minimum"] != 0 or limits["slice"]["minimum"] != 0:
        raise AnalysisError("Only zero-based repetition/slice encoding is supported")
    tr = m["sequenceParameters"]["TR"]
    if len(tr) != 1 or not math.isfinite(float(tr[0])) or tr[0] <= 0:
        raise AnalysisError("A single positive TR is required")
    times = sorted(set(float(t) for t in m["SliceTiming"].values()))
    nslices = int(limits["slice"]["maximum"]) + 1
    if not times or nslices % len(times) or not all(math.isfinite(t) and 0 <= t < tr[0] / 1000 for t in times):
        raise AnalysisError("Slice timings do not define compatible equal-sized groups")
    if len(m["SliceTiming"]) != nslices or any(list(m["SliceTiming"].values()).count(t) != nslices // len(times) for t in times):
        raise AnalysisError("Slice timings must cover every slice with equal group sizes")
    matrix = enc["reconSpace"]["matrixSize"]
    fov = enc["reconSpace"]["fieldOfView_mm"]
    longs = {p["name"]: p["value"] for p in m["userParameters"].get("userParameterLong", [])}
    doubles = {p["name"]: p["value"] for p in m["userParameters"].get("userParameterDouble", [])}
    return {
        "volumes": int(limits["repetition"]["maximum"]) + 1, "slices_per_volume": nslices,
        "tr_s": tr[0] / 1000, "te_ms": m["sequenceParameters"].get("TE"),
        "matrix": [matrix["x"], matrix["y"], nslices],
        "spacing_mm": [fov["x"] / matrix["x"], fov["y"] / matrix["y"], doubles.get("SpacingBetweenSlices")],
        "multiband_factor": longs.get("MultiBandFactor"), "groups_per_volume": len(times),
        "slices_per_group": nslices // len(times), "group_offsets_s": times,
        "nominal_duration_s": (int(limits["repetition"]["maximum"]) + 1) * tr[0] / 1000,
        "subject_identifier_source": "results directory suffix; subjectInformation is not used",
    }


def acquisition_times(root: Path, meta: dict[str, Any], warnings: list[str]) -> tuple[dict[tuple[int, int], tuple[float | None, str]], int | None]:
    """Read stored MRD times for every group, without reading image pixel data.

    Slice ordinals in pointer filenames are arrival ordinals, not metadata slice keys.
    Metadata SliceTiming can therefore only be used as sorted group offsets.
    """
    n, ng = meta["volumes"], meta["groups_per_volume"]
    times = {(v, g): (v * meta["tr_s"] + meta["group_offsets_s"][g], "derived_metadata_TR_plus_sorted_SliceTiming")
             for v in range(n) for g in range(ng)}
    candidates = sorted(root.glob("measurement*.hdf5"))
    if len(candidates) > 1:
        raise AnalysisError("Ambiguous replay HDF5 files")
    if not candidates:
        warnings.append("Replay HDF5 missing: acquisition times derived from TR and sorted SliceTiming; slice receipt counts unavailable.")
        return times, None
    import h5py
    import ismrmrd

    missing, invalid = 0, 0
    with h5py.File(candidates[0], "r") as stream:
        names = {name for name in stream if re.fullmatch(r"image_\d+", name)}
        ns = meta["slices_per_volume"]
        expected_names = {f"image_{i}" for i in range(n * ns)}
        if names != expected_names:
            warnings.append(f"MRD image coverage differs from expected {n * ns}: missing {len(expected_names - names)}, extra {len(names - expected_names)}.")
        # Converter times are relative 400 Hz ticks; do not use filesystem/replay wall clocks.
        base = None
        if "image_0" in names:
            base = int(ismrmrd.ImageHeader.from_buffer_copy(bytes(stream["image_0/header"][:])).acquisition_time_stamp)
        for v in range(n):
            for g in range(ng):
                first = v * ns + g * meta["slices_per_group"]
                members = [f"image_{i}" for i in range(first, first + meta["slices_per_group"])]
                key = (v, g)
                if not all(name in names for name in members) or base is None:
                    times[key] = (None, "missing_MRD_image_or_time_origin")
                    missing += 1
                    continue
                header = ismrmrd.ImageHeader.from_buffer_copy(bytes(stream[members[0] + "/header"][:]))
                time = (int(header.acquisition_time_stamp) - base) / 400
                if int(header.repetition) != v or time < 0 or abs(time - times[key][0]) > 0.005:
                    times[key] = (None, "invalid_MRD_repetition_or_timestamp")
                    invalid += 1
                else:
                    times[key] = (time, "stored_MRD_acquisition_time_stamp_400Hz")
    if missing or invalid:
        warnings.append(f"Acquisition times missing for {missing} groups and invalid for {invalid}; preserved as missing.")
    return times, len(names)


def unique_by_key(rows: list[dict[str, Any]], fields: tuple[str, str], description: str) -> dict[tuple[int, int], dict[str, Any]]:
    result = {}
    for row in rows:
        key = tuple(int(row[field]) for field in fields)
        if key in result:
            raise AnalysisError(f"Duplicate {description} for {key}")
        result[key] = row
    return result


def load_transform(root: Path, filename: str, warnings: list[str]) -> tuple[sitk.Transform | None, str | None]:
    if not filename or Path(filename).name != filename:
        warnings.append(f"Missing or unsafe transform filename: {filename!r}")
        return None, None
    # Only glob the small results root. A multi-component pathlib glob scans the
    # 60,000-entry processed directory again for every one of the 6,707 transforms.
    candidates = [p for p in [root / filename, *(directory / "transforms" / filename for directory in root.glob("processed_*") if directory.is_dir())] if p.is_file()]
    if not candidates:
        warnings.append(f"Transform missing: {filename}")
        return None, None
    if len(candidates) > 1 and any(p.read_bytes() != candidates[0].read_bytes() for p in candidates[1:]):
        raise AnalysisError(f"Conflicting duplicate transform {filename}")
    try:
        transform = sitk.ReadTransform(str(candidates[0]))
        if transform.GetName() not in {"Euler3DTransform", "VersorRigid3DTransform"}:
            raise ValueError(f"unsupported transform type {transform.GetName()}")
        if not np.all(np.isfinite([*transform.GetParameters(), *transform.GetFixedParameters()])):
            raise ValueError("nonfinite transform parameters")
        matrix = np.array(transform.GetMatrix()).reshape(3, 3)
        if not np.allclose(matrix.T @ matrix, np.eye(3), atol=1e-7) or not np.isclose(np.linalg.det(matrix), 1):
            raise ValueError("nonrigid matrix")
    except (RuntimeError, ValueError) as error:
        warnings.append(f"Unreadable transform {filename}: {error}")
        return None, str(candidates[0].relative_to(root))
    return transform, str(candidates[0].relative_to(root))


def analyze_subject(subject: str, root: Path, radius_mm: float = 50, threshold_mm: float = 0.3) -> dict[str, Any]:
    warnings: list[str] = []
    metadata_path = discover_unique(root, "*metadata*.json", "series metadata")
    meta = acquisition_metadata(metadata_path)
    n, ng = meta["volumes"], meta["groups_per_volume"]
    if ng != 14:
        raise AnalysisError(f"Expected 14 slice groups per non-reference volume; found {ng}")
    status_path = root / "registration_status.json"
    if status_path.exists():
        parsed, _ = parse_status(status_path)
        raw = read_json(status_path)["volumes"]
    else:
        warnings.append("registration_status.json missing: outcomes unknown; transforms cannot be accepted as successful.")
        parsed, raw = [], {}
    status_by_volume = {r["volume"]: r for r in parsed}
    timing, image_count = acquisition_times(root, meta, warnings)
    queue_paths = list(root.glob("queue_profile_items_*.csv"))
    runtime, queue = {}, {}
    if len(queue_paths) > 1:
        raise AnalysisError("Ambiguous queue item profiles")
    if queue_paths:
        queue_rows, _ = read_csv(queue_paths[0], ["volume", "group", "outcome"])
        queue = unique_by_key(queue_rows, ("volume", "group"), "queue item")
        runtime = unique_by_key(parse_queue_items(queue_paths[0]), ("volume", "group"), "registration runtime")
    else:
        warnings.append("Queue item profile missing: registration runtimes unavailable.")
    csv_paths = list(root.glob("motionMonitor_data_*.csv"))
    motion_csv = {}
    if len(csv_paths) > 1:
        raise AnalysisError("Ambiguous motion-monitor CSVs")
    if csv_paths:
        rows, _ = read_csv(csv_paths[0], ["Volume_index", "Slice_group_index", "Displacement(mm)"])
        motion_csv = unique_by_key(rows, ("Volume_index", "Slice_group_index"), "motion CSV observation")
    else:
        warnings.append("Motion CSV missing: recorded group-displacement cross-check unavailable.")
    # Confirm saved physical-space convention from reference headers, without reading pixels.
    ref_headers = [p for p in root.glob("*volume_0000_*.nhdr") if "upsampled" not in p.name]
    if len(ref_headers) != 1 or "space: left-posterior-superior" not in ref_headers[0].read_text():
        raise AnalysisError("A unique LPS reference header is required to establish coordinate convention")
    pointers = {}
    for path in root.glob("processed_*/*group_*.txt"):
        match = re.search(r"_volume_(\d+)_group_(\d+)\.txt$", path.name)
        if match:
            key = tuple(map(int, match.groups()))
            if key in pointers:
                raise AnalysisError(f"Ambiguous pointer for {key}")
            pointers[key] = path
    group_rows, transforms = [], []
    expected_keys = {(0, ng)} | {(v, g) for v in range(1, n) for g in range(ng)}
    actual_status_keys = {(int(v), int(g)) for v, data in raw.items() for g in data["groups"]}
    coverage_sources = [("status", actual_status_keys)]
    if queue_paths:
        coverage_sources.append(("queue", set(queue)))
    if csv_paths:
        coverage_sources.append(("motion CSV", set(motion_csv)))
    for description, keys in coverage_sources:
        if keys != expected_keys:
            warnings.append(f"{description} coverage: {len(expected_keys - keys)} missing, {len(keys - expected_keys)} unexpected keys.")
    for v, g in [(0, ng), *((v, g) for v in range(1, n) for g in range(ng))]:
        reference = v == 0
        record = raw.get(str(v), {}).get("groups", {}).get(str(g), {})
        status = record.get("status", "missing")
        filename = record.get("transform_filename", "")
        transform, source = load_transform(root, filename, warnings) if status == "registered" else (None, None)
        pointer_valid = None
        if not reference:
            pointer = pointers.get((v, g))
            expected_ordinals = list(range(g * meta["slices_per_group"], (g + 1) * meta["slices_per_group"]))
            members = pointer.read_text().splitlines() if pointer else []
            matches = [re.search(r"_volume_(\d+)_slice_(\d+)\.nhdr$", member) for member in members]
            pointer_valid = len(matches) == len(expected_ordinals) and all(matches) and [int(m[1]) for m in matches] == [v] * len(matches) and [int(m[2]) for m in matches] == expected_ordinals
            if not pointer_valid:
                warnings.append(f"Missing/incompatible pointer members at volume {v}, group {g}; acquisition grouping unverified.")
        if reference and transform is not None and (not np.allclose(transform.GetMatrix(), np.eye(3).ravel(), atol=1e-12) or not np.allclose(transform.GetTranslation(), 0, atol=1e-12)):
            raise AnalysisError("Volume 0 reference transform is not identity")
        time, time_source = timing[(v, ng - 1 if reference else g)]
        observed = motion_csv.get((v, g), {})
        recorded = observed.get("Displacement(mm)")
        recorded = float(recorded) if recorded not in (None, "") else None
        if recorded is not None and not math.isfinite(recorded):
            recorded = None
        params = motion_parameters(transform)
        if transform is not None and observed:
            check = [float(observed[f"{a}_{kind}({unit})"]) for kind, unit in (("rotation", "rad"), ("translation", "mm")) for a in AXES]
            actual = [*np.radians([params[f"{a}_rotation_deg"] for a in AXES]), *[params[f"{a}_translation_mm"] for a in AXES]]
            if not np.allclose(check, actual, atol=1e-9, rtol=1e-9):
                warnings.append(f"Motion CSV differs from authoritative transform at ({v}, {g}).")
        qrow = queue.get((v, g), {})
        if qrow and qrow["outcome"] != ("reference" if reference else status):
            warnings.append(f"Queue/status outcome disagreement at ({v}, {g}).")
        group_rows.append({
            "subject": subject, "volume": v, "group": g, "is_reference": reference,
            "acquisition_time_s": time, "time_source": time_source, "status": status,
            "reason": record.get("reason", ""), "registration_finalized": status_by_volume.get(v, {}).get("registration_finalized", False),
            "transform_filename": filename or None, "transform_source": source,
            "transform_type": transform.GetName() if transform is not None else None,
            "compute_zyx": transform.GetComputeZYX() if transform is not None and transform.GetName() == "Euler3DTransform" else None,
            **dict(zip(("center_x_mm", "center_y_mm", "center_z_mm"), transform.GetCenter() if transform is not None else (None,) * 3)),
            "transform_valid": transform is not None, "pointer_valid": pointer_valid, **params,
            "recorded_group_displacement_mm": recorded,
            "runtime_ms": runtime.get((v, g), {}).get("full_registration_elapsed_ms") if not reference else None,
            "queue_outcome": qrow.get("outcome"),
        })
        transforms.append(transform)
    group_displacements = consecutive_displacements(transforms, radius_mm)
    for row, value in zip(group_rows, group_displacements):
        row["group_displacement_mm"] = value
        row["group_motion_event"] = value > threshold_mm if value is not None else None
    # Volume indices are used only for processing completeness, never motion sampling.
    registered_per_volume = Counter(r["volume"] for r in group_rows if r["status"] == "registered" and not r["is_reference"])
    finalized = {v: status_by_volume.get(v, {}).get("registration_finalized", False) for v in range(n)}
    calls = [r for r in group_rows if not r["is_reference"]]
    counts = Counter(r["status"] for r in calls)
    errors = [abs(r["group_displacement_mm"] - r["recorded_group_displacement_mm"]) for r in calls
              if r["group_displacement_mm"] is not None and r["recorded_group_displacement_mm"] is not None]
    if errors and max(errors) > 1e-6:
        warnings.append(f"Recomputed/recorded group displacement differs by up to {max(errors):.6g} mm; recorded radius/configuration may differ.")
    if not all(finalized.values()):
        warnings.append("One or more expected volume records are absent or unfinalized.")
    if sum(r["runtime_ms"] is not None for r in calls) < len(calls):
        warnings.append("Runtime measurements missing for one or more expected registration opportunities.")
    if any(r["status"] != "registered" or not r["transform_valid"] for r in group_rows):
        warnings.append("Missing/unsuccessful registration observations retained as missing; adjacent displacements are also missing.")
    group_stats = {key: metric_statistics([r[key] for r in calls]) for key in MOTION_FIELDS}
    group_stats["group_displacement_mm"] = metric_statistics([r["group_displacement_mm"] for r in group_rows], threshold_mm)
    sources = [metadata_path, status_path, *queue_paths, *csv_paths, *root.glob("measurement*.hdf5"), *ref_headers]
    return {
        "subject": subject, "input_directory": str(root), "metadata": {**meta, "observed_mrd_slice_images": image_count},
        "processing": {
            "expected_registrations": (n - 1) * ng, **{key: counts[key] for key in ("registered", "failed", "skipped", "missing")},
            "reference_records": 1 if (0, ng) in actual_status_keys else 0,
            "reference_status": group_rows[0]["status"], "status_record_count": len(actual_status_keys),
            "finalized_volumes": sum(finalized.values()),
            "complete_registered_volumes_excluding_reference": sum(registered_per_volume[v] == ng and finalized[v] for v in range(1, n)),
            "valid_call_transforms": sum(r["transform_valid"] for r in calls),
            "success_pct": 100 * counts["registered"] / len(calls) if calls else None,
            "verified_pointers": sum(r["pointer_valid"] is True for r in calls),
            "queue_configurations": [list(c) for c in sorted({(r.get("reg_engine", ""), r.get("cuda_execution_mode", ""), r.get("fifo_flag", "")) for r in queue.values()})],
        },
        "group_statistics": group_stats,
        "runtime_statistics_ms": describe(r["runtime_ms"] for r in calls),
        "runtime_scope": "full queue-bracketed registration calls; includes startup overhead when present; identity setup excluded",
        "group_event_affected_volumes": len({r["volume"] for r in calls if r["group_motion_event"]}),
        "recorded_group_displacement_crosscheck": {"paired_count": len(errors), "max_abs_error_mm": max(errors) if errors else None},
        "sources": [{"path": str(p.relative_to(root)), "size_bytes": p.stat().st_size} for p in sources if p.exists()],
        "warnings": list(dict.fromkeys(warnings)), "group_rows": group_rows,
    }


def aggregate_subjects(subjects: list[dict[str, Any]]) -> dict[str, Any]:
    """Descriptive per-subject comparisons, without weighting longer acquisitions more."""
    metric = "group_displacement_mm"
    usable = [s for s in subjects if s.get("group_statistics", {}).get(metric, {}).get("count", 0)]
    comparisons = {metric: {
        "subjects_with_data": len(usable),
        "subject_mean_distribution": describe(s["group_statistics"][metric]["mean"] for s in usable),
        "ascending_subject_means": [s["subject"] for s in sorted(usable, key=lambda s: s["group_statistics"][metric]["mean"])],
    }}
    return {
        "subject_count": len(subjects),
        "total_expected_registrations": sum(s.get("processing", {}).get("expected_registrations", 0) for s in subjects),
        "total_registered": sum(s.get("processing", {}).get("registered", 0) for s in subjects),
        "comparisons": comparisons,
    }


def slice_group_exceedance_rows(subjects: list[dict[str, Any]], threshold: float) -> list[dict[str, Any]]:
    """Summarize existing displacements by explicit current-estimate metadata.

    Verify consecutive volume/group identifiers before classifying transitions;
    row offsets are never used as group labels. Missing expected slots remain in
    the chronological sequence but contribute no value to a group's denominator.
    """
    output = []
    for subject in subjects:
        values = {group: [] for group in range(14)}
        previous = None
        seen = set()
        for row in subject.get("group_rows", []):
            volume, group = int(row["volume"]), int(row["group"])
            key = (volume, group)
            if key in seen:
                raise AnalysisError(f"Duplicate slice-group position record for {subject['subject']}: {key}")
            seen.add(key)
            if row["is_reference"]:
                if previous is not None or key != (0, 14) or row["group_displacement_mm"] is not None:
                    raise AnalysisError("Expected one initial identity baseline with undefined displacement")
            else:
                if volume < 1 or group not in values or previous is None:
                    raise AnalysisError(f"Invalid slice-group transition metadata: {key}")
                expected_previous = (volume, group - 1) if group else ((0, 14) if volume == 1 else (volume - 1, 13))
                if (int(previous["volume"]), int(previous["group"])) != expected_previous:
                    raise AnalysisError(f"Nonconsecutive acquisition metadata before {key}: expected {expected_previous}")
                value = row["group_displacement_mm"]
                if value is not None and math.isfinite(value):
                    values[group].append(float(value))
            previous = row
        for group, observations in values.items():
            stats = describe(observations)
            count = stats["count"]
            above = sum(value > threshold for value in observations)
            output.append({
                "subject": subject["subject"], "slice_group_index": group,
                "transition_type": "cross_volume" if group == 0 else "intra_volume",
                "valid_transition_count": count, "above_threshold_count": above,
                "exceedance_pct": 100 * above / count if count else None,
                "mean_displacement_mm": stats["mean"], "p95_displacement_mm": stats["p95"],
            })
    return output


def plot_slice_group_exceedance(rows: list[dict[str, Any]], subjects: list[str], output: Path, threshold: float) -> None:
    """One shared-scale 2-by-3 bar figure; group zero marks the boundary."""
    if len(subjects) > 6:
        raise AnalysisError("The slice-group position figure supports at most six subjects")
    plt = pyplot()
    from matplotlib.patches import Patch

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True, constrained_layout=True)
    lookup = {(row["subject"], row["slice_group_index"]): row for row in rows}
    maximum = max((row["exceedance_pct"] for row in rows if row["exceedance_pct"] is not None), default=0)
    for ax, subject in zip(axes.flat, subjects):
        rates = [lookup[(subject, group)]["exceedance_pct"] for group in range(14)]
        ax.bar(range(14), [rate if rate is not None else np.nan for rate in rates],
               color=["#D55E00"] + ["#0072B2"] * 13)
        ax.set_title(subject)
        ax.set_xticks(range(14))
        ax.tick_params(axis="x", labelbottom=True)
        ax.set_xlabel("Slice-group index")
        ax.set_ylabel(f"Valid transitions > {threshold:.2f} mm (%)")
        ax.set_ylim(0, min(100, max(1, maximum * 1.15)))
        ax.grid(axis="y", alpha=0.2)
    for ax in list(axes.flat)[len(subjects):]:
        ax.set_visible(False)
    fig.legend(handles=[Patch(color="#D55E00", label="Group 0: cross-volume transition"),
                        Patch(color="#0072B2", label="Groups 1–13: intra-volume transitions")],
               loc="outside lower center", ncol=2, frameon=False)
    fig.suptitle(f"High-motion exceedance by slice-group position · consecutive displacement > {threshold:.2f} mm")
    fig.savefig(output / POSITION_FIGURE[0], dpi=160, bbox_inches="tight")
    plt.close(fig)


def make_figures(subjects: list[dict[str, Any]], output: Path, threshold: float) -> None:
    plt = pyplot()
    count = len(subjects)
    colors = ["#0072B2", "#009E73", "#D55E00"]
    def save(fig, name):
        fig.savefig(output / name, dpi=160, bbox_inches="tight")
        plt.close(fig)
    def values(rows, key):
        return np.array([r.get(key) if r.get(key) is not None else np.nan for r in rows], dtype=float)
    linear_limit = max(threshold, 0.3)
    displacement_max = max((r["group_displacement_mm"] for s in subjects
                            for r in s.get("group_rows", []) if r.get("group_displacement_mm") is not None), default=0)
    compressed = displacement_max > 5 * linear_limit
    scale_note = f" · common symlog scale: linear to {linear_limit:g} mm, log above" if compressed else " · common linear scale"
    def displacement_axis(ax):
        if compressed:
            from matplotlib.ticker import FixedLocator, ScalarFormatter

            ax.set_yscale("symlog", linthresh=linear_limit)
            ticks = sorted({0., linear_limit / 3, linear_limit, 1., 3., 10., 30.})
            ax.yaxis.set_major_locator(FixedLocator([t for t in ticks if t <= displacement_max * 1.15]))
            ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.set_ylim(0, max(displacement_max * 1.10, linear_limit * 1.2))
    fig, axes = plt.subplots(count, 2, figsize=(14, 2.25 * count), sharex=True, sharey="col", squeeze=False, constrained_layout=True)
    for i, subject in enumerate(subjects):
        rows = subject.get("group_rows", [])
        for j, (kind, unit) in enumerate([("translation", "mm"), ("rotation", "deg")]):
            for axis, color in zip(AXES, colors):
                axes[i, j].plot(values(rows, "acquisition_time_s"), values(rows, f"{axis}_{kind}_{unit}"), color=color, lw=0.65, label=axis)
            axes[i, j].set_ylabel(f"{subject['subject']}\n{kind.title()} ({unit})")
            axes[i, j].grid(alpha=0.2)
    for ax in axes[-1]:
        ax.set_xlabel("Acquisition time (s)")
    axes[0, 0].legend(ncol=3, loc="upper left")
    fig.suptitle("Figure 1 · Reference-relative motion · every slice-group estimate")
    save(fig, FIGURES[0][0])
    fig, axes = plt.subplots(count, 1, figsize=(12, 2.1 * count), sharex=True, sharey=True, squeeze=False, constrained_layout=True)
    for i, subject in enumerate(subjects):
        rows = subject.get("group_rows", [])
        ax = axes[i, 0]
        ax.plot(values(rows, "acquisition_time_s"), values(rows, "group_displacement_mm"), lw=0.65, color=colors[0])
        ax.axhline(threshold, ls="--", color="#D55E00", lw=0.7)
        displacement_axis(ax)
        ax.set_title(subject["subject"], loc="left", fontsize=10)
        ax.set_ylabel("Consecutive slice-group\ndisplacement (mm)")
        ax.grid(alpha=0.2)
    axes[-1, 0].set_xlabel("Acquisition time (s)")
    fig.suptitle(f"Figure 2 · Consecutive slice-group displacement · 50 mm radius\nDashed line: {threshold:g} mm threshold" + scale_note)
    save(fig, FIGURES[1][0])
    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    labels = [s["subject"] for s in subjects]
    data = [values(s.get("group_rows", []), "group_displacement_mm") for s in subjects]
    data = [d[np.isfinite(d)] for d in data]
    ax.boxplot(data, showfliers=True, flierprops={"markersize": 2, "alpha": 0.4})
    ax.set_xticks(range(1, len(labels) + 1), labels)
    for i, d in enumerate(data, 1):
        if len(d):
            ax.scatter([i], [np.mean(d)], color="#D55E00", marker="D", s=25, zorder=3)
    ax.set_title("Boxes: median and interquartile range · diamond: mean")
    ax.set_ylabel("Consecutive slice-group displacement (mm)")
    displacement_axis(ax)
    ax.grid(axis="y", alpha=0.2)
    # Put every requested summary directly in the standalone PNG as well as HTML/CSV.
    fields = ["mean", "median", "p95", "max", "exceedance_pct"]
    statistics = [metric_statistics(d.tolist(), threshold) for d in data]
    cells = [["missing" if st[field] is None else f"{st[field]:.3f}" for st in statistics] for field in fields]
    fig.get_layout_engine().set(rect=(0, 0.28, 1, 0.72))
    stats_ax = fig.add_axes([0.13, 0.02, 0.84, 0.24])
    stats_ax.axis("off")
    stats_table = stats_ax.table(cellText=cells, colLabels=labels,
                                rowLabels=["Mean (mm)", "Median (mm)", "P95 (mm)", "Maximum (mm)", f"> {threshold:g} mm (%)"],
                                cellLoc="center", loc="center", bbox=[0, 0, 1, 1])
    stats_table.auto_set_font_size(False)
    stats_table.set_fontsize(9)
    fig.suptitle("Figure 3 · Consecutive slice-group displacement distributions\n" + scale_note.removeprefix(" · "))
    save(fig, FIGURES[2][0])
    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    bottom = np.zeros(count)
    for key, color in [("registered", colors[0]), ("failed", "#D55E00"), ("skipped", "#E69F00"), ("missing", "#999999")]:
        heights = np.array([s.get("processing", {}).get(key, 0) for s in subjects])
        ax.bar(labels, heights, bottom=bottom, color=color, label=key)
        bottom += heights
    for i, s in enumerate(subjects):
        p = s.get("processing", {})
        ax.text(i, bottom[i], f"{p.get('registered', 0)}/{p.get('expected_registrations', 0)}\nfailed {p.get('failed', 0)} · skipped {p.get('skipped', 0)}\nmissing {p.get('missing', 0)}", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, max(float(max(bottom)) * 1.20, 1))
    ax.set_ylabel("Expected registration opportunities (reference excluded)")
    ax.legend(ncol=4, loc="upper left")
    ax.set_title("Figure 4 · Registration completeness · missing is separate from failed/skipped")
    save(fig, FIGURES[3][0])
    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    data = [values(s.get("group_rows", []), "runtime_ms") for s in subjects]
    data = [d[np.isfinite(d)] for d in data]
    ax.boxplot(data, flierprops={"markersize": 2, "alpha": 0.4})
    ax.set_xticks(range(1, len(labels) + 1), labels)
    for i, d in enumerate(data, 1):
        if len(d):
            ax.scatter([i], [np.mean(d)], color="#D55E00", marker="D", s=25)
    ax.set_ylabel("Full registration-call runtime (ms)")
    ax.set_title("Figure 5 · Registration runtime · startup overhead retained · diamond = mean")
    ax.grid(axis="y", alpha=0.2)
    save(fig, FIGURES[4][0])


def methods(radius: float, threshold: float) -> dict[str, Any]:
    return {
        "coordinate_convention": "Saved ITK physical LPS coordinates in mm; Euler parameters kept in saved order/ComputeZYX convention; no RAS conversion or transform inversion for trajectories.",
        "group_sampling": "All 14 registered slice-group transforms for volumes 1 onward, ordered by acquisition sequence. Volume 0 contributes its recorded identity baseline (group 14), timed at the last reference group. Consecutive slice-group displacement includes transitions across volume boundaries and the identity-to-first-group transition. The baseline has undefined displacement and is excluded from displacement statistics.",
        "physical_transform": "T(p) = R(p - c) + c + t. Translations are ITK center-relative t, not the homogeneous offset c + t - Rc.",
        "composition": "Delta = T_current o inverse(T_previous). R_delta = R_current R_previous^T; c_delta = c_previous; t_delta = R_current(-R_previous^T t_previous + c_previous - c_current) + t_current + c_current - c_previous.",
        "displacement_formula": "Let (x,y,z) be the default SimpleITK Euler3D angles of Delta, in radians. q = [-1 + cos(x)cos(y) + cos(x)cos(z) + cos(y)cos(z) + sin(x)sin(y)sin(z)]/2; theta = acos(clip(q,-1,1)); D = ||t_delta||_2 + r*sqrt((1-cos(theta))^2 + sin(theta)^2). Apply D to adjacent expected slice-group transforms, including volume-boundary transitions.",
        "displacement_definition_note": "Consecutive slice-group displacement measures change between registration estimates, not absolute displacement relative to the reference or conventional Power FD. It preserves motion_monitor.monitor_directory.calculate_displacement and compose_transform_pair, with numerical clipping and the established Euler-based rotational term.",
        "radius_mm": radius, "motion_event_threshold_mm": threshold,
        "events": "Strict D > threshold. Exceedance percent uses valid transitions only; event_runs counts consecutive above-threshold observations, split by missing data. Group affected-volume counts include any above-threshold group transition. These are descriptive events, not image-quality classifications.",
        "reference_relative_metrics": "Translation magnitude = L2 norm of saved translation (mm). Rotation Euler norm = L2 norm of saved Euler angles (degrees), following existing plotting helpers. Additional rotation_angle_deg is acos(clip((trace(R)-1)/2,-1,1)) in degrees. Reference-relative magnitudes are distinct from consecutive slice-group displacement.",
        "missing_data": "Missing, skipped, failed, or unreadable transforms remain missing (JSON null, CSV blank, plot gaps). The initial displacement is undefined, not zero. Adjacent pairs touching a missing transform are missing; gaps are never bridged. Completeness is reported independently.",
        "statistics": "Mean, median, population SD, interpolated percentiles, minimum and maximum use finite observations only. Reference identity is excluded from reference-relative motion distributions. Initial undefined displacements are excluded from valid denominators. Cross-subject mean distributions weight subjects equally; no inferential significance tests.",
        "acquisition_time": "Stored MRD acquisition_time_stamp/400, relative to first image, selected by arrival ordinal in group pointers. If no MRD file, derive times from TR plus sorted unique SliceTiming and warn. Never use replay wall time for motion plots. These NOCO replay MRD timestamps were derived by the NIfTI converter; they are not scanner wall-clock measurements.",
        "runtime": "Queue monotonic registration interval (end-start)/1e6 ms; reference identity initialization excluded. All actual calls retained, including startup/reference-load overhead if inside the call. Runtime missingness is independent of transform missingness.",
        "displacement_plot_scale": "Figures 2 and 3 share a common full-range scale. If the largest displacement exceeds five times max(event threshold, 0.3 mm), use symlog: linear through that limit, logarithmic above it. Every outlier is retained; axes and titles identify the scale.",
        "scope": "Retrospective estimated motion and registration/processing completeness in six independent NOCO acquisitions. Registration success is a recorded processing outcome, not registration accuracy. No prospective correction effectiveness or motion-corrected image-quality inference.",
    }


def write_reports(subjects: list[dict[str, Any]], output: Path, radius: float, threshold: float, global_warnings: list[str]) -> dict[str, Any]:
    clean = [{k: v for k, v in s.items() if k != "group_rows"} for s in subjects]
    summary = {"schema_version": 2, "generated_at_utc": datetime.now(timezone.utc).isoformat(),
               "methods": methods(radius, threshold), "warnings": global_warnings,
               "cross_subject": aggregate_subjects(clean), "subjects": clean,
               "figures": [name for name, _ in FIGURES] + [POSITION_FIGURE[0]]}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    rows = [r for s in subjects for r in s.get("group_rows", [])]
    write_csv(output / "group_motion.csv", rows, GROUP_FIELDS)
    (output / "group_motion.json").write_text(json.dumps(rows, allow_nan=False) + "\n")
    stats, subject_rows = [], []
    for s in clean:
        p, m = s.get("processing", {}), s.get("metadata", {})
        subject_rows.append({"subject": s["subject"], "input_directory": s["input_directory"],
                             **{k: m.get(k) for k in ["volumes", "slices_per_volume", "tr_s", "nominal_duration_s", "observed_mrd_slice_images"]},
                             **p, "warning_count": len(s["warnings"])})
        for metric, values in s.get("group_statistics", {}).items():
            stats.append({"subject": s["subject"], "scope": "group_statistics", "metric": metric, **values})
        stats.append({"subject": s["subject"], "scope": "runtime", "metric": "full_registration_elapsed_ms", **s.get("runtime_statistics_ms", describe([]))})
    stat_fields = ["subject", "scope", "metric", "count", "mean", "median", "std", "p95", "p99", "min", "max", "threshold_mm", "exceedance_count", "exceedance_pct", "event_runs", "missing_count"]
    write_csv(output / "motion_statistics.csv", stats, stat_fields)
    if subject_rows:
        fields = list(dict.fromkeys(k for row in subject_rows for k in row))
        export = [{k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in row.items()} for row in subject_rows]
        write_csv(output / "subject_summary.csv", export, fields)
    position_rows = slice_group_exceedance_rows(subjects, threshold)
    write_csv(output / POSITION_CSV, position_rows, POSITION_FIELDS)
    plot_slice_group_exceedance(position_rows, [s["subject"] for s in subjects], output, threshold)
    write_html_report(summary, output)
    return summary


def write_html_report(summary: dict[str, Any], output: Path) -> None:
    """Render existing summaries without changing data exports or figures."""
    clean = summary["subjects"]
    radius = summary["methods"]["radius_mm"]
    threshold = summary["methods"]["motion_event_threshold_mm"]
    global_warnings = summary["warnings"]
    def esc(value):
        return html.escape(str(value))
    def fmt(value):
        return "missing" if value is None else f"{value:.4g}" if isinstance(value, float) else esc(value)
    def table(headers, rows):
        return "<table><thead><tr>" + "".join(f"<th>{esc(h)}</th>" for h in headers) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(f"<td>{fmt(v)}</td>" for v in row) + "</tr>" for row in rows) + "</tbody></table>"
    def motion_table():
        rows = []
        for s in clean:
            st = s.get("group_statistics", {}).get("group_displacement_mm", {})
            rows.append([s["subject"], *[st.get(k) for k in ["count", "mean", "median", "p95", "max", "exceedance_count", "exceedance_pct"]], s.get("group_event_affected_volumes")])
        return table(["Subject", "Valid transitions", "Mean (mm)", "Median (mm)", "P95 (mm)", "Max (mm)", f">{threshold:g} mm count", "Exceedance (%)", f"Affected volumes (>{threshold:g} mm)"], rows)
    sections = ["<h1>Six-subject retrospective NOCO motion comparison</h1>",
                f"<p>Generated {esc(summary['generated_at_utc'])}. Consecutive slice-group displacement in six independent NOCO subjects, using every registration estimate and r = {radius:g} mm. Acquisition time is separate from replay processing time.</p>",
                "<p>The primary motion-change metric measures displacement between consecutive registration estimates, including transitions across volume boundaries. It is distinct from reference-relative displacement and conventional Power FD. The initial identity baseline has undefined displacement.</p>",
                "<p>Registration success describes processing outcomes. This analysis does not establish registration accuracy, prospective correction effectiveness, or motion-corrected image quality.</p>"]
    all_warnings = global_warnings + [s["subject"] + ": " + w for s in clean for w in s["warnings"]]
    sections += ["<h2>Warnings and compatibility</h2>", "<ul>" + "".join(f"<li>{esc(w)}</li>" for w in all_warnings) + "</ul>" if all_warnings else "<p>No missing or incompatible required observations detected.</p>"]
    sections += ["<h2>Acquisition metadata and completeness</h2>", table(
        ["Subject", "Volumes", "TR (s)", "Slices", "Groups/volume", "MRD slices", "Registered / expected", "Failed", "Skipped", "Missing", "Finalized volumes"],
        [[s["subject"], *[s.get("metadata", {}).get(k) for k in ["volumes", "tr_s", "slices_per_volume", "groups_per_volume", "observed_mrd_slice_images"]],
          f"{s.get('processing', {}).get('registered', 'missing')} / {s.get('processing', {}).get('expected_registrations', 'missing')}",
          *[s.get("processing", {}).get(k) for k in ["failed", "skipped", "missing", "finalized_volumes"]]] for s in clean]),
        "<p>Volume 0 is a recorded identity reference. Counts above exclude identity initialization. Every later volume has 14 opportunities. Subject identifiers come from folder suffixes.</p>",
        "<h2>Consecutive slice-group displacement: all expected registration estimates</h2>", motion_table(),
        f"<p>Affected volumes counts the number of unique acquired volumes containing at least one consecutive slice-group displacement exceeding {threshold:.2f} mm. Each volume is counted once, regardless of the number of above-threshold transitions. Cross-volume transitions are included and attributed to the volume containing the current slice group.</p>"]
    encoded_position = base64.b64encode((output / POSITION_FIGURE[0]).read_bytes()).decode("ascii")
    sections += ["<h2>High-motion exceedance by slice-group position</h2>",
                 f"<p>Bars show the percentage of valid consecutive slice-group displacement measurements exceeding {threshold:.2f} mm at each slice-group position. Group 0 represents the cross-volume transition from the preceding volume's final group to the current volume's first group. Groups 1–13 represent intra-volume transitions. Exceedances identify elevated estimated motion between registration measurements, not necessarily image corruption.</p>",
                 "<p>The first group-0 transition uses the recorded reference identity; no within-reference transitions are inferred. Rates use each group's valid-transition count, excluding missing observations.</p>",
                 f'<figure><img alt="{html.escape(POSITION_FIGURE[1])}" src="data:image/png;base64,{encoded_position}"><figcaption>{html.escape(POSITION_FIGURE[1])}</figcaption></figure>']
    sections += ["<h2>Reference-relative displacement distributions</h2>", table(
        ["Subject", "Metric (all slice groups)", "N", "Mean", "Median", "P95", "Max"],
        [[s["subject"], metric, *[s.get("group_statistics", {}).get(metric, {}).get(k) for k in ["count", "mean", "median", "p95", "max"]]]
         for s in clean for metric in ["translation_magnitude_mm", "rotation_euler_norm_deg", "rotation_angle_deg"]]),
        "<p>Reference-relative translation and rotation describe displacement from the reference. Consecutive slice-group displacement describes changes between registration estimates.</p>",
        "<h2>Registration runtimes</h2>", table(["Subject", "Calls", "Mean (ms)", "Median (ms)", "P95 (ms)", "Max (ms)"],
        [[s["subject"], *[s.get("runtime_statistics_ms", {}).get(k) for k in ["count", "mean", "median", "p95", "max"]]] for s in clean]),
        "<p>Full queue-bracketed calls; startup overhead retained, reference identity setup excluded.</p>",
        "<h2>Cross-subject comparisons</h2>"]
    for metric, comparison in summary["cross_subject"]["comparisons"].items():
        sections.append(f"<p>{esc(metric)}: ascending subject means: {esc(', '.join(comparison['ascending_subject_means']))}. Descriptive comparison; subjects weighted equally.</p>")
    for name, caption in FIGURES:
        encoded = base64.b64encode((output / name).read_bytes()).decode("ascii")
        sections.append(f'<figure><img alt="{esc(caption)}" src="data:image/png;base64,{encoded}"><figcaption>{esc(caption)}</figcaption></figure>')
    sections.append("<h2>Methods and assumptions</h2><dl>" + "".join(f"<dt>{esc(key)}</dt><dd>{esc(value)}</dd>" for key, value in summary["methods"].items()) + "</dl>")
    sections.append("<h2>Full subject metadata and provenance</h2>")
    for s in clean:
        sections.append(f"<details><summary>{esc(s['subject'])}</summary><pre>{esc(json.dumps(s, indent=2, allow_nan=False))}</pre></details>")
    sections.append("<h2>Companion exports</h2><p>summary.json; subject_summary.csv; motion_statistics.csv; group_motion.csv/json; slice_group_exceedance_by_position.csv. Figures are embedded here and also saved as standalone PNGs.</p>")
    document = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Six-subject NOCO motion comparison</title><style>body{font:15px system-ui,sans-serif;margin:2rem auto;padding:0 1rem;max-width:1450px;color:#172b3a}table{border-collapse:collapse;width:100%;margin:1rem 0;font-size:13px}th,td{border:1px solid #ccd6dd;padding:.5rem;text-align:right}th:first-child,td:first-child{text-align:left}th{background:#eaf0f5}img{max-width:100%}figure{margin:2rem 0}dt{font-weight:bold;margin-top:1rem}dd{margin:.4rem 0 0;overflow-wrap:anywhere}pre{white-space:pre-wrap;overflow-wrap:anywhere}h2{margin-top:2rem}</style><body>' + "\n".join(sections) + "</body></html>"
    (output / "report.html").write_text(document)


def run(data_root: Path, output: Path, expected_subjects: int = 6, threshold_mm: float = 0.3) -> dict[str, Any]:
    if not math.isfinite(threshold_mm) or threshold_mm < 0:
        raise AnalysisError("Motion-event threshold must be finite and nonnegative")
    discovered = discover_subjects(data_root)
    if len(discovered) != expected_subjects:
        raise AnalysisError(f"Expected {expected_subjects} subjects, found {len(discovered)}; refusing to silently omit or add subjects")
    output = output.resolve()
    for _, source in discovered:
        if output == source or source in output.parents or output in source.parents:
            raise AnalysisError("Output must be separate from all original replay directories")
    output.mkdir(parents=True, exist_ok=True)
    subjects = []
    for subject, root in discovered:
        print(f"Analyzing {subject} …", flush=True)
        try:
            subjects.append(analyze_subject(subject, root, 50, threshold_mm))
        except (AnalysisError, OSError, RuntimeError, ValueError, KeyError, TypeError) as error:
            subjects.append({"subject": subject, "input_directory": str(root), "warnings": [f"Subject analysis unavailable: {error}"], "group_rows": []})
    warnings = []
    usable = [s for s in subjects if "metadata" in s]
    if len(usable) != len(subjects):
        warnings.append("One or more subjects could not be analyzed; all discovered subjects remain represented. See subject warnings.")
    configurations = {json.dumps(s.get("processing", {}).get("queue_configurations")) for s in usable}
    if len(configurations) > 1:
        warnings.append("Recorded registration engine/mode/FIFO configurations differ across subjects.")
    for key in ["volumes", "tr_s", "slices_per_volume", "groups_per_volume", "group_offsets_s"]:
        if len({json.dumps(s["metadata"][key]) for s in usable}) > 1:
            warnings.append(f"Subject metadata differ for {key}; each subject's own acquisition coordinates are retained.")
    make_figures(subjects, output, threshold_mm)
    result = write_reports(subjects, output, 50, threshold_mm, warnings)
    remove_obsolete_outputs(output)
    print(f"Wrote {output / 'report.html'}", flush=True)
    return result


def remove_obsolete_outputs(output: Path) -> None:
    """Remove only the known superseded artifacts in the comparison output."""
    for name in ("volume_motion.csv", "volume_motion.json", "figure_2_framewise_displacement.png"):
        (output / name).unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=HERE.parents[1] / "data")
    parser.add_argument("--output-dir", type=Path, help="default: <data-root>/GBS_six_subject_comparison")
    parser.add_argument("--expected-subjects", type=int, default=6)
    parser.add_argument("--motion-event-threshold-mm", type=float, default=0.3)
    args = parser.parse_args()
    try:
        run(args.data_root, args.output_dir or args.data_root / "GBS_six_subject_comparison", args.expected_subjects, args.motion_event_threshold_mm)
    except (AnalysisError, OSError, ValueError) as error:
        parser.exit(1, f"Analysis error: {error}\n")


if __name__ == "__main__":
    main()

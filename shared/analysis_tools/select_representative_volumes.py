#!/usr/bin/env python3
"""Select MRI inspection volumes without changing replay data or comparison reports.

Run from the repository root:
  python3 shared/analysis_tools/select_representative_volumes.py \
      --image-root /path/to/source/scan_data

Only representative_volume_selection.{csv,json,md} are written under the new
representative_volumes directory. No image extraction or image writing occurs.
"""
from __future__ import annotations

import argparse
import copy
from collections import Counter
from itertools import combinations
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

import h5py
import nibabel as nib
import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))
import compare_subject_motion as motion  # noqa: E402
from analysis_common import AnalysisError, describe, discover_unique, read_json, write_csv  # noqa: E402
from shared.myhelpers_jauger.nifti_to_mrd import read_image  # noqa: E402

SUBJECTS = ("GBS024", "GBS029", "GBS061", "GBS067", "GBS069", "GBS071")
THRESHOLD_MM = 0.30
RADIUS_MM = 50.0
CATEGORIES = ("minimal", "low_moderate", "high")
GLOBAL_MAXIMUM_CATEGORY = "Global maximum displacement"
INTRA_MAXIMUM_CATEGORY = "Strictly intra-volume maximum displacement"
MAXIMUM_COUNT_CATEGORY = "Maximum slice-group exceedance count"
CSV_FIELDS = [
    "subject", "volume_index", "category", "mean_intra_volume_displacement_mm",
    "maximum_intra_volume_displacement_mm", "p95_intra_volume_displacement_mm",
    "exceedance_count", "exceedance_group_indices", "cross_volume_displacement_mm",
    "maximum_reference_relative_translation_mm", "maximum_reference_relative_rotation_deg",
    "selection_rationale", "source_directory", "nifti_path", "nifti_volume_index",
    "mrd_path", "mrd_image_number_first", "mrd_image_number_last", "nrrd_directory",
    "maximum_transition_slice_group_index", "maximum_transition_displacement_mm", "maximum_transition_type",
    "valid_all_group_transition_count", "all_group_exceedance_count", "all_group_exceedance_pct",
    "all_group_exceedance_indices", "mean_all_group_displacement_mm", "p95_all_group_displacement_mm",
    "maximum_all_group_displacement_mm", "group_0_exceedance_count", "intra_volume_exceedance_count",
]


def category(maximum: float | None) -> str:
    if maximum is None:
        return "incomplete"
    if maximum <= .10:
        return "minimal"
    if .15 <= maximum <= .30:
        return "low_moderate"
    if maximum > .30:
        return "high"
    return "unclassified"


def volume_characteristics(subject: dict) -> list[dict]:
    grouped: dict[int, dict[int, dict]] = {}
    for row in subject["group_rows"]:
        if row["is_reference"]:
            continue
        volume, group = int(row["volume"]), int(row["group"])
        groups = grouped.setdefault(volume, {})
        if group in groups:
            raise AnalysisError(f"Duplicate volume/group observation: {volume}/{group}")
        groups[group] = row
    output = []
    for volume in range(1, subject["metadata"]["volumes"]):
        groups = grouped.get(volume, {})
        if set(groups) != set(range(14)):
            raise AnalysisError(f"Expected explicit groups 0–13 for {subject['subject']} volume {volume}")
        rows = [groups[g] for g in range(14)]
        intra = [groups[g]["group_displacement_mm"] for g in range(1, 14)]
        stats = describe(intra)
        complete = stats["count"] == 13 and all(
            r["status"] == "registered" and r["transform_valid"] and r["registration_finalized"] and r["pointer_valid"]
            for r in rows)
        output.append({
            "subject": subject["subject"], "volume_index": volume,
            "category": category(stats["max"]) if complete else "incomplete",
            "complete_registration_records": complete, "valid_intra_volume_transitions": stats["count"],
            "mean_intra_volume_displacement_mm": stats["mean"],
            "maximum_intra_volume_displacement_mm": stats["max"],
            "p95_intra_volume_displacement_mm": stats["p95"],
            "exceedance_count": sum(value is not None and value > THRESHOLD_MM for value in intra),
            "exceedance_group_indices": [g for g in range(1, 14) if intra[g - 1] is not None and intra[g - 1] > THRESHOLD_MM],
            "moderately_elevated_transition_count": sum(value is not None and value >= .15 for value in intra),
            "cross_volume_displacement_mm": groups[0]["group_displacement_mm"],
            "maximum_reference_relative_translation_mm": describe(r["translation_magnitude_mm"] for r in rows)["max"],
            "maximum_reference_relative_rotation_deg": describe(r["rotation_angle_deg"] for r in rows)["max"],
            "source_directory": subject["input_directory"],
            "_trajectory": rows,
        })
    lookup = {r["volume_index"]: r for r in output}
    for record in output:
        record["adjacent_large_event_volumes"] = [
            v for v in [record["volume_index"] - 1, record["volume_index"] + 1]
            if v in lookup and (lookup[v]["exceedance_count"] or (lookup[v]["cross_volume_displacement_mm"] or 0) > THRESHOLD_MM)
        ]
    return output


def rank(record: dict, requested: str) -> tuple:
    boundary_large = (record["cross_volume_displacement_mm"] or 0) > THRESHOLD_MM
    if requested == "minimal":
        values = (bool(record["adjacent_large_event_volumes"]), boundary_large,
                  record["mean_intra_volume_displacement_mm"], record["maximum_intra_volume_displacement_mm"])
    elif requested == "low_moderate":
        values = (-record["moderately_elevated_transition_count"], boundary_large,
                  -record["mean_intra_volume_displacement_mm"], -record["p95_intra_volume_displacement_mm"])
    else:
        values = (-record["maximum_intra_volume_displacement_mm"], -record["p95_intra_volume_displacement_mm"],
                  -record["exceedance_count"], -record["mean_intra_volume_displacement_mm"])
    return (*values, record["subject"], record["volume_index"])


def select_primary(records: list[dict]) -> list[dict]:
    """Prefer six distinct subjects; preserve the supplied category boundaries."""
    best = {cat: {} for cat in CATEGORIES}
    for record in records:
        cat = record["category"]
        if cat in best:
            previous = best[cat].get(record["subject"])
            if previous is None or rank(record, cat) < rank(previous, cat):
                best[cat][record["subject"]] = record
    choices = []
    # Lexicographic priority: clear high-amplitude examples first, then quiet
    # minimal examples, then multi-transition moderate examples, with stable ties.
    for high_subjects in combinations(sorted(best["high"]), 2):
        if "GBS071" in best["high"] and "GBS071" not in high_subjects:
            continue
        for minimal_subjects in combinations(sorted(set(best["minimal"]) - set(high_subjects)), 2):
            remaining = set(SUBJECTS) - set(high_subjects) - set(minimal_subjects)
            if len(remaining) != 2 or not remaining <= set(best["low_moderate"]):
                continue
            chosen = {"high": [best["high"][s] for s in high_subjects],
                      "minimal": [best["minimal"][s] for s in minimal_subjects],
                      "low_moderate": [best["low_moderate"][s] for s in sorted(remaining)]}
            key = tuple(tuple(sorted(rank(r, cat) for r in chosen[cat])) for cat in ("high", "minimal", "low_moderate"))
            choices.append((key, chosen))
    if not choices:
        raise AnalysisError("Cannot obtain two examples per category covering six distinct subjects; review eligibility before changing selection policy")
    chosen = min(choices, key=lambda item: item[0])[1]
    return [r for cat in CATEGORIES for r in sorted(chosen[cat], key=lambda row: rank(row, cat))]


def select_diagnostics(records: list[dict], primary: list[dict]) -> list[dict]:
    excluded = {(r["subject"], r["volume_index"]) for r in primary}
    selected = []
    for subject in ("GBS061", "GBS067", "GBS071"):
        candidates = [r for r in records if r["subject"] == subject and r["complete_registration_records"]
                      and {3, 4} & set(r["exceedance_group_indices"])
                      and (subject, r["volume_index"]) not in excluded]
        def diagnostic_rank(r):
            d = {row["group"]: row["group_displacement_mm"] for row in r["_trajectory"]}
            return (not {3, 4} <= set(r["exceedance_group_indices"]),
                    (r["cross_volume_displacement_mm"] or 0) > THRESHOLD_MM,
                    len(set(r["exceedance_group_indices"]) - {3, 4}),
                    -min(d[3], d[4]), -max(d[3], d[4]), r["volume_index"])
        if candidates:
            selected.append(min(candidates, key=diagnostic_rank))
    return selected


def largest_transition(records: list[dict], *, intra_only: bool = False) -> tuple[dict, dict]:
    """Find a finite consecutive displacement, using current-estimate group IDs."""
    candidates = [(record, row) for record in records for row in record["_trajectory"]
                  if not row.get("is_reference", False) and 0 <= row["group"] <= 13
                  and (not intra_only or row["group"] >= 1)
                  and row["group_displacement_mm"] is not None and math.isfinite(row["group_displacement_mm"])]
    if not candidates:
        raise AnalysisError("No valid consecutive transitions for maximum selection")
    return min(candidates, key=lambda pair: (-pair[1]["group_displacement_mm"], pair[0]["subject"], pair[0]["volume_index"], pair[1]["group"]))


def maximum_selection(record: dict, row: dict, category_label: str, existing: list[dict], nifti: Path | None) -> dict:
    previous = next((r for r in existing if r["subject"] == record["subject"] and r["volume_index"] == record["volume_index"]), None)
    result = copy.deepcopy(previous) if previous is not None else decorate_selection(record, False, nifti)
    result["category"] = category_label
    transition_type = "cross_volume" if row["group"] == 0 else "intra_volume"
    result.update(maximum_transition_slice_group_index=row["group"], maximum_transition_displacement_mm=row["group_displacement_mm"],
                  maximum_transition_type=transition_type)
    result["selection_rationale"] = (
        f"{category_label}: {row['group_displacement_mm']:.9f} mm at current slice group {row['group']} ({transition_type}). "
        + (f"The transition ends in volume {record['volume_index']} and starts at the preceding volume's final registered group; "
           "it is excluded from all intra-volume mean/P95/maximum and exceedance counts. " if row["group"] == 0 else
           "This maximum is restricted to groups 1–13 and excludes every group-0 boundary transition. ")
        + (f"This is the same acquired volume as the existing {previous['category']} selection, repeated only to identify its distinct maximum-transition role."
           if previous is not None else "The original representative and diagnostic selections are retained unchanged."))
    if previous is not None:
        result["existing_selection_reference"] = {"subject": previous["subject"], "volume_index": previous["volume_index"], "category": previous["category"]}
    return result


def all_group_exceedance_statistics(record: dict) -> dict:
    """Count groups 0–13, keeping group 0's cross-volume contribution explicit."""
    rows = record["_trajectory"]
    if len(rows) != 14 or {r["group"] for r in rows} != set(range(14)):
        raise AnalysisError("Exceedance-count selection requires explicit, unique groups 0–13")
    values = [r["group_displacement_mm"] for r in rows]
    stats = describe(values)
    groups = sorted(r["group"] for r in rows if r["group_displacement_mm"] is not None
                    and math.isfinite(r["group_displacement_mm"]) and r["group_displacement_mm"] > THRESHOLD_MM)
    return {
        "subject": record["subject"], "volume_index": record["volume_index"],
        "valid_all_group_transition_count": stats["count"], "all_group_exceedance_count": len(groups),
        "all_group_exceedance_pct": 100 * len(groups) / stats["count"] if stats["count"] else None,
        "all_group_exceedance_indices": groups,
        "mean_all_group_displacement_mm": stats["mean"], "p95_all_group_displacement_mm": stats["p95"],
        "maximum_all_group_displacement_mm": stats["max"],
        "group_0_exceedance_count": int(0 in groups),
        "intra_volume_exceedance_count": sum(group > 0 for group in groups),
    }


def most_exceedances(records: list[dict]) -> tuple[dict, dict, list[dict]]:
    if not records:
        raise AnalysisError("No non-reference volumes for exceedance-count selection")
    ranked = sorted(((record, all_group_exceedance_statistics(record)) for record in records),
                    key=lambda pair: (-pair[1]["all_group_exceedance_count"], pair[0]["subject"], pair[0]["volume_index"]))
    record, stats = ranked[0]
    others = [s for _, s in ranked[1:] if s["all_group_exceedance_count"] == stats["all_group_exceedance_count"]]
    return record, stats, others


def exceedance_count_selection(record: dict, stats: dict, existing: list[dict], nifti: Path | None) -> dict:
    previous = next((r for r in existing if r["subject"] == record["subject"] and r["volume_index"] == record["volume_index"]), None)
    result = copy.deepcopy(previous) if previous is not None else decorate_selection(record, False, nifti)
    result.update(stats)
    result["category"] = MAXIMUM_COUNT_CATEGORY
    result["selection_rationale"] = (
        f"Greatest number of strict >{THRESHOLD_MM:.2f} mm consecutive slice-group displacements across all six subjects: "
        f"{stats['all_group_exceedance_count']}/{stats['valid_all_group_transition_count']} "
        f"({stats['all_group_exceedance_pct']:.6f}%) at groups {stats['all_group_exceedance_indices']}. "
        f"Group 0 contributes {stats['group_0_exceedance_count']} cross-volume exceedance, attributed to the current volume; "
        f"groups 1–13 contribute {stats['intra_volume_exceedance_count']} intra-volume exceedances. "
        "All-group mean/P95/maximum include group 0; the existing intra-volume fields continue to exclude it. "
        "This selects event frequency rather than the largest displacement amplitude.")
    if previous is not None:
        result["existing_selection_reference"] = {"subject": previous["subject"], "volume_index": previous["volume_index"], "category": previous["category"]}
    return result


def net_transform_change(record: dict, start: int, end: int) -> float:
    groups = {r["group"]: r for r in record["_trajectory"]}
    transforms = []
    for group in [start, end]:
        transform, _ = motion.load_transform(Path(record["source_directory"]), groups[group]["transform_filename"], [])
        if transform is None:
            raise AnalysisError("Selected diagnostic transform is unavailable")
        transforms.append(transform)
    return motion.displacement(motion.relative_transform(*transforms), RADIUS_MM)


def find_nifti_sources(image_root: Path | None, subjects: list[dict]) -> dict[str, Path | None]:
    expected = {}
    for subject in subjects:
        meta = read_json(discover_unique(Path(subject["input_directory"]), "*metadata*.json", "metadata"))
        names = [p["value"] for p in meta["userParameters"]["userParameterString"] if p["name"] == "SourceNIfTI"]
        expected[subject["subject"]] = names[0] if len(names) == 1 else None
    matches = {subject: [] for subject in expected}
    if image_root is not None:
        names = set(expected.values()) - {None}
        for directory, _, files in os.walk(image_root):
            for name in names & set(files):
                path = Path(directory) / name
                for subject, expected_name in expected.items():
                    if name == expected_name and any(re.search(rf"(?:^|_){subject}(?:_|$)", part) for part in path.parts):
                        matches[subject].append(path.resolve())
    # Missing or ambiguous source images remain unresolved, never guessed.
    return {subject: paths[0] if len(paths) == 1 else None for subject, paths in matches.items()}


def image_mapping(record: dict, nifti_path: Path | None) -> dict:
    """Verify selected repetition/slice headers, image values, and NRRD ordinals."""
    root = Path(record["source_directory"])
    mrd = discover_unique(root, "measurement*.hdf5", "replay MRD")
    processed = sorted(p for p in root.glob("processed_*") if p.is_dir())
    if len(processed) != 1:
        raise AnalysisError("Ambiguous processed-image directory")
    v = record["volume_index"]
    first, last = v * 56, (v + 1) * 56 - 1
    image = nib.load(nifti_path) if nifti_path is not None else None
    if image is not None and image.shape != (90, 90, 56, 480):
        raise AnalysisError(f"Unexpected source NIfTI dimensions: {image.shape}")
    source = np.asanyarray(image.dataobj[:, :, :, v]) if image is not None else None
    seen, headers, max_error = set(), [], 0.0
    with h5py.File(mrd, "r") as stream:
        for ordinal, number in enumerate(range(first, last + 1)):
            header, _, data = read_image(stream[f"image_{number}"])
            if int(header.repetition) != v or int(header.channels) != 1 or tuple(header.matrix_size) != (90, 90, 1):
                raise AnalysisError(f"MRD header does not map to selected volume {v}: image_{number}")
            spatial_slice = int(header.slice)
            if spatial_slice in seen or not 0 <= spatial_slice < 56:
                raise AnalysisError("Duplicate/out-of-range MRD spatial slice index")
            seen.add(spatial_slice)
            paths = list(processed[0].glob(f"*volume_{v:04d}_slice_{ordinal:04d}.nhdr"))
            if len(paths) != 1:
                raise AnalysisError(f"Selected NRRD slice is missing/ambiguous: {v}/{ordinal}")
            path = paths[0]
            fields = dict(line.split(": ", 1) for line in path.read_text().splitlines() if ": " in line)
            raw_path = path.parent / fields["data file"]
            if raw_path.read_bytes() != data.tobytes():
                raise AnalysisError(f"NRRD/MRD pixel mismatch: {path}")
            if source is not None:
                error = float(np.max(np.abs(data[0, 0].T.astype(float) - source[:, :, spatial_slice].astype(float))))
                max_error = max(max_error, error)
                if error != 0:
                    raise AnalysisError(f"Source NIfTI/MRD voxel mismatch for {v}/{spatial_slice}: {error}")
            headers.append({"arrival_ordinal": ordinal, "mrd_image_key": f"image_{number}", "mrd_repetition": int(header.repetition),
                            "spatial_slice_index": spatial_slice, "nrrd_header": str(path), "raw_path": str(raw_path)})
    return {
        "nifti_path": str(nifti_path) if nifti_path is not None else None,
        "nifti_volume_index": v if nifti_path is not None else None,
        "mrd_path": str(mrd), "mrd_image_number_first": first, "mrd_image_number_last": last,
        "nrrd_directory": str(processed[0]), "verified_slice_count": len(seen),
        "nifti_mrd_maximum_voxel_error": max_error if source is not None else None,
        "image_mapping_verified": True, "slice_mapping": headers,
        "nifti_resolution": "Exact SourceNIfTI basename, subject directory, 4D dimensions, and all selected-volume voxels verified against MRD" if source is not None else "Unresolved: missing or ambiguous source NIfTI; use verified MRD/NRRD slices",
    }


def decorate_selection(record: dict, diagnostic: bool, nifti: Path | None) -> dict:
    selected = {k: v for k, v in record.items() if not k.startswith("_")}
    selected["category"] = "diagnostic" if diagnostic else record["category"]
    selected["motion_category"] = record["category"]
    if diagnostic:
        trajectory = {r["group"]: r for r in record["_trajectory"]}
        net = net_transform_change(record, 2, 4)
        net5 = net_transform_change(record, 2, 5)
        d3, d4 = (trajectory[g]["group_displacement_mm"] for g in [3, 4])
        selected["diagnostic_window"] = {"group_3_displacement_mm": d3, "group_4_displacement_mm": d4,
                                         "net_group_2_to_4_change_mm": net, "net_group_2_to_5_change_mm": net5}
        rationale = (f"Paired group-3/4 estimates: D3={d3:.6f}, D4={d4:.6f} mm; net group-2-to-4 transform change={net:.6f} mm "
                     f"and group-2-to-5={net5:.6f} mm. "
                     + ("The trajectory returns toward its pre-spike estimate, consistent with a transient estimate excursion. " if net < .5 * (d3 + d4) else "The trajectory retains a substantial change after the pair. ")
                     + "Group 0 is separate; this example is not selected for boundary motion. Estimated-transform behavior alone cannot distinguish real movement from registration error.")
    elif record["category"] == "minimal":
        rationale = ("All 13 intra-volume transitions are <=0.10 mm; consistently low mean displacement. "
                     f"Adjacent volumes with >0.30 mm events: {record['adjacent_large_event_volumes']}. "
                     "One of two quiet examples, with subject diversity preserved.")
    elif record["category"] == "low_moderate":
        rationale = (f"Maximum intra-volume displacement is within 0.15–0.30 mm, with {record['moderately_elevated_transition_count']} "
                     "of 13 transitions >=0.15 mm and no intra-volume >0.30 mm event; multiple moderate estimates rather than a single outlier.")
    else:
        rationale = (f"A high-amplitude intra-volume example with {record['exceedance_count']} above-threshold transitions at groups "
                     f"{record['exceedance_group_indices']}; selected by high maximum/P95, not group-0 boundary displacement.")
    selected["selection_rationale"] = rationale
    trajectory_fields = ["group", "acquisition_time_s", "group_displacement_mm", *motion.MOTION_FIELDS,
                         "transform_filename", "transform_source", "center_x_mm", "center_y_mm", "center_z_mm", "compute_zyx"]
    selected["reference_relative_trajectory"] = [{key: row[key] for key in trajectory_fields} for row in record["_trajectory"]]
    selected.update(image_mapping(record, nifti))
    return selected


def write_manifest(document: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    selections = document["selections"]
    rows = [{**r, "exceedance_group_indices": json.dumps(r["exceedance_group_indices"]),
             "all_group_exceedance_indices": json.dumps(r["all_group_exceedance_indices"]) if "all_group_exceedance_indices" in r else None}
            for r in selections]
    write_csv(output / "representative_volume_selection.csv", rows, CSV_FIELDS)
    (output / "representative_volume_selection.json").write_text(json.dumps(document, indent=2, allow_nan=False) + "\n")
    lines = ["# Representative acquired volumes for manual inspection", "",
             "All volume indices are **zero-based** acquisition/MRD repetition indices (reference = 0). "
             "Classification uses only the 13 consecutive transitions at groups **1–13**. Group 0 is listed separately. "
             "This is a deterministic selection for inspection, not a population sample or a claim about image quality.", "",
             "## Eligibility (criteria unchanged)", "",
             "| Subject | Minimal (max ≤0.10) | Low–moderate (0.15≤max≤0.30) | High (max >0.30) | Gap (0.10<max<0.15) |",
             "|---|---:|---:|---:|---:|"]
    for subject, counts in document["eligibility_by_subject"].items():
        lines.append(f"| {subject} | {counts.get('minimal', 0)} | {counts.get('low_moderate', 0)} | {counts.get('high', 0)} | {counts.get('unclassified', 0)} |")
    lines += ["", f"Total eligible: {document['eligibility_totals']}. All 2,874 non-reference volumes have complete records; no threshold adjustments were necessary.", "",
              "## Selection policy", "", document["methods"]["selection_policy"], "",
              "## Primary and separate diagnostic selections", "",
              "Displacements and translations are mm; reference-relative matrix rotation angles are degrees. Exceedance groups exclude group 0.", "",
              "| Category | Subject | Volume | Mean | Maximum | P95 | >0.30 count | Groups | Group 0 | Max ref translation | Max ref angle |",
              "|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|"]
    for r in selections:
        if r["category"] in (GLOBAL_MAXIMUM_CATEGORY, INTRA_MAXIMUM_CATEGORY, MAXIMUM_COUNT_CATEGORY):
            continue
        lines.append(f"| {r['category']} | {r['subject']} | {r['volume_index']} | {r['mean_intra_volume_displacement_mm']:.6f} | {r['maximum_intra_volume_displacement_mm']:.6f} | {r['p95_intra_volume_displacement_mm']:.6f} | {r['exceedance_count']} | {r['exceedance_group_indices']} | {r['cross_volume_displacement_mm']:.6f} | {r['maximum_reference_relative_translation_mm']:.6f} | {r['maximum_reference_relative_rotation_deg']:.6f} |")
    maxima = [r for r in selections if r["category"] in (GLOBAL_MAXIMUM_CATEGORY, INTRA_MAXIMUM_CATEGORY)]
    if maxima:
        lines += ["", "## Global and strictly intra-volume maximum displacement", "",
                  "These are consecutive slice-group displacement measurements with a 50 mm rotational radius, not conventional volume-to-volume FD. "
                  "Each transition belongs to the volume containing its current slice group. Original selections above are unchanged; "
                  "a repeated volume below has an explicitly distinct maximum-transition role.", "",
                  "| Selection | Subject | Zero-based volume | Slice group | Displacement (mm) | Transition type | Intra mean (mm) | Intra P95 (mm) | Intra max (mm) | Intra >0.30 count |",
                  "|---|---|---:|---:|---:|---|---:|---:|---:|---:|"]
        for r in maxima:
            lines.append(f"| {r['category']} | {r['subject']} | {r['volume_index']} | {r['maximum_transition_slice_group_index']} | {r['maximum_transition_displacement_mm']:.9f} | {r['maximum_transition_type']} | {r['mean_intra_volume_displacement_mm']:.6f} | {r['p95_intra_volume_displacement_mm']:.6f} | {r['maximum_intra_volume_displacement_mm']:.6f} | {r['exceedance_count']} |")
        lines += ["", "The global maximum is checked against the maximum of the six subjects' saved comparison statistics. "
                  "A group-0 maximum spans two acquired volumes; its larger boundary value must not be interpreted as intra-volume motion."]
    count_selection = next((r for r in selections if r["category"] == MAXIMUM_COUNT_CATEGORY), None)
    if count_selection is not None:
        r = count_selection
        count_summary = document["maximum_exceedance_count"]
        other = count_summary["gbs071_volume_162_comparison"]
        lines += ["", "## Maximum slice-group exceedance count", "",
                  "This comparison uses **all 14 groups (0–13)** and strict displacement >0.30 mm, with the established 50 mm rotation radius. "
                  "Group 0 is a cross-volume transition assigned to its current volume; groups 1–13 are intra-volume transitions.", "",
                  "| Subject | Zero-based volume | >0.30 count / valid | Percentage | Exceedance groups | All-group mean (mm) | All-group P95 (mm) | All-group max (mm) | Group 0 count | Groups 1–13 count |",
                  "|---|---:|---:|---:|---|---:|---:|---:|---:|---:|"]
        table_rows = [r]
        if (other["subject"], other["volume_index"]) != (r["subject"], r["volume_index"]):
            table_rows.append(other)
        for entry in table_rows:
            lines.append(f"| {entry['subject']} | {entry['volume_index']} | {entry['all_group_exceedance_count']} / {entry['valid_all_group_transition_count']} | {entry['all_group_exceedance_pct']:.6f}% | {entry['all_group_exceedance_indices']} | {entry['mean_all_group_displacement_mm']:.9f} | {entry['p95_all_group_displacement_mm']:.9f} | {entry['maximum_all_group_displacement_mm']:.9f} | {entry['group_0_exceedance_count']} | {entry['intra_volume_exceedance_count']} |")
        lines += ["",
                  "Other volumes tied at the maximum count: " + (json.dumps(count_summary["other_tied_volumes"]) if count_summary["other_tied_volumes"] else "**none** (unique maximum)."), "",
                  f"GBS071 volume 162 has {other['all_group_exceedance_count']}/14 exceedances ({other['all_group_exceedance_pct']:.6f}%), "
                  f"at groups {other['all_group_exceedance_indices']}; group 0 contributes {other['group_0_exceedance_count']}. "
                  + ("It also holds the maximum-count distinction." if (r["subject"], r["volume_index"]) == ("GBS071", 162) else
                     "It retains the largest strictly intra-volume displacement amplitude, but a separate volume has the highest event count. "
                     "Frequency and amplitude identify different examples.")]
    lines += ["", "## Verified image indexing and access", "",
              "Source NIfTI uses `dataobj[:, :, :, volume_index]`; a viewer with one-based frame labels shows frame `volume_index + 1`. "
              "MRD `image_N` numbering was checked against all 56 image headers per selection: `header.repetition` equals the selected volume and "
              "`header.slice` supplies its spatial slice index. NRRD filenames use arrival ordinals, **not** spatial slice indices. "
              "The JSON contains the complete MRD/NRRD slice mapping. All selected NRRD raw pixels matched MRD pixels. "
              "For resolved NIfTI sources, every selected-volume voxel was checked against MRD; unresolved NIfTI paths remain null. "
              "Images have not been copied, extracted, or modified.", "",
              "For a verified 4D NIfTI, inspect `nibabel.load(path).dataobj[..., volume_index]` in memory. "
              "The existing `shared/myhelpers_jauger/nifti_to_mrd.py:read_image` returns `(header, attribute, data)` for each MRD image; "
              "place `data[0,0].T` into spatial slice `header.slice` for the verified repetition. Do not use arrival ordinal as anatomical z."]
    for r in selections:
        nifti_note = (f"- Verified 4D NIfTI: `{r['nifti_path']}`; fourth-axis index **{r['volume_index']}**."
                      if r["nifti_path"] is not None else "- Source NIfTI unresolved; use the verified MRD/NRRD slice mapping below.")
        lines += ["", f"### {r['category']}: {r['subject']} volume {r['volume_index']}", "", r["selection_rationale"], "",
                  f"- Replay source: `{r['source_directory']}`", nifti_note,
                  f"- MRD: `{r['mrd_path']}`; keys `image_{r['mrd_image_number_first']}` through `image_{r['mrd_image_number_last']}` inclusive.",
                  f"- NRRD slices: `{r['nrrd_directory']}`; `*_volume_{r['volume_index']:04d}_slice_0000.nhdr` through `..._slice_0055.nhdr`.", "",
                  "| Group | Consecutive D (mm) | Tx | Ty | Tz | Rx (deg) | Ry (deg) | Rz (deg) | Ref angle (deg) |",
                  "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for row in r["reference_relative_trajectory"]:
            vals = [row["group_displacement_mm"], *[row[f"{a}_translation_mm"] for a in "XYZ"], *[row[f"{a}_rotation_deg"] for a in "XYZ"], row["rotation_angle_deg"]]
            lines.append(f"| {row['group']} | " + " | ".join(f"{value:.6f}" for value in vals) + " |")
    lines += ["", "## Validation and limitations", "", json.dumps(document["validation"], indent=2), "",
              "Group-3/4 paired excursions can reflect registration-position dependence or real movement; these estimates do not establish image corruption. "
              "Large boundary measurements do not change an intra-volume category. Net diagnostic group-2-to-4/5 changes are explicitly separate, nonconsecutive "
              "transform comparisons, used only to assess return toward a pre-spike estimate. The original comparison report and all original replay data remain unchanged.", "",
              "Reproduce with:", "", "```bash", document["reproduce_command"], "```", ""]
    (output / "representative_volume_selection.md").write_text("\n".join(lines))


def run(comparison: Path, output: Path, image_root: Path | None) -> dict:
    comparison, output = comparison.resolve(), output.resolve()
    if output == comparison or not output.is_relative_to(comparison):
        raise AnalysisError("Selection outputs must be in a new subdirectory of the comparison directory")
    summary_path = comparison / "summary.json"
    summary = read_json(summary_path)
    if summary["methods"]["radius_mm"] != RADIUS_MM or summary["methods"]["motion_event_threshold_mm"] != THRESHOLD_MM:
        raise AnalysisError("Existing comparison must use 50 mm radius and 0.30 mm threshold")
    if {s["subject"] for s in summary["subjects"]} != set(SUBJECTS):
        raise AnalysisError("Expected the established six-subject cohort")
    subjects, records, checks = [], [], {}
    for saved in summary["subjects"]:
        root = Path(saved["input_directory"]).resolve()
        if not root.name.endswith("_" + saved["subject"]) or output.is_relative_to(root):
            raise AnalysisError("Subject/source mismatch or unsafe output directory")
        print(f"Rechecking transforms for {saved['subject']} …", flush=True)
        current = motion.analyze_subject(saved["subject"], root, RADIUS_MM, THRESHOLD_MM)
        for key in ("group_statistics", "processing", "group_event_affected_volumes"):
            if current[key] != saved[key]:
                raise AnalysisError(f"Recomputed {key} differs from existing comparison for {saved['subject']}")
        if current["warnings"]:
            raise AnalysisError(f"Unexpected input warnings: {current['warnings']}")
        subjects.append(current)
        subject_records = volume_characteristics(current)
        records.extend(subject_records)
        intra_count = sum(r["valid_intra_volume_transitions"] for r in subject_records)
        boundary_count = sum(r["cross_volume_displacement_mm"] is not None for r in subject_records)
        if intra_count != 479 * 13 or boundary_count != 479:
            raise AnalysisError("Unexpected per-volume registration counts")
        partition = [row["group_displacement_mm"] for record in subject_records for row in record["_trajectory"]]
        observed_stats = current["group_statistics"]["group_displacement_mm"]
        if any(value != observed_stats[key] for key, value in describe(partition).items()):
            raise AnalysisError("Per-volume partition does not reproduce the comparison displacement statistics")
        intra_events = sum(r["exceedance_count"] for r in subject_records)
        boundary_events = sum((r["cross_volume_displacement_mm"] or 0) > THRESHOLD_MM for r in subject_records)
        affected = sum(r["exceedance_count"] > 0 or (r["cross_volume_displacement_mm"] or 0) > THRESHOLD_MM for r in subject_records)
        if intra_events + boundary_events != observed_stats["exceedance_count"] or affected != current["group_event_affected_volumes"]:
            raise AnalysisError("Per-volume event counts do not reproduce the comparison")
        checks[saved["subject"]] = {"group_statistics_exactly_match_comparison": True, "complete_volumes": sum(r["complete_registration_records"] for r in subject_records),
                                    "valid_intra_volume_transitions": intra_count, "valid_boundary_transitions": boundary_count,
                                    "total_transitions": intra_count + boundary_count,
                                    "per_volume_partition_exactly_matches_all_descriptive_statistics": True,
                                    "intra_volume_exceedances": intra_events, "boundary_exceedances": boundary_events,
                                    "affected_volumes": affected}
    primary = select_primary(records)
    diagnostic = select_diagnostics(records, primary)
    sources = find_nifti_sources(image_root, subjects)
    selections = []
    for record in [*primary, *diagnostic]:
        is_diagnostic = (record["subject"], record["volume_index"]) in {(r["subject"], r["volume_index"]) for r in diagnostic}
        print(f"Verifying image mapping: {record['subject']} volume {record['volume_index']} …", flush=True)
        selections.append(decorate_selection(record, is_diagnostic, sources[record["subject"]]))
    unique = {(r["subject"], r["volume_index"]) for r in selections}
    if len(unique) != len(selections) or not all(r["complete_registration_records"] for r in selections):
        raise AnalysisError("Duplicate/incomplete selections")
    original_selections = copy.deepcopy(selections)
    global_record, global_row = largest_transition(records)
    saved_global_maximum = max(s["group_statistics"]["group_displacement_mm"]["max"] for s in summary["subjects"])
    if global_row["group_displacement_mm"] != saved_global_maximum:
        raise AnalysisError("Global maximum differs from the saved six-subject comparison")
    print(f"Verifying global maximum image mapping: {global_record['subject']} volume {global_record['volume_index']} …", flush=True)
    selections.append(maximum_selection(global_record, global_row, GLOBAL_MAXIMUM_CATEGORY, original_selections, sources[global_record["subject"]]))
    intra_maximum = None
    if global_row["group"] == 0:
        intra_record, intra_row = largest_transition(records, intra_only=True)
        intra_maximum = maximum_selection(intra_record, intra_row, INTRA_MAXIMUM_CATEGORY, original_selections, sources[intra_record["subject"]])
        selections.append(intra_maximum)
    prior_roles = copy.deepcopy(selections)
    count_record, count_stats, count_ties = most_exceedances(records)
    if count_stats["valid_all_group_transition_count"] != 14:
        raise AnalysisError("Maximum-count selection does not have all 14 valid transitions")
    print(f"Verifying maximum exceedance-count image mapping: {count_record['subject']} volume {count_record['volume_index']} …", flush=True)
    selections.append(exceedance_count_selection(count_record, count_stats, prior_roles, sources[count_record["subject"]]))
    if selections[:-1] != prior_roles:
        raise AnalysisError("Existing maximum-displacement selection roles changed")
    comparison_162 = all_group_exceedance_statistics(next(r for r in records if r["subject"] == "GBS071" and r["volume_index"] == 162))
    if selections[:len(original_selections)] != original_selections:
        raise AnalysisError("Existing representative selections changed")
    if not all(r["complete_registration_records"] for r in selections):
        raise AnalysisError("Maximum-displacement selection has incomplete registration records")
    # The complete per-volume summary is included in JSON; detailed trajectories
    # for the original selections and maximum-transition additions are retained.
    counts = {s: dict(Counter(r["category"] for r in records if r["subject"] == s)) for s in SUBJECTS}
    command = f"python3 shared/analysis_tools/select_representative_volumes.py --comparison-dir {comparison} --output-dir {output}"
    if image_root is not None:
        command += f" --image-root {image_root}"
    document = {
        "schema_version": 1, "volume_indexing": "zero-based; reference volume is 0; candidates are 1–479",
        "comparison_summary_path": str(summary_path), "comparison_summary_sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        "methods": {"radius_mm": RADIUS_MM, "threshold_mm": THRESHOLD_MM,
                    "intra_volume_groups": list(range(1, 14)), "boundary_group": 0,
                    "displacement": summary["methods"]["displacement_formula"], "composition": summary["methods"]["composition"],
                    "rotation_magnitude": "Matrix axis-angle from rotation_angle_deg, not Euler-vector norm",
                    "selection_policy": "Criteria unchanged: minimal max<=0.10; low_moderate 0.15<=max<=0.30; high max>0.30 mm. "
                    "Require complete 14-group records. Prefer six distinct subjects, two per category, including GBS071 high. "
                    "Enumerate feasible subject assignments and prioritize the best high maxima/P95s, then minimal examples away from adjacent >0.30 mm events "
                    "and with low mean/max, then moderate examples with many >=0.15 mm transitions, quiet boundaries, and high mean/P95. "
                    "Subject and volume index break ties. Diagnostics: one each from GBS061/067/071, favor both groups 3/4 above threshold, quiet boundary, "
                    "few other exceedances, then large paired excursions. Do not duplicate primary selections.",
                    "maximum_displacement_search": "Search every finite consecutive slice-group displacement across all subjects, including group 0; assign the result to its current volume. Verify against the saved comparison maximum. If it is a boundary transition, separately report the largest strictly intra-volume (groups 1–13) transition. Original selections stay unchanged; an already selected volume is repeated only with an explicit maximum-transition role and reference.",
                    "maximum_exceedance_count_search": "Examine every non-reference volume across all subjects. Count strict >0.30 mm transitions across groups 0–13, including the boundary at group 0 assigned to the current volume. Report all tied maximum-count volumes; use subject identifier then zero-based volume index to choose a single deterministic representative. Preserve existing intra-volume fields and store all-group statistics separately.",
                    "missing_data": "Missing transforms/observations are never replaced by zero; incomplete volumes cannot be selected.",
                    "trajectory_source_for_all_volumes": str(comparison / "group_motion.json")},
        "criteria_adjustments": [], "eligibility_by_subject": counts,
        "eligibility_totals": dict(Counter(r["category"] for r in records)),
        "all_non_reference_volume_characteristics": [{k: v for k, v in r.items() if not k.startswith("_")} for r in records],
        "selections": selections,
        "global_maximum": {k: selections[len(original_selections)][k] for k in ["subject", "volume_index", "maximum_transition_slice_group_index", "maximum_transition_displacement_mm", "maximum_transition_type", "nifti_path", "nifti_volume_index"]},
        "strictly_intra_volume_maximum": {k: intra_maximum[k] for k in ["subject", "volume_index", "maximum_transition_slice_group_index", "maximum_transition_displacement_mm", "maximum_transition_type", "nifti_path", "nifti_volume_index"]} if intra_maximum is not None else None,
        "maximum_exceedance_count": {**count_stats, "other_tied_volumes": count_ties,
                                      "nifti_path": selections[-1]["nifti_path"], "nifti_volume_index": selections[-1]["nifti_volume_index"],
                                      "gbs071_volume_162_comparison": comparison_162},
        "validation": {"subjects": checks, "all_primary_subjects_distinct": len({r["subject"] for r in primary}) == 6,
                                                 "original_selections_unchanged": True, "global_maximum_exactly_matches_comparison": True,
                                                 "saved_comparison_global_maximum_mm": saved_global_maximum,
                                                 "existing_maximum_displacement_roles_unchanged": True,
                                                 "exceedance_count_search_volumes": len(records),
                                                 "original_selection_pairs_unique": True,
                                                 "all_selected_subject_volume_pairs_unique": len({(r['subject'], r['volume_index']) for r in selections}) == len(selections),
                                                 "repeated_volume_roles": [r["existing_selection_reference"] for r in selections[len(original_selections):] if "existing_selection_reference" in r],
                                                 "selected_volume_count": len(selections),
                                                 "selected_nifti_voxels_match_mrd": all(r["nifti_mrd_maximum_voxel_error"] == 0 for r in selections),
                                                 "selected_mrd_nrrd_slice_counts": [r["verified_slice_count"] for r in selections]},
        "reproduce_command": command,
    }
    write_manifest(document, output)
    print(f"Wrote selection manifests to {output}", flush=True)
    return document


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-dir", type=Path, default=HERE.parents[1] / "data/GBS_six_subject_comparison")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--image-root", type=Path, help="read-only search root for exact SourceNIfTI basenames")
    args = parser.parse_args()
    try:
        run(args.comparison_dir, args.output_dir or args.comparison_dir / "representative_volumes", args.image_root)
    except (AnalysisError, OSError, ValueError) as error:
        parser.exit(1, f"Selection error: {error}\n")


if __name__ == "__main__":
    main()

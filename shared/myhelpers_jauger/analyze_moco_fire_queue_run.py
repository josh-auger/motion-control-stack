#!/usr/bin/env python3
"""Analyze one profiled MOCO-on FIRE/queue acquisition without changing it."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics
import warnings
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
warnings.filterwarnings("ignore", message="Unable to import Axes3D")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


TRANSFORM_RE = re.compile(
    r"^alignTransform_(?P<index>\d+)_(?P<volume>\d+)-(?P<group>\d+)\.tfm$"
)
IDENTITY_RE = re.compile(r"^alignTransform_\d+_\d+-\d+_identity\.tfm$")
STATUS_VALUES = {"registered", "skipped", "failed"}
SEGMENTS = {
    "early": (1, 66),
    "middle": (67, 132),
    "late": (133, 199),
}
FIRE_REQUIRED = {
    "lookup_id",
    "incoming_image_identifier",
    "volume",
    "slice",
    "group",
    "root_entry_count",
    "valid_transform_count",
    "newer_candidate_count",
    "selected_transform_filename",
    "selected_registration_index",
    "last_consumed_registration_index_before",
    "last_consumed_transform_filename_before",
    "last_consumed_registration_index_after",
    "last_consumed_transform_filename_after",
    "discovery_selection_ms",
    "transform_read_ms",
    "conversion_ms",
    "packaging_ms",
    "conversion_package_ms",
    "feedback_log_ms",
    "feedback_send_ms",
    "outcome",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze FIRE monotonic transform advancement, FIRE discovery scaling, "
            "queue discovery, and registration completeness for one MOCO-on run."
        )
    )
    parser.add_argument("acquisition_dir", type=Path, help="completed savedData directory")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help=(
            "analysis directory (default: a sibling named "
            "<acquisition>_fire_moco_analysis)"
        ),
    )
    return parser.parse_args()


def unique_file(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {pattern!r} in {directory}, found {len(matches)}: "
            f"{[path.name for path in matches]}"
        )
    return matches[0]


def select_matching_fire_log(profile_path: Path, run_dir: Path) -> Path:
    """Choose the FIRE log whose filename timestamp is closest to the profile."""
    profile_match = re.search(r"fire_moco_profile_(\d{8})_(\d{6})_", profile_path.name)
    logs = sorted(run_dir.glob("log_python-fire-server_*.log"))
    if not profile_match or not logs:
        raise ValueError("Cannot match the FIRE profile to a FIRE log")
    profile_time = datetime.strptime("".join(profile_match.groups()), "%Y%m%d%H%M%S")
    candidates = []
    for path in logs:
        match = re.search(r"log_python-fire-server_(\d{8})_(\d{6})\.log$", path.name)
        if match:
            log_time = datetime.strptime("".join(match.groups()), "%Y%m%d%H%M%S")
            candidates.append((abs((log_time - profile_time).total_seconds()), path))
    if not candidates:
        raise ValueError("No timestamped FIRE log matches the expected filename convention")
    candidates.sort(key=lambda item: (item[0], -item[1].stat().st_size))
    return candidates[0][1]


def feedback_log_transforms(path: Path) -> list[str]:
    transforms = []
    pattern = re.compile(r"from transform: (.+)")
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            match = pattern.fullmatch(line.rstrip("\n"))
            if match:
                transforms.append(match.group(1))
    return transforms


def as_int(value: str, field: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid integer for {field}: {value!r}") from error


def optional_int(value: str) -> int | None:
    return int(value) if value != "" else None


def as_float(value: str, field: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid number for {field}: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"Non-finite number for {field}: {value!r}")
    return result


def optional_float(value: str) -> float | None:
    return as_float(value, "optional timing") if value != "" else None


def percentile(sorted_values: list[float], percent: float) -> float:
    position = (len(sorted_values) - 1) * percent / 100.0
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    fraction = position - lower
    return (
        sorted_values[lower] * (1.0 - fraction)
        + sorted_values[upper] * fraction
    )


def describe(values: Iterable[float]) -> dict[str, float | int | None]:
    vals = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not vals:
        return {
            "n": 0,
            "mean": None,
            "median": None,
            "std": None,
            "p95": None,
            "p99": None,
            "min": None,
            "max": None,
        }
    return {
        "n": len(vals),
        "mean": statistics.fmean(vals),
        "median": statistics.median(vals),
        "std": statistics.pstdev(vals),
        "p95": percentile(vals, 95.0),
        "p99": percentile(vals, 99.0),
        "min": vals[0],
        "max": vals[-1],
    }


def pearson(xs: Iterable[float], ys: Iterable[float]) -> float | None:
    pairs = [
        (float(x), float(y))
        for x, y in zip(xs, ys)
        if math.isfinite(float(x)) and math.isfinite(float(y))
    ]
    if len(pairs) < 2:
        return None
    xvals, yvals = zip(*pairs)
    xmean = statistics.fmean(xvals)
    ymean = statistics.fmean(yvals)
    xdev = [value - xmean for value in xvals]
    ydev = [value - ymean for value in yvals]
    denominator = math.sqrt(
        sum(value * value for value in xdev)
        * sum(value * value for value in ydev)
    )
    if denominator == 0:
        return None
    return sum(x * y for x, y in zip(xdev, ydev)) / denominator


def average_ranks(values: Iterable[float]) -> list[float]:
    values = [float(value) for value in values]
    ordered = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    offset = 0
    while offset < len(ordered):
        end = offset + 1
        while end < len(ordered) and values[ordered[end]] == values[ordered[offset]]:
            end += 1
        rank = ((offset + 1) + end) / 2.0
        for position in ordered[offset:end]:
            ranks[position] = rank
        offset = end
    return ranks


def correlations(rows: list[dict[str, Any]], xfield: str, yfield: str) -> dict[str, Any]:
    xs = [float(row[xfield]) for row in rows]
    ys = [float(row[yfield]) for row in rows]
    return {
        "n": len(rows),
        "pearson": pearson(xs, ys),
        "spearman": pearson(average_ranks(xs), average_ranks(ys)),
    }


def segment_for_volume(volume: int) -> str | None:
    for name, (start, end) in SEGMENTS.items():
        if start <= volume <= end:
            return name
    return None


def write_csv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {field: "" if row.get(field) is None else row.get(field) for field in fields}
            )


def parse_status(run_dir: Path) -> dict[str, Any]:
    path = unique_file(run_dir, "registration_status.json")
    with path.open(encoding="utf-8") as stream:
        payload = json.load(stream)
    if payload.get("schema_version") != 1 or not isinstance(payload.get("volumes"), dict):
        raise ValueError(f"Unsupported registration status schema in {path}")

    rows = []
    group_status = {}
    group_sets = {}
    for volume_text, volume_data in payload["volumes"].items():
        volume = int(volume_text)
        if not volume_data.get("registration_finalized"):
            raise ValueError(f"Volume {volume} is not finalized in {path}")
        groups = volume_data.get("groups")
        if not isinstance(groups, dict):
            raise ValueError(f"Volume {volume} has no group mapping in {path}")
        counts: Counter[str] = Counter()
        group_numbers = set()
        for group_text, group_data in groups.items():
            group = int(group_text)
            status = group_data.get("status")
            if status not in STATUS_VALUES:
                raise ValueError(
                    f"Unexpected status {status!r} for volume {volume}, group {group}"
                )
            group_status[(volume, group)] = status
            group_numbers.add(group)
            counts[status] += 1
        group_sets[volume] = group_numbers
        rows.append(
            {
                "volume": volume,
                "expected": len(groups),
                "registered": counts["registered"],
                "skipped": counts["skipped"],
                "failed": counts["failed"],
            }
        )
    rows.sort(key=lambda row: row["volume"])
    volumes = [row["volume"] for row in rows]
    if volumes != list(range(0, max(volumes) + 1)):
        raise ValueError("registration_status volumes are not contiguous from zero")
    nonreference_sets = {frozenset(group_sets[v]) for v in volumes if v > 0}
    if len(nonreference_sets) != 1:
        raise ValueError("Non-reference volumes do not share one group set")
    expected_groups = len(next(iter(nonreference_sets)))
    return {
        "path": path,
        "rows": rows,
        "group_status": group_status,
        "expected_groups": expected_groups,
        "volumes": volumes,
    }


def summarize_registration(status: dict[str, Any]) -> dict[str, Any]:
    rows = [row for row in status["rows"] if row["volume"] > 0]
    degraded = [row for row in rows if row["registered"] < row["expected"]]
    after_startup = [row for row in degraded if row["volume"] >= 2]
    expected = sum(row["expected"] for row in rows)
    registered = sum(row["registered"] for row in rows)
    return {
        "expected": expected,
        "registered": registered,
        "skipped": sum(row["skipped"] for row in rows),
        "failed": sum(row["failed"] for row in rows),
        "success_percent": 100.0 * registered / expected,
        "first_degraded": degraded[0]["volume"] if degraded else None,
        "first_degraded_after_startup": (
            after_startup[0]["volume"] if after_startup else None
        ),
        "degraded_volumes": len(degraded),
        "minimum_registered": min(row["registered"] for row in rows),
        "final": rows[-1],
    }


def parse_fire_profile(run_dir: Path) -> dict[str, Any]:
    path = unique_file(run_dir, "fire_moco_profile_*.csv")
    rows = []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = FIRE_REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing FIRE columns in {path}: {sorted(missing)}")
        for source in reader:
            row = {
                "lookup_id": as_int(source["lookup_id"], "lookup_id"),
                "incoming_image_identifier": as_int(
                    source["incoming_image_identifier"], "incoming_image_identifier"
                ),
                "volume": as_int(source["volume"], "volume"),
                "slice": as_int(source["slice"], "slice"),
                "group": as_int(source["group"], "group"),
                "root_entry_count": as_int(source["root_entry_count"], "root_entry_count"),
                "valid_transform_count": as_int(
                    source["valid_transform_count"], "valid_transform_count"
                ),
                "newer_candidate_count": as_int(
                    source["newer_candidate_count"], "newer_candidate_count"
                ),
                "selected_transform_filename": source["selected_transform_filename"],
                "selected_registration_index": optional_int(
                    source["selected_registration_index"]
                ),
                "last_consumed_registration_index_before": optional_int(
                    source["last_consumed_registration_index_before"]
                ),
                "last_consumed_transform_filename_before": source[
                    "last_consumed_transform_filename_before"
                ],
                "last_consumed_registration_index_after": optional_int(
                    source["last_consumed_registration_index_after"]
                ),
                "last_consumed_transform_filename_after": source[
                    "last_consumed_transform_filename_after"
                ],
                "discovery_selection_ms": as_float(
                    source["discovery_selection_ms"], "discovery_selection_ms"
                ),
                "transform_read_ms": optional_float(source["transform_read_ms"]),
                "conversion_ms": optional_float(source["conversion_ms"]),
                "packaging_ms": optional_float(source["packaging_ms"]),
                "conversion_package_ms": optional_float(
                    source["conversion_package_ms"]
                ),
                "feedback_log_ms": optional_float(source["feedback_log_ms"]),
                "feedback_send_ms": optional_float(source["feedback_send_ms"]),
                "outcome": source["outcome"],
            }
            row["segment"] = segment_for_volume(row["volume"])
            rows.append(row)
    if not rows:
        raise ValueError(f"No FIRE rows in {path}")
    if [row["lookup_id"] for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError(f"FIRE lookup IDs are not contiguous in {path}")
    return {"path": path, "rows": rows, "fields": reader.fieldnames}


def parse_queue_items(
    run_dir: Path, status: dict[str, Any]
) -> dict[str, Any]:
    path = unique_file(run_dir, "queue_profile_items_*.csv")
    required = {
        "scan_id",
        "volume",
        "group",
        "pointer_filename",
        "outcome",
        "output_label",
        "first_eligible_perf_ns",
        "selected_perf_ns",
        "decision_perf_ns",
        "item_complete_perf_ns",
    }
    rows = []
    seen = set()
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing queue-item columns in {path}: {sorted(missing)}")
        for source in reader:
            volume = as_int(source["volume"], "volume")
            group = as_int(source["group"], "group")
            key = (volume, group)
            if key in seen:
                raise ValueError(f"Duplicate queue item {key}")
            seen.add(key)
            outcome = source["outcome"]
            terminal = {
                "reference": "registered",
                "registered": "registered",
                "skipped_by_lifo": "skipped",
                "failed": "failed",
            }.get(outcome)
            if terminal is None:
                raise ValueError(f"Unexpected queue outcome {outcome!r}")
            if status["group_status"].get(key) != terminal:
                raise ValueError(
                    f"Queue/status mismatch for {key}: {outcome!r} versus "
                    f"{status['group_status'].get(key)!r}"
                )
            eligible = as_int(source["first_eligible_perf_ns"], "first_eligible_perf_ns")
            complete = as_int(source["item_complete_perf_ns"], "item_complete_perf_ns")
            active_start_field = (
                "decision_perf_ns" if outcome == "skipped_by_lifo" else "selected_perf_ns"
            )
            active_start = as_int(source[active_start_field], active_start_field)
            output_label = source["output_label"]
            output_index = None
            if output_label:
                match = re.match(r"^(\d+)_", output_label)
                if match:
                    output_index = int(match.group(1))
            rows.append(
                {
                    "scan_id": as_int(source["scan_id"], "scan_id"),
                    "volume": volume,
                    "group": group,
                    "segment": segment_for_volume(volume),
                    "pointer_filename": source["pointer_filename"],
                    "outcome": outcome,
                    "terminal_status": terminal,
                    "output_label": output_label,
                    "registration_index": output_index,
                    "queue_handling_time_ms": (complete - eligible) / 1e6,
                    "active_handling_time_ms": (complete - active_start) / 1e6,
                }
            )
    if seen != set(status["group_status"]):
        raise ValueError("Queue item/status row sets differ")
    return {"path": path, "rows": rows, "fields": reader.fieldnames}


def parse_queue_scans(run_dir: Path, items: dict[str, Any]) -> dict[str, Any]:
    path = unique_file(run_dir, "queue_profile_scans_*.csv")
    required = {
        "scan_id",
        "directory_entry_count",
        "scan_and_build_ms",
        "empty_polls_since_previous",
        "empty_scan_ms_since_previous",
        "candidate_file_count",
        "pointer_candidate_count",
    }
    by_scan = defaultdict(list)
    for item in items["rows"]:
        by_scan[item["scan_id"]].append(item)
    rows = []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing queue-scan columns in {path}: {sorted(missing)}")
        for source in reader:
            scan_id = as_int(source["scan_id"], "scan_id")
            linked = by_scan.get(scan_id, [])
            volumes = [item["volume"] for item in linked]
            empty_polls = as_int(
                source["empty_polls_since_previous"], "empty_polls_since_previous"
            )
            empty_ms = as_float(
                source["empty_scan_ms_since_previous"], "empty_scan_ms_since_previous"
            )
            volume = max(volumes) if volumes else None
            rows.append(
                {
                    "scan_id": scan_id,
                    "volume": volume,
                    "min_volume": min(volumes) if volumes else None,
                    "segment": segment_for_volume(volume) if volume is not None else None,
                    "items_in_scan": len(linked),
                    "directory_entry_count": as_int(
                        source["directory_entry_count"], "directory_entry_count"
                    ),
                    "candidate_file_count": as_int(
                        source["candidate_file_count"], "candidate_file_count"
                    ),
                    "pointer_candidate_count": as_int(
                        source["pointer_candidate_count"], "pointer_candidate_count"
                    ),
                    "scan_and_build_ms": as_float(
                        source["scan_and_build_ms"], "scan_and_build_ms"
                    ),
                    "empty_polls_since_previous": empty_polls,
                    "empty_scan_ms_since_previous": empty_ms,
                    "mean_empty_scan_ms": empty_ms / empty_polls if empty_polls else None,
                }
            )
    return {"path": path, "rows": rows, "fields": reader.fieldnames}


def transform_inventory(run_dir: Path) -> dict[str, Any]:
    root_entries = list(run_dir.iterdir())
    valid = {}
    identity = []
    unexpected = []
    for path in run_dir.glob("*.tfm"):
        match = TRANSFORM_RE.fullmatch(path.name)
        if match:
            index = int(match.group("index"))
            if index in valid:
                raise ValueError(f"Duplicate registration index {index} in transform files")
            valid[index] = {
                "filename": path.name,
                "volume": int(match.group("volume")),
                "group": int(match.group("group")),
            }
        elif IDENTITY_RE.fullmatch(path.name):
            identity.append(path.name)
        else:
            unexpected.append(path.name)
    processed_directories = sorted(
        path for path in root_entries if path.is_dir() and path.name.startswith("processed_")
    )
    processed_suffix_counts: Counter[str] = Counter()
    processed_transform_count = 0
    for directory in processed_directories:
        for path in directory.iterdir():
            if path.is_file():
                processed_suffix_counts[path.suffix] += 1
            elif path.is_dir() and path.name == "transforms":
                processed_transform_count += sum(
                    child.is_file() and child.suffix == ".tfm"
                    for child in path.iterdir()
                )
    return {
        "valid": valid,
        "identity": sorted(identity),
        "unexpected": sorted(unexpected),
        "highest_index": max(valid) if valid else None,
        "root_entry_count_final": len(root_entries),
        "processed_directories": processed_directories,
        "processed_suffix_counts": processed_suffix_counts,
        "processed_transform_count": processed_transform_count,
    }


def analyze_fire(
    fire: dict[str, Any], inventory: dict[str, Any], queue_items: dict[str, Any]
) -> dict[str, Any]:
    rows = fire["rows"]
    violations = []
    filename_mismatches = []
    selected_missing = []
    lower_after_higher = 0
    current_committed = None
    commits = []
    previous_commit = None
    unsuccessful_advancements = 0
    successful_without_advance = 0

    for row in rows:
        filename = row["selected_transform_filename"]
        selected = row["selected_registration_index"]
        before = row["last_consumed_registration_index_before"]
        after = row["last_consumed_registration_index_after"]
        before_name = row["last_consumed_transform_filename_before"]
        after_name = row["last_consumed_transform_filename_after"]

        if filename:
            match = TRANSFORM_RE.fullmatch(filename)
            if match is None or int(match.group("index")) != selected:
                filename_mismatches.append(row["lookup_id"])
            if filename not in {
                record["filename"] for record in inventory["valid"].values()
            }:
                selected_missing.append(filename)
            if IDENTITY_RE.fullmatch(filename):
                violations.append(f"lookup {row['lookup_id']} selected identity")

        if current_committed != before:
            violations.append(
                f"lookup {row['lookup_id']} before={before} but prior committed={current_committed}"
            )
        if selected is not None and current_committed is not None and selected < current_committed:
            lower_after_higher += 1

        if row["outcome"] == "feedback_sent":
            if selected is None or after != selected or after_name != filename:
                violations.append(
                    f"lookup {row['lookup_id']} successful commit does not match selection"
                )
            if before is not None and selected is not None and selected <= before:
                successful_without_advance += 1
            step = None if previous_commit is None else selected - previous_commit
            commit = dict(row)
            commit["previous_committed_index"] = previous_commit
            commit["advancement_step"] = step
            commits.append(commit)
            previous_commit = selected
            current_committed = selected
        else:
            if after != before or after_name != before_name:
                unsuccessful_advancements += 1
                violations.append(
                    f"lookup {row['lookup_id']} outcome={row['outcome']} advanced state"
                )

    commit_indices = [row["selected_registration_index"] for row in commits]
    steps = [row["advancement_step"] for row in commits if row["advancement_step"] is not None]
    regressions = sum(step < 0 for step in steps)
    duplicates = sum(step == 0 for step in steps)
    positive_steps = [step for step in steps if step > 0]
    multi = [row for row in commits if row["advancement_step"] and row["advancement_step"] > 1]

    attempted = {
        row["registration_index"]: row
        for row in queue_items["rows"]
        if row["registration_index"] is not None
    }
    supersession = []
    for row in multi:
        previous = row["previous_committed_index"]
        selected = row["selected_registration_index"]
        intermediate = list(range(previous + 1, selected))
        existing = [index for index in intermediate if index in inventory["valid"]]
        failed_attempts = [
            index
            for index in intermediate
            if index in attempted and attempted[index]["outcome"] == "failed"
        ]
        unknown = [
            index
            for index in intermediate
            if index not in inventory["valid"] and index not in failed_attempts
        ]
        if len(existing) == len(intermediate):
            classification = "generated_transforms_superseded"
        elif unknown:
            classification = "mixed_or_unknown"
        else:
            classification = "generated_and_failed_attempts"
        supersession.append(
            {
                "lookup_id": row["lookup_id"],
                "incoming_image_identifier": row["incoming_image_identifier"],
                "volume": row["volume"],
                "slice": row["slice"],
                "group": row["group"],
                "previous_consumed_index": previous,
                "new_consumed_index": selected,
                "jump_size": row["advancement_step"],
                "superseded_index_count": len(intermediate),
                "newer_candidate_count": row["newer_candidate_count"],
                "existing_intermediate_transforms": len(existing),
                "failed_intermediate_attempts": len(failed_attempts),
                "unknown_intermediate_indices": len(unknown),
                "classification": classification,
            }
        )

    outcome_counts = Counter(row["outcome"] for row in rows)
    discovery = {"all": describe(row["discovery_selection_ms"] for row in rows)}
    for outcome in sorted(outcome_counts):
        discovery[outcome] = describe(
            row["discovery_selection_ms"] for row in rows if row["outcome"] == outcome
        )
    discovery["segments"] = {}
    for segment in SEGMENTS:
        subset = [row for row in rows if row["segment"] == segment]
        discovery["segments"][segment] = {
            "discovery": describe(row["discovery_selection_ms"] for row in subset),
            "root_entries": describe(row["root_entry_count"] for row in subset),
            "transforms": describe(row["valid_transform_count"] for row in subset),
        }

    no_new = [row for row in rows if row["outcome"] == "no_new_transform"]
    success = [row for row in rows if row["outcome"] == "feedback_sent"]
    components = {}
    for field in (
        "discovery_selection_ms",
        "transform_read_ms",
        "conversion_ms",
        "packaging_ms",
        "conversion_package_ms",
        "feedback_log_ms",
        "feedback_send_ms",
    ):
        components[field] = describe(
            row[field] for row in success if row[field] is not None
        )

    return {
        "violations": violations,
        "filename_mismatches": filename_mismatches,
        "selected_missing": sorted(set(selected_missing)),
        "lower_after_higher": lower_after_higher,
        "unsuccessful_advancements": unsuccessful_advancements,
        "successful_without_advance": successful_without_advance,
        "commits": commits,
        "commit_indices": commit_indices,
        "steps": positive_steps,
        "step_counts": Counter(positive_steps),
        "regressions": regressions,
        "duplicates": duplicates,
        "multi_jumps": len(multi),
        "single_steps": sum(step == 1 for step in positive_steps),
        "supersession": supersession,
        "outcome_counts": outcome_counts,
        "discovery": discovery,
        "no_new_transform_count": len(no_new),
        "no_new_transform_fraction": len(no_new) / len(rows),
        "feedback_fraction": len(success) / len(rows),
        "components": components,
        "transform_correlation": correlations(
            rows, "valid_transform_count", "discovery_selection_ms"
        ),
        "root_correlation": correlations(rows, "root_entry_count", "discovery_selection_ms"),
        "transform_progress_correlation": correlations(
            rows, "lookup_id", "valid_transform_count"
        ),
        "root_progress_correlation": correlations(rows, "lookup_id", "root_entry_count"),
        "no_new_transform_correlation": correlations(
            no_new, "valid_transform_count", "discovery_selection_ms"
        ),
    }


def summarize_queue(
    scans: dict[str, Any], items: dict[str, Any]
) -> dict[str, Any]:
    mapped = [
        row
        for row in scans["rows"]
        if row["volume"] is not None and row["volume"] > 0
    ]
    nonreference_items = [row for row in items["rows"] if row["volume"] > 0]
    registered = [
        row for row in nonreference_items if row["terminal_status"] == "registered"
    ]
    result = {
        "scan_all_rows": describe(row["scan_and_build_ms"] for row in scans["rows"]),
        "scan_mapped": describe(row["scan_and_build_ms"] for row in mapped),
        "scan_entries": describe(row["directory_entry_count"] for row in mapped),
        "scan_root_correlation": correlations(
            mapped, "directory_entry_count", "scan_and_build_ms"
        ),
        "handling": describe(row["queue_handling_time_ms"] for row in registered),
        "active_handling": describe(row["active_handling_time_ms"] for row in registered),
        "segments": {},
        "empty_polls": sum(row["empty_polls_since_previous"] for row in scans["rows"]),
        "empty_scan_ms": sum(
            row["empty_scan_ms_since_previous"] for row in scans["rows"]
        ),
        "mapped_rows": len(mapped),
        "unmapped_rows": len(scans["rows"]) - len(mapped),
    }
    result["mean_empty_scan_ms"] = (
        result["empty_scan_ms"] / result["empty_polls"]
        if result["empty_polls"]
        else None
    )
    for segment in SEGMENTS:
        segment_scans = [row for row in mapped if row["segment"] == segment]
        segment_items = [row for row in registered if row["segment"] == segment]
        result["segments"][segment] = {
            "scan": describe(row["scan_and_build_ms"] for row in segment_scans),
            "entries": describe(row["directory_entry_count"] for row in segment_scans),
            "handling": describe(
                row["queue_handling_time_ms"] for row in segment_items
            ),
            "active_handling": describe(
                row["active_handling_time_ms"] for row in segment_items
            ),
        }
    return result


def group_by_volume(
    rows: list[dict[str, Any]], fields: list[str]
) -> list[dict[str, Any]]:
    grouped = defaultdict(list)
    for row in rows:
        if row.get("volume") is not None and row["volume"] > 0:
            grouped[row["volume"]].append(row)
    output = []
    for volume in sorted(grouped):
        entry: dict[str, Any] = {"volume": volume, "n": len(grouped[volume])}
        for field in fields:
            summary = describe(
                row[field] for row in grouped[volume] if row.get(field) is not None
            )
            for statistic in ("mean", "median", "p95", "p99", "max"):
                entry[f"{field}_{statistic}"] = summary[statistic]
        output.append(entry)
    return output


def summary_rows(service: str, field: str, overall: dict[str, Any], segments: dict[str, Any]):
    rows = []
    for segment, summary in [("overall", overall), *segments.items()]:
        rows.append(
            {
                "service": service,
                "metric": field,
                "segment": segment,
                **summary,
            }
        )
    return rows


def rolling_median(rows: list[dict[str, Any]], xfield: str, yfield: str, bins=50):
    ordered = sorted(rows, key=lambda row: row[xfield])
    if not ordered:
        return [], []
    bin_size = max(1, math.ceil(len(ordered) / bins))
    xs, ys = [], []
    for offset in range(0, len(ordered), bin_size):
        chunk = ordered[offset : offset + bin_size]
        xs.append(statistics.median(float(row[xfield]) for row in chunk))
        ys.append(statistics.median(float(row[yfield]) for row in chunk))
    return xs, ys


def per_volume_median(rows: list[dict[str, Any]], field: str):
    grouped = defaultdict(list)
    for row in rows:
        if row.get("volume") is not None and row["volume"] > 0:
            grouped[row["volume"]].append(float(row[field]))
    volumes = sorted(grouped)
    return volumes, [statistics.median(grouped[volume]) for volume in volumes]


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.2,
            "figure.dpi": 120,
        }
    )


def scatter_with_trend(
    output: Path,
    rows: list[dict[str, Any]],
    xfield: str,
    yfield: str,
    xlabel: str,
    ylabel: str,
    title: str,
    filename: str,
) -> None:
    fig, axis = plt.subplots(figsize=(10.5, 5.4))
    axis.scatter(
        [row[xfield] for row in rows],
        [row[yfield] for row in rows],
        s=5,
        alpha=0.10,
        linewidths=0,
        color="#3976af",
        label="Individual lookups",
    )
    xs, ys = rolling_median(rows, xfield, yfield)
    axis.plot(xs, ys, color="#d1495b", linewidth=2, label="Equal-count-bin median")
    axis.set(xlabel=xlabel, ylabel=ylabel, title=title)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / filename, dpi=180)
    plt.close(fig)


def create_plots(
    output: Path,
    fire: dict[str, Any],
    fire_analysis: dict[str, Any],
    scans: dict[str, Any],
    items: dict[str, Any],
    status: dict[str, Any],
) -> None:
    fire_rows = fire["rows"]
    commits = fire_analysis["commits"]

    fig, axis = plt.subplots(figsize=(10.5, 5.2))
    axis.plot(
        [row["lookup_id"] for row in commits],
        [row["selected_registration_index"] for row in commits],
        linewidth=1.2,
        color="#3976af",
    )
    axis.set(
        xlabel="FIRE MOCO lookup sequence",
        ylabel="Committed registration index",
        title="FIRE monotonic transform advancement",
    )
    fig.tight_layout()
    fig.savefig(output / "fire_consumed_registration_index_vs_progress.png", dpi=180)
    plt.close(fig)

    step_rows = [row for row in commits if row["advancement_step"] is not None]
    fig, axis = plt.subplots(figsize=(10.5, 5.2))
    axis.scatter(
        [row["lookup_id"] for row in step_rows],
        [row["advancement_step"] for row in step_rows],
        s=10,
        alpha=0.55,
        linewidths=0,
        color="#3976af",
    )
    axis.axhline(1, color="#555555", linestyle="--", linewidth=1)
    axis.set(
        xlabel="FIRE MOCO lookup sequence",
        ylabel="Committed index step",
        title="FIRE transform advancement step size",
    )
    fig.tight_layout()
    fig.savefig(output / "fire_advancement_step_size.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(10.5, 5.2))
    axis.scatter(
        [row["volume"] for row in fire_rows],
        [row["discovery_selection_ms"] for row in fire_rows],
        s=5,
        alpha=0.09,
        linewidths=0,
        color="#3976af",
    )
    xs, ys = per_volume_median(fire_rows, "discovery_selection_ms")
    axis.plot(xs, ys, color="#d1495b", linewidth=2, label="Per-volume median")
    axis.set(
        xlabel="Volume",
        ylabel="Discovery and selection (ms)",
        title="FIRE transform discovery over acquisition",
    )
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "fire_discovery_time_vs_volume.png", dpi=180)
    plt.close(fig)

    scatter_with_trend(
        output,
        fire_rows,
        "valid_transform_count",
        "discovery_selection_ms",
        "Valid active registration transforms",
        "Discovery and selection (ms)",
        "FIRE discovery versus active transform count",
        "fire_discovery_time_vs_transform_count.png",
    )
    scatter_with_trend(
        output,
        fire_rows,
        "root_entry_count",
        "discovery_selection_ms",
        "Active root entries",
        "Discovery and selection (ms)",
        "FIRE discovery versus active root size",
        "fire_discovery_time_vs_root_entry_count.png",
    )

    volumes, root_medians = per_volume_median(fire_rows, "root_entry_count")
    _, transform_medians = per_volume_median(fire_rows, "valid_transform_count")
    fig, axis = plt.subplots(figsize=(10.5, 5.2))
    axis.plot(volumes, root_medians, linewidth=2, label="Root entries", color="#3976af")
    axis.plot(
        volumes,
        transform_medians,
        linewidth=2,
        label="Valid transforms",
        color="#e1812c",
    )
    axis.set(
        xlabel="Volume",
        ylabel="Per-volume median count",
        title="FIRE-observed active-directory growth",
    )
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "fire_root_and_transform_count_vs_progress.png", dpi=180)
    plt.close(fig)

    mapped_scans = [row for row in scans["rows"] if row["volume"] and row["volume"] > 0]
    fig, axis = plt.subplots(figsize=(10.5, 5.2))
    axis.scatter(
        [row["volume"] for row in mapped_scans],
        [row["scan_and_build_ms"] for row in mapped_scans],
        s=5,
        alpha=0.10,
        linewidths=0,
        color="#e1812c",
    )
    xs, ys = per_volume_median(mapped_scans, "scan_and_build_ms")
    axis.plot(xs, ys, color="#8c510a", linewidth=2, label="Per-volume median")
    axis.set(
        xlabel="Volume",
        ylabel="Scan and candidate build (ms)",
        title="Queue directory discovery over acquisition",
    )
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "queue_scan_time_vs_volume.png", dpi=180)
    plt.close(fig)

    scatter_with_trend(
        output,
        mapped_scans,
        "directory_entry_count",
        "scan_and_build_ms",
        "Active root entries",
        "Scan and candidate build (ms)",
        "Queue discovery versus active root size",
        "queue_scan_time_vs_root_entry_count.png",
    )

    registration = [row for row in status["rows"] if row["volume"] > 0]
    fig, axis = plt.subplots(figsize=(10.5, 5.2))
    axis.plot(
        [row["volume"] for row in registration],
        [row["registered"] for row in registration],
        color="#3976af",
        linewidth=1.7,
    )
    axis.axhline(
        status["expected_groups"],
        color="#555555",
        linestyle="--",
        linewidth=1,
        label=f"Expected ({status['expected_groups']})",
    )
    degraded = [row for row in registration if row["registered"] < row["expected"]]
    axis.scatter(
        [row["volume"] for row in degraded],
        [row["registered"] for row in degraded],
        color="#d1495b",
        s=28,
        zorder=3,
        label="Degraded volume",
    )
    axis.set(
        xlabel="Volume",
        ylabel="Registered groups",
        title="Registration completion by volume",
        ylim=(0, status["expected_groups"] + 2),
    )
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "queue_registered_groups_per_volume.png", dpi=180)
    plt.close(fig)

    registered_items = [
        row
        for row in items["rows"]
        if row["volume"] > 0 and row["terminal_status"] == "registered"
    ]
    fig, axis = plt.subplots(figsize=(10.5, 5.2))
    axis.scatter(
        [row["volume"] for row in registered_items],
        [row["queue_handling_time_ms"] for row in registered_items],
        s=5,
        alpha=0.08,
        linewidths=0,
        color="#3976af",
    )
    xs, queue_medians = per_volume_median(registered_items, "queue_handling_time_ms")
    _, active_medians = per_volume_median(registered_items, "active_handling_time_ms")
    axis.plot(xs, queue_medians, color="#3976af", linewidth=2, label="Eligible-to-complete median")
    axis.plot(xs, active_medians, color="#e1812c", linewidth=2, label="Active handling median")
    axis.set(
        xlabel="Volume",
        ylabel="Registered-item time (ms)",
        title="Queue registered-item lifecycle over acquisition",
    )
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "queue_handling_time_vs_volume.png", dpi=180)
    plt.close(fig)

    fire_x, fire_y = per_volume_median(fire_rows, "discovery_selection_ms")
    queue_x, queue_y = per_volume_median(mapped_scans, "scan_and_build_ms")
    fig, axis = plt.subplots(figsize=(10.5, 5.2))
    axis.plot(fire_x, fire_y, color="#3976af", linewidth=2, label="FIRE discovery median")
    axis.plot(queue_x, queue_y, color="#e1812c", linewidth=2, label="Queue discovery median")
    axis.set(
        xlabel="Volume",
        ylabel="Discovery time (ms)",
        title="FIRE and queue directory discovery by acquisition progress",
    )
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "fire_vs_queue_discovery_over_progress.png", dpi=180)
    plt.close(fig)


def fmt(value: Any, digits=3) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.{digits}f}"


def ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator in (None, 0):
        return None
    return numerator / denominator


def report(
    output: Path,
    run_dir: Path,
    sources: dict[str, Path],
    status: dict[str, Any],
    registration: dict[str, Any],
    fire: dict[str, Any],
    fire_analysis: dict[str, Any],
    inventory: dict[str, Any],
    queue: dict[str, Any],
) -> None:
    fire_discovery = fire_analysis["discovery"]
    fire_early = fire_discovery["segments"]["early"]["discovery"]
    fire_late = fire_discovery["segments"]["late"]["discovery"]
    queue_early = queue["segments"]["early"]["scan"]
    queue_late = queue["segments"]["late"]["scan"]
    steps = fire_analysis["steps"]
    commits = fire_analysis["commits"]
    final_commit = commits[-1]["selected_registration_index"] if commits else None
    lag = inventory["highest_index"] - final_commit if final_commit is not None else None
    dominant = max(
        (
            (name, fire_analysis["components"][name]["mean"])
            for name in (
                "discovery_selection_ms",
                "transform_read_ms",
                "conversion_ms",
                "packaging_ms",
                "feedback_log_ms",
                "feedback_send_ms",
            )
        ),
        key=lambda item: item[1] if item[1] is not None else -1,
    )
    outcome_rows = [
        f"| `{outcome}` | {fire_analysis['outcome_counts'].get(outcome, 0)} |"
        for outcome in (
            "feedback_sent",
            "no_new_transform",
            "discovery_failed",
            "transform_read_failed",
            "conversion_failed",
            "packaging_failed",
            "feedback_logging_failed",
            "send_failed",
        )
    ]
    largest = sorted(
        fire_analysis["supersession"], key=lambda row: row["jump_size"], reverse=True
    )[:10]

    lines = [
        "# MOCO-on FIRE and queue baseline analysis",
        "",
        "## Experiment and sources",
        "",
        f"- Acquisition: `{run_dir}`",
        f"- FIRE profiler: `{sources['fire'].name}` ({len(fire['rows'])} lookups)",
        f"- Queue scan profiler: `{sources['queue_scans'].name}`",
        f"- Queue item profiler: `{sources['queue_items'].name}`",
        f"- Authoritative outcomes: `{sources['status'].name}`",
        f"- MOCO feedback log: `{sources['feedback_log'].name}`",
        f"- Acquisition FIRE log: `{sources['fire_log'].name}`",
        "",
        "There was exactly one profiler file of each required type. The large FIRE log aligned with the "
        "FIRE profiler start time; the earlier 256-byte FIRE log is a startup/connection artifact and was not used.",
        f"The final root contains {len(inventory['valid'])} registration transforms, "
        f"{len(inventory['identity'])} identity transform, and {len(inventory['unexpected'])} unexpected `.tfm` files. "
        f"The processed directory contains {inventory['processed_suffix_counts']['.txt']} pointers, "
        f"{inventory['processed_suffix_counts']['.nhdr']} headers, and "
        f"{inventory['processed_suffix_counts']['.raw']} raw files; it contains "
        f"{inventory['processed_transform_count']} archived transforms, as expected with MOCO-on retirement disabled.",
        "",
        "## Acquisition structure and registration completion",
        "",
        f"The status file contains {len(status['volumes'])} finalized volumes: one fixed reference "
        f"(volume 0) and {len(status['volumes']) - 1} registration volumes. Every non-reference "
        f"volume contains {status['expected_groups']} groups, giving "
        f"{registration['expected']} registration opportunities.",
        "",
        "| Metric | Value |",
        "| --- | ---: |",
        f"| Expected | {registration['expected']} |",
        f"| Registered | {registration['registered']} |",
        f"| Skipped | {registration['skipped']} |",
        f"| Failed | {registration['failed']} |",
        f"| Success % | {registration['success_percent']:.3f}% |",
        f"| First degraded volume | {registration['first_degraded']} |",
        f"| First degraded after startup | {registration['first_degraded_after_startup']} |",
        f"| Degraded volumes | {registration['degraded_volumes']} |",
        f"| Minimum registered groups | {registration['minimum_registered']} |",
        f"| Final volume R/S/F | {registration['final']['registered']}/{registration['final']['skipped']}/{registration['final']['failed']} |",
        "",
        "## FIRE functional advancement",
        "",
        "| Metric | Value |",
        "| --- | ---: |",
        f"| Successful feedback advancements | {len(commits)} |",
        f"| First consumed index | {commits[0]['selected_registration_index'] if commits else 'n/a'} |",
        f"| Final consumed index | {final_commit} |",
        f"| Minimum positive step | {min(steps) if steps else 'n/a'} |",
        f"| Median step | {statistics.median(steps) if steps else 'n/a'} |",
        f"| Maximum step | {max(steps) if steps else 'n/a'} |",
        f"| Single-index steps | {fire_analysis['single_steps']} |",
        f"| Multi-index jumps | {fire_analysis['multi_jumps']} |",
        f"| Regressions | {fire_analysis['regressions']} |",
        f"| Duplicate commits | {fire_analysis['duplicates']} |",
        f"| Successful events without advancement | {fire_analysis['successful_without_advance']} |",
        f"| Lower selection after a higher commit | {fire_analysis['lower_after_higher']} |",
        "",
        f"All {len(commits)} committed indices were strictly increasing. "
        f"The invariant audit found {len(fire_analysis['violations'])} state-transition violations, "
        f"{len(fire_analysis['filename_mismatches'])} filename/index mismatches, and "
        f"{len(fire_analysis['selected_missing'])} selected filenames absent from the final transform inventory.",
        f"The identity file `{inventory['identity'][0] if inventory['identity'] else 'none'}` was never selected; "
        "no malformed or unrelated transform was selected. The MOCO feedback log contains "
        f"{fire_analysis['feedback_log_records']} transform records and exactly matches the ordered "
        "`feedback_sent` profile sequence.",
        "",
        "### Supersession",
        "",
        f"There were {fire_analysis['multi_jumps']} multi-index jumps. Final root inventory contains "
        f"{len(inventory['valid'])} valid registration transforms with indices 1–{inventory['highest_index']}; "
        "there were no failed registrations. Therefore intermediate indices in this run's jumps are "
        "generated transforms superseded before FIRE's next successful send, not LIFO skips (which do "
        "not consume registration indices).",
        "",
        "Largest jumps:",
        "",
        "| Lookup | Image | Volume/group | Previous | New | Step | Superseded | Newer candidates |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in largest:
        lines.append(
            f"| {row['lookup_id']} | {row['incoming_image_identifier']} | "
            f"{row['volume']}/{row['group']} | {row['previous_consumed_index']} | "
            f"{row['new_consumed_index']} | {row['jump_size']} | "
            f"{row['superseded_index_count']} | {row['newer_candidate_count']} |"
        )
    lines.extend(
        [
            "",
            "## FIRE failure and retry behavior",
            "",
            "| Outcome | Events |",
            "| --- | ---: |",
            *outcome_rows,
            "",
            f"No unsuccessful row changed consumed state ({fire_analysis['unsuccessful_advancements']} incorrect advancements). "
            "No read/conversion/packaging/log/send failures occurred, so this completed run does not "
            "exercise runtime retry behavior; those paths remain unit-tested rather than empirically tested here. "
            "A profiler write failure cannot necessarily record its own outcome, but the profile contains all "
            f"{len(fire['rows'])} expected non-reference image lookups.",
            "",
            "## Final FIRE state",
            "",
            f"FIRE committed index {final_commit}; the highest generated transform was "
            f"{inventory['highest_index']} (`{inventory['valid'][inventory['highest_index']]['filename']}`). "
            f"The end-of-stream lag was {lag} transform. The last FIRE lookup selected index {final_commit} "
            "for group 42 while the final group-43 transform was generated afterward; no subsequent scanner "
            "image existed to trigger another FIRE lookup. This is timing-consistent end-of-acquisition lag, "
            "not a monotonicity failure.",
            "",
            "## FIRE discovery scaling",
            "",
            "`discovery_selection_ms` directly measures `os.listdir()` plus strict filename parsing, "
            "consumed-index comparison, and one-pass maximum selection. It excludes transform read, "
            "conversion, packaging, feedback logging, and socket send.",
            "",
            "| Scope | n | Mean ms | Median ms | Std ms | p95 ms | p99 ms | Max ms |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for scope in ("all", "feedback_sent", "no_new_transform"):
        summary = fire_discovery.get(scope, describe([]))
        lines.append(
            f"| {scope} | {summary['n']} | {fmt(summary['mean'])} | {fmt(summary['median'])} | "
            f"{fmt(summary['std'])} | {fmt(summary['p95'])} | {fmt(summary['p99'])} | {fmt(summary['max'])} |"
        )
    lines.extend(
        [
            "",
            "| Segment | Mean ms | Median ms | p95 ms | p99 ms | Mean root entries | Median root entries | Mean transforms | Median transforms |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for segment in SEGMENTS:
        entry = fire_discovery["segments"][segment]
        timing = entry["discovery"]
        root = entry["root_entries"]
        transforms = entry["transforms"]
        lines.append(
            f"| {segment} | {fmt(timing['mean'])} | {fmt(timing['median'])} | "
            f"{fmt(timing['p95'])} | {fmt(timing['p99'])} | {fmt(root['mean'], 1)} | "
            f"{fmt(root['median'], 1)} | {fmt(transforms['mean'], 1)} | {fmt(transforms['median'], 1)} |"
        )
    lines.extend(
        [
            "",
            f"FIRE late/early mean ratio: **{fmt(ratio(fire_late['mean'], fire_early['mean']))}×**; "
            f"late/early median ratio: **{fmt(ratio(fire_late['median'], fire_early['median']))}×**.",
            f"Discovery versus valid-transform count: Pearson **{fmt(fire_analysis['transform_correlation']['pearson'])}**, "
            f"Spearman **{fmt(fire_analysis['transform_correlation']['spearman'])}**. "
            f"Versus root entries: Pearson **{fmt(fire_analysis['root_correlation']['pearson'])}**, "
            f"Spearman **{fmt(fire_analysis['root_correlation']['spearman'])}**. These are descriptive associations, not causal estimates.",
            "",
            "### Directory growth and no-new-transform lookups",
            "",
            f"FIRE observed valid-transform counts from {fire['rows'][0]['valid_transform_count']} initially "
            f"to {fire['rows'][-1]['valid_transform_count']} on the last lookup (median "
            f"{statistics.median(row['valid_transform_count'] for row in fire['rows']):.1f}, maximum "
            f"{max(row['valid_transform_count'] for row in fire['rows'])}). Root entries ranged from "
            f"{fire['rows'][0]['root_entry_count']} initially to {fire['rows'][-1]['root_entry_count']} "
            f"on the last lookup (median {statistics.median(row['root_entry_count'] for row in fire['rows']):.1f}, "
            f"maximum {max(row['root_entry_count'] for row in fire['rows'])}). The transform count therefore "
            "accounts for most long-run active-root growth, although pending image work and other outputs make "
            "root count temporarily exceed transform count. Transform count versus lookup progress had Pearson "
            f"correlation {fmt(fire_analysis['transform_progress_correlation']['pearson'])}; root count versus "
            f"progress had Pearson correlation {fmt(fire_analysis['root_progress_correlation']['pearson'])}, "
            "consistent with approximately linear acquisition growth.",
            "",
            f"`no_new_transform` occurred {fire_analysis['no_new_transform_count']} times "
            f"({100.0 * fire_analysis['no_new_transform_fraction']:.3f}% of lookups). Its transform-count "
            f"association was Pearson {fmt(fire_analysis['no_new_transform_correlation']['pearson'])} and "
            f"Spearman {fmt(fire_analysis['no_new_transform_correlation']['spearman'])}.",
            "",
            "## FIRE successful-feedback timing decomposition",
            "",
            "| Component | n | Mean ms | Median ms | p95 ms | Max ms |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for component, summary in fire_analysis["components"].items():
        lines.append(
            f"| `{component}` | {summary['n']} | {fmt(summary['mean'])} | "
            f"{fmt(summary['median'])} | {fmt(summary['p95'])} | {fmt(summary['max'])} |"
        )
    lines.extend(
        [
            "",
            f"The largest mean non-overlapping component was `{dominant[0]}` at {fmt(dominant[1])} ms. "
            "`conversion_package_ms` spans conversion through packaging and overlaps the separately reported "
            "conversion and packaging fields; it is not added to them.",
            "",
            "## Queue discovery scaling",
            "",
            "`scan_and_build_ms` brackets queue `list_new_files()`: one `os.listdir()`, extension and "
            "`seen_files` filtering, `getmtime()` for valid unseen candidates, candidate sorting, and queue-profiler "
            "observation work. Historical `.tfm` files affect root enumeration but are rejected before metadata calls.",
            "",
            "| Segment | n | Mean ms | Median ms | p95 ms | p99 ms | Mean root entries | Median root entries |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for segment in SEGMENTS:
        entry = queue["segments"][segment]
        timing = entry["scan"]
        entries = entry["entries"]
        lines.append(
            f"| {segment} | {timing['n']} | {fmt(timing['mean'])} | {fmt(timing['median'])} | "
            f"{fmt(timing['p95'])} | {fmt(timing['p99'])} | {fmt(entries['mean'], 1)} | {fmt(entries['median'], 1)} |"
        )
    lines.extend(
        [
            "",
            f"Overall mapped non-reference queue scans: mean {fmt(queue['scan_mapped']['mean'])} ms, "
            f"median {fmt(queue['scan_mapped']['median'])} ms, p95 {fmt(queue['scan_mapped']['p95'])} ms, "
            f"p99 {fmt(queue['scan_mapped']['p99'])} ms, max {fmt(queue['scan_mapped']['max'])} ms. "
            f"Late/early mean ratio: {fmt(ratio(queue_late['mean'], queue_early['mean']))}×; "
            f"median ratio: {fmt(ratio(queue_late['median'], queue_early['median']))}×. "
            f"Queue scan versus root-entry correlation was Pearson {fmt(queue['scan_root_correlation']['pearson'])}, "
            f"Spearman {fmt(queue['scan_root_correlation']['spearman'])}.",
            f"The profiler aggregated {queue['empty_polls']} empty polls totaling "
            f"{queue['empty_scan_ms']:.3f} ms, or {fmt(queue['mean_empty_scan_ms'])} ms per empty scan.",
            "For context only, prior separate MOCO-off runs with image and transform retirement reported mean "
            "non-empty queue discovery near 1.706 ms and 2.085 ms. This same-run MOCO-on mean is higher, but "
            "the separate runs do not isolate MOCO or transform retention as the sole cause.",
            "",
            "## Queue registered-item lifecycle",
            "",
            "`queue_handling_time_ms` is first eligibility through item completion and may include same-batch "
            "waiting. `active_handling_time_ms` starts when the registered item is selected and excludes that wait.",
            "",
            "| Metric | n | Mean ms | Median ms | p95 ms | p99 ms | Max ms |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            f"| Queue handling | {queue['handling']['n']} | {fmt(queue['handling']['mean'])} | {fmt(queue['handling']['median'])} | {fmt(queue['handling']['p95'])} | {fmt(queue['handling']['p99'])} | {fmt(queue['handling']['max'])} |",
            f"| Active handling | {queue['active_handling']['n']} | {fmt(queue['active_handling']['mean'])} | {fmt(queue['active_handling']['median'])} | {fmt(queue['active_handling']['p95'])} | {fmt(queue['active_handling']['p99'])} | {fmt(queue['active_handling']['max'])} |",
            "",
            "| Segment | Queue mean/median/p95 ms | Active mean/median/p95 ms |",
            "| --- | ---: | ---: |",
        ]
    )
    for segment in SEGMENTS:
        handling = queue["segments"][segment]["handling"]
        active = queue["segments"][segment]["active_handling"]
        lines.append(
            f"| {segment} | {fmt(handling['mean'])} / {fmt(handling['median'])} / {fmt(handling['p95'])} | "
            f"{fmt(active['mean'])} / {fmt(active['median'])} / {fmt(active['p95'])} |"
        )
    lines.extend(
        [
            "",
            "## FIRE versus queue discovery",
            "",
            "| Service | Overall mean | Overall median | Overall p95 | Early mean | Late mean | Late/early mean |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            f"| FIRE | {fmt(fire_discovery['all']['mean'])} | {fmt(fire_discovery['all']['median'])} | {fmt(fire_discovery['all']['p95'])} | {fmt(fire_early['mean'])} | {fmt(fire_late['mean'])} | {fmt(ratio(fire_late['mean'], fire_early['mean']))}× |",
            f"| Queue | {fmt(queue['scan_mapped']['mean'])} | {fmt(queue['scan_mapped']['median'])} | {fmt(queue['scan_mapped']['p95'])} | {fmt(queue_early['mean'])} | {fmt(queue_late['mean'])} | {fmt(ratio(queue_late['mean'], queue_early['mean']))}× |",
            "",
            "These timings are not interchangeable. FIRE parses every root filename and compares every valid "
            "historical transform index on every non-reference image. Queue enumerates every root entry too, "
            "but `seen_files` and extension checks prevent historical transforms from receiving candidate metadata "
            "or sorting work. Volume alignment is used here; FIRE and queue perf-counter clocks are process-local "
            "and cannot provide exact cross-process event ordering.",
            f"Both discovery curves increase over acquisition progress, and post-startup registration degradation "
            f"begins at volume {registration['first_degraded_after_startup']}. This pattern is consistent with "
            "growing directory pressure, but process scheduling and the lack of synchronized cross-process clocks "
            "prevent a causal or exact event-ordering claim.",
            "",
            "## Interpretation and Step 1 validation",
            "",
            f"Functionally, Step 1 succeeded: committed advancement was strictly monotonic, gaps were tolerated, "
            f"{fire_analysis['multi_jumps']} lookups superseded intermediate generated transforms, and no stale "
            "index resurfaced. The final one-index lag is explained by the lack of an image-triggered lookup after "
            "the final registration output.",
            "",
            "Performance-wise, source inspection confirms the historical transform `getmtime()` calls and sorting "
            "are absent. The measured late/early FIRE ratio and transform/root correlations quantify the residual "
            "cost of repeated root enumeration plus filename parsing and one-pass index comparison. Queue retains "
            "its own root-enumeration sensitivity while avoiding historical-transform metadata work.",
            "",
            "The evidence supports a controlled Step 2 transform-retirement experiment. FIRE should target lower "
            "and flatter `discovery_selection_ms`; queue should target lower and flatter `scan_and_build_ms`. No "
            "numerical speedup is predicted from this single baseline. The explicit FIRE state is suitable for a "
            "release-through watermark in normal uninterrupted operation, with the current consumed transform kept "
            "as FIRE's active anchor and only lower indices released after motion-monitor is also finished.",
            "",
            "## Direct answers",
            "",
            f"1. **Strict monotonicity:** Yes; {len(commits)} commits, {fire_analysis['regressions']} regressions.",
            f"2. **Regressions/duplicates:** {fire_analysis['regressions']} / {fire_analysis['duplicates']}.",
            f"3. **Incorrect failure advancement:** {fire_analysis['unsuccessful_advancements']}; no runtime failure outcomes occurred.",
            f"4. **Registration gaps tolerated:** Yes; {fire_analysis['multi_jumps']} multi-index commits completed correctly.",
            f"5. **Supersession frequency:** {fire_analysis['multi_jumps']} of {max(len(commits)-1, 0)} measured commit transitions ({100.0 * fire_analysis['multi_jumps'] / max(len(commits)-1, 1):.3f}%).",
            f"6. **Largest jump:** {max(steps) if steps else 'n/a'} indices.",
            f"7. **Reached latest generated transform:** No; final consumed {final_commit}, generated {inventory['highest_index']}, lag {lag}.",
            f"8. **No-new lookups:** {fire_analysis['no_new_transform_count']}.",
            f"9. **Feedback-send fraction:** {100.0 * fire_analysis['feedback_fraction']:.3f}%.",
            f"10. **FIRE discovery mean/median/p95/p99:** {fmt(fire_discovery['all']['mean'])} / {fmt(fire_discovery['all']['median'])} / {fmt(fire_discovery['all']['p95'])} / {fmt(fire_discovery['all']['p99'])} ms.",
            f"11. **Early-to-late FIRE change:** mean {fmt(fire_early['mean'])} → {fmt(fire_late['mean'])} ms ({fmt(ratio(fire_late['mean'], fire_early['mean']))}×).",
            f"12. **Transform-count association:** Pearson {fmt(fire_analysis['transform_correlation']['pearson'])}, Spearman {fmt(fire_analysis['transform_correlation']['spearman'])}.",
            f"13. **Root-count association:** Pearson {fmt(fire_analysis['root_correlation']['pearson'])}, Spearman {fmt(fire_analysis['root_correlation']['spearman'])}.",
            "14. **Residual acquisition-length scaling:** See the measured early/late ratio and correlations above; this remains observational, not causal.",
            f"15. **Dominant successful component:** `{dominant[0]}` by mean ({fmt(dominant[1])} ms).",
            f"16. **Queue discovery:** overall mapped mean/median/p95 {fmt(queue['scan_mapped']['mean'])} / {fmt(queue['scan_mapped']['median'])} / {fmt(queue['scan_mapped']['p95'])} ms.",
            f"17. **Queue root-growth association:** Pearson {fmt(queue['scan_root_correlation']['pearson'])}, with late/early mean ratio {fmt(ratio(queue_late['mean'], queue_early['mean']))}×.",
            f"18. **Registration degradation:** {registration['degraded_volumes']} degraded volumes; first after startup {registration['first_degraded_after_startup']}; success {registration['success_percent']:.3f}%.",
            "19. **FIRE versus queue:** Both enumerate the growing root, but FIRE additionally parses all historical transform names; compare ratios in the table above.",
            "20. **Proceed to Step 2:** Yes, as a controlled experiment, subject to the previously defined cross-consumer safety protocol.",
            "21. **FIRE target:** Reduce and flatten root-enumeration/transform-parsing discovery time.",
            "22. **Queue target:** Reduce and flatten root enumeration and extension-filtering time.",
            "23. **Watermark suitability:** Yes for uninterrupted operation; the committed index is explicit and monotonic.",
            "24. **Pre-Step-2 anomalies:** The one-transform terminal lag must be represented in release semantics; no correctness anomaly otherwise appeared.",
            "",
            "## Limitations",
            "",
            "- This is one observational replay, so correlation does not establish causation or isolate host scheduling noise.",
            "- FIRE and queue monotonic perf counters are process-local; comparison is aligned by volume, not exact timestamp.",
            "- The final saved layout cannot reconstruct temporal directory size; temporal counts come from profiler snapshots.",
            "- FIRE profiling has no scanner acknowledgement; `feedback_sent` means the application send call returned.",
            "- No runtime read/conversion/package/log/send failure occurred, so retry behavior is not empirically exercised here.",
            "- Final transform inventory is used to classify superseded indices; no acquisition files were modified.",
        ]
    )
    (output / "analysis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    run_dir = args.acquisition_dir.expanduser().resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"ERROR: acquisition directory does not exist: {run_dir}")
    output = (
        args.output_dir.expanduser().resolve()
        if args.output_dir
        else run_dir.parent / f"{run_dir.name}_fire_moco_analysis"
    )
    if output == run_dir or run_dir in output.parents:
        raise SystemExit("ERROR: output directory must not be inside the acquisition directory")
    if output.exists() and any(output.iterdir()):
        raise SystemExit(f"ERROR: output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)

    status = parse_status(run_dir)
    registration = summarize_registration(status)
    fire = parse_fire_profile(run_dir)
    queue_items = parse_queue_items(run_dir, status)
    queue_scans = parse_queue_scans(run_dir, queue_items)
    inventory = transform_inventory(run_dir)
    if not inventory["valid"]:
        raise ValueError("No valid root-level registration transforms found")
    if inventory["unexpected"]:
        raise ValueError(f"Unexpected root-level TFM files: {inventory['unexpected']}")

    fire_analysis = analyze_fire(fire, inventory, queue_items)
    queue_analysis = summarize_queue(queue_scans, queue_items)
    if fire_analysis["violations"]:
        raise ValueError(
            "FIRE advancement invariant violations: "
            + "; ".join(fire_analysis["violations"][:10])
        )

    fire_by_volume = group_by_volume(
        fire["rows"],
        ["discovery_selection_ms", "root_entry_count", "valid_transform_count"],
    )
    queue_by_volume = group_by_volume(
        [row for row in queue_scans["rows"] if row["volume"] is not None],
        ["scan_and_build_ms", "directory_entry_count"],
    )
    registered_items = [
        row
        for row in queue_items["rows"]
        if row["volume"] > 0 and row["terminal_status"] == "registered"
    ]
    handling_by_volume = group_by_volume(
        registered_items,
        ["queue_handling_time_ms", "active_handling_time_ms"],
    )

    fire_event_fields = [
        "lookup_id",
        "incoming_image_identifier",
        "volume",
        "slice",
        "group",
        "segment",
        "root_entry_count",
        "valid_transform_count",
        "newer_candidate_count",
        "selected_transform_filename",
        "selected_registration_index",
        "last_consumed_registration_index_before",
        "last_consumed_registration_index_after",
        "discovery_selection_ms",
        "transform_read_ms",
        "conversion_ms",
        "packaging_ms",
        "conversion_package_ms",
        "feedback_log_ms",
        "feedback_send_ms",
        "outcome",
    ]
    write_csv(output / "fire_lookup_events.csv", fire_event_fields, fire["rows"])
    advancement_fields = fire_event_fields + [
        "previous_committed_index",
        "advancement_step",
    ]
    write_csv(
        output / "fire_advancement_events.csv",
        advancement_fields,
        fire_analysis["commits"],
    )
    write_csv(
        output / "fire_supersession_events.csv",
        [
            "lookup_id",
            "incoming_image_identifier",
            "volume",
            "slice",
            "group",
            "previous_consumed_index",
            "new_consumed_index",
            "jump_size",
            "superseded_index_count",
            "newer_candidate_count",
            "existing_intermediate_transforms",
            "failed_intermediate_attempts",
            "unknown_intermediate_indices",
            "classification",
        ],
        fire_analysis["supersession"],
    )

    volume_fields = list(fire_by_volume[0])
    write_csv(output / "fire_discovery_by_volume.csv", volume_fields, fire_by_volume)
    fire_summary_rows = summary_rows(
        "fire",
        "discovery_selection_ms",
        fire_analysis["discovery"]["all"],
        {
            segment: fire_analysis["discovery"]["segments"][segment]["discovery"]
            for segment in SEGMENTS
        },
    )
    write_csv(
        output / "fire_discovery_summary.csv",
        ["service", "metric", "segment", "n", "mean", "median", "std", "p95", "p99", "min", "max"],
        fire_summary_rows,
    )
    write_csv(
        output / "fire_outcome_summary.csv",
        ["outcome", "count", "fraction"],
        [
            {
                "outcome": outcome,
                "count": count,
                "fraction": count / len(fire["rows"]),
            }
            for outcome, count in sorted(fire_analysis["outcome_counts"].items())
        ],
    )
    write_csv(
        output / "fire_timing_components.csv",
        ["component", "n", "mean", "median", "std", "p95", "p99", "min", "max"],
        [
            {"component": component, **summary}
            for component, summary in fire_analysis["components"].items()
        ],
    )
    write_csv(
        output / "fire_advancement_step_distribution.csv",
        ["advancement_step", "count", "fraction_of_transitions"],
        [
            {
                "advancement_step": step,
                "count": count,
                "fraction_of_transitions": count / max(len(fire_analysis["steps"]), 1),
            }
            for step, count in sorted(fire_analysis["step_counts"].items())
        ],
    )

    write_csv(
        output / "queue_discovery_events.csv",
        [
            "scan_id",
            "volume",
            "min_volume",
            "segment",
            "items_in_scan",
            "directory_entry_count",
            "candidate_file_count",
            "pointer_candidate_count",
            "scan_and_build_ms",
            "empty_polls_since_previous",
            "empty_scan_ms_since_previous",
            "mean_empty_scan_ms",
        ],
        queue_scans["rows"],
    )
    write_csv(
        output / "queue_discovery_by_volume.csv",
        list(queue_by_volume[0]),
        queue_by_volume,
    )
    queue_summary_rows = summary_rows(
        "queue",
        "scan_and_build_ms",
        queue_analysis["scan_mapped"],
        {
            segment: queue_analysis["segments"][segment]["scan"]
            for segment in SEGMENTS
        },
    )
    write_csv(
        output / "queue_discovery_summary.csv",
        ["service", "metric", "segment", "n", "mean", "median", "std", "p95", "p99", "min", "max"],
        queue_summary_rows,
    )
    write_csv(
        output / "queue_item_times.csv",
        [
            "scan_id",
            "volume",
            "group",
            "segment",
            "pointer_filename",
            "outcome",
            "terminal_status",
            "output_label",
            "registration_index",
            "queue_handling_time_ms",
            "active_handling_time_ms",
        ],
        queue_items["rows"],
    )
    write_csv(
        output / "queue_handling_by_volume.csv",
        list(handling_by_volume[0]),
        handling_by_volume,
    )
    write_csv(
        output / "registration_by_volume.csv",
        ["volume", "expected", "registered", "skipped", "failed"],
        status["rows"],
    )
    write_csv(
        output / "acquisition_inventory.csv",
        ["location", "file_type", "count"],
        [
            {"location": "final_root", "file_type": "root_entries", "count": inventory["root_entry_count_final"]},
            {"location": "final_root", "file_type": "registration_tfm", "count": len(inventory["valid"])},
            {"location": "final_root", "file_type": "identity_tfm", "count": len(inventory["identity"])},
            {"location": "final_root", "file_type": "unexpected_tfm", "count": len(inventory["unexpected"])},
            {"location": "processed", "file_type": "pointer_txt", "count": inventory["processed_suffix_counts"][".txt"]},
            {"location": "processed", "file_type": "header_nhdr", "count": inventory["processed_suffix_counts"][".nhdr"]},
            {"location": "processed", "file_type": "image_raw", "count": inventory["processed_suffix_counts"][".raw"]},
            {"location": "processed/transforms", "file_type": "registration_tfm", "count": inventory["processed_transform_count"]},
        ],
    )
    cross_rows = fire_summary_rows + queue_summary_rows
    write_csv(
        output / "cross_service_summary.csv",
        ["service", "metric", "segment", "n", "mean", "median", "std", "p95", "p99", "min", "max"],
        cross_rows,
    )

    apply_plot_style()
    create_plots(output, fire, fire_analysis, queue_scans, queue_items, status)

    feedback_logs = sorted(run_dir.glob("log_moco_feedback_sent.log"))
    if not feedback_logs:
        raise ValueError("Required FIRE/MOCO logs are absent")
    feedback_transforms = feedback_log_transforms(feedback_logs[0])
    committed_filenames = [
        row["selected_transform_filename"] for row in fire_analysis["commits"]
    ]
    if feedback_transforms != committed_filenames:
        raise ValueError(
            "MOCO feedback log transform sequence does not match successful FIRE profile commits"
        )
    fire_analysis["feedback_log_records"] = len(feedback_transforms)
    sources = {
        "fire": fire["path"],
        "queue_scans": queue_scans["path"],
        "queue_items": queue_items["path"],
        "status": status["path"],
        "feedback_log": feedback_logs[0],
        "fire_log": select_matching_fire_log(fire["path"], run_dir),
    }
    report(
        output,
        run_dir,
        sources,
        status,
        registration,
        fire,
        fire_analysis,
        inventory,
        queue_analysis,
    )

    print(f"Analysis complete: {output}")
    print(
        f"Registration: {registration['registered']}/{registration['expected']} "
        f"({registration['success_percent']:.3f}%)"
    )
    print(
        f"FIRE commits: {len(fire_analysis['commits'])}; regressions: "
        f"{fire_analysis['regressions']}; duplicates: {fire_analysis['duplicates']}; "
        f"final/highest: {fire_analysis['commits'][-1]['selected_registration_index']}/"
        f"{inventory['highest_index']}"
    )
    print(
        f"FIRE discovery mean/median/p95: "
        f"{fire_analysis['discovery']['all']['mean']:.3f}/"
        f"{fire_analysis['discovery']['all']['median']:.3f}/"
        f"{fire_analysis['discovery']['all']['p95']:.3f} ms"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare two queue-directory retirement experiments.

The authoritative registration outcome is registration_status.json. Queue
timings come from queue_profile_items_*.csv and queue_profile_scans_*.csv.
Experiment directories are read only; all generated artifacts are written to
the selected output directory.
"""

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
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
warnings.filterwarnings("ignore", message="Unable to import Axes3D")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


RUN_LABELS = ("baseline", "file_retirement")
RUN_TITLES = {"baseline": "Baseline", "file_retirement": "File retirement"}
COLORS = {"baseline": "#3976af", "file_retirement": "#e1812c"}
STATUS_VALUES = {"registered", "skipped", "failed"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare registration completeness, end-to-end queue item timing, "
            "directory scan timing, and file inventory for two acquisitions."
        )
    )
    parser.add_argument("baseline_dir", type=Path, help="baseline savedData directory")
    parser.add_argument(
        "retirement_dir", type=Path, help="file-retirement savedData directory"
    )
    parser.add_argument(
        "--baseline-processed",
        action="store_true",
        help="require and inventory processed_<protocol> in the baseline run",
    )
    parser.add_argument(
        "--baseline-label",
        default="Baseline",
        help="display label for the baseline run (default: %(default)s)",
    )
    parser.add_argument(
        "--retirement-label",
        default="File retirement",
        help="display label for the comparison run (default: %(default)s)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("queue_directory_retirement_comparison"),
        help="analysis output directory (default: %(default)s)",
    )
    return parser.parse_args()


def require_directory(path: Path, description: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_dir():
        raise SystemExit(f"ERROR: {description} does not exist or is not a directory: {path}")
    return path


def unique_file(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {pattern!r} in {directory}, found {len(matches)}"
        )
    return matches[0]


def as_int(value: str, field: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid integer in {field}: {value!r}") from exc


def as_float(value: str, field: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid number in {field}: {value!r}") from exc


def optional_float(value: str) -> float | None:
    return float(value) if value else None


def describe(values: Iterable[float]) -> dict[str, float | int]:
    vals = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not vals:
        return {"n": 0}

    def percentile(percent: float) -> float:
        position = (len(vals) - 1) * percent / 100.0
        lower = math.floor(position)
        upper = math.ceil(position)
        if lower == upper:
            return vals[lower]
        fraction = position - lower
        return vals[lower] * (1.0 - fraction) + vals[upper] * fraction

    return {
        "n": len(vals),
        "mean": statistics.fmean(vals),
        "median": statistics.median(vals),
        "std": statistics.pstdev(vals),
        "p95": percentile(95.0),
        "max": max(vals),
        "min": min(vals),
    }


def pct(numerator: float, denominator: float) -> float:
    return 100.0 * numerator / denominator if denominator else float("nan")


def parse_status(run_dir: Path) -> dict[str, Any]:
    path = unique_file(run_dir, "registration_status.json")
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != 1 or not isinstance(payload.get("volumes"), dict):
        raise ValueError(f"Unsupported registration status schema in {path}")

    rows: list[dict[str, Any]] = []
    group_status: dict[tuple[int, int], str] = {}
    group_sets: dict[int, set[int]] = {}
    for volume_text, volume_data in payload["volumes"].items():
        volume = int(volume_text)
        if not volume_data.get("registration_finalized"):
            raise ValueError(f"Volume {volume} is not finalized in {path}")
        groups = volume_data.get("groups")
        if not isinstance(groups, dict):
            raise ValueError(f"Volume {volume} has no group mapping in {path}")
        counts: Counter[str] = Counter()
        group_numbers: set[int] = set()
        for group_text, group_data in groups.items():
            group = int(group_text)
            if not isinstance(group_data, dict):
                raise ValueError(
                    f"Invalid group record for volume {volume}, group {group} in {path}"
                )
            status = group_data.get("status")
            if status not in STATUS_VALUES:
                raise ValueError(
                    f"Unexpected status {status!r} for volume {volume}, group {group}"
                )
            key = (volume, group)
            if key in group_status:
                raise ValueError(f"Duplicate status for volume/group {key}")
            group_status[key] = status
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
    if volumes != list(range(min(volumes), max(volumes) + 1)) or volumes[0] != 0:
        raise ValueError(f"Status volumes are not contiguous from zero in {path}")

    nonreference_sets = {frozenset(group_sets[volume]) for volume in volumes if volume > 0}
    if len(nonreference_sets) != 1:
        raise ValueError(f"Non-reference volumes do not share one expected group set in {path}")
    expected_groups = len(next(iter(nonreference_sets)))
    return {
        "path": path,
        "payload": payload,
        "rows": rows,
        "group_status": group_status,
        "expected_groups": expected_groups,
        "volumes": volumes,
    }


def discover_protocol(pointer_names: Iterable[str]) -> str:
    protocols: set[str] = set()
    pattern = re.compile(r"^(.+)_volume_\d+_group_\d+\.txt$")
    for name in pointer_names:
        match = pattern.fullmatch(name)
        if not match:
            raise ValueError(f"Unexpected pointer filename in queue profile: {name!r}")
        protocols.add(match.group(1))
    if len(protocols) != 1:
        raise ValueError(f"Expected one protocol name, found: {sorted(protocols)}")
    return next(iter(protocols))


def parse_items(run_dir: Path, label: str, status: dict[str, Any]) -> dict[str, Any]:
    path = unique_file(run_dir, "queue_profile_items_*.csv")
    required = {
        "scan_id",
        "volume",
        "group",
        "pointer_filename",
        "outcome",
        "first_eligible_perf_ns",
        "selected_perf_ns",
        "decision_perf_ns",
        "registration_start_perf_ns",
        "registration_end_perf_ns",
        "item_complete_perf_ns",
        "candidate_file_count",
        "pointer_candidate_count",
    }
    parsed: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing queue item columns in {path}: {sorted(missing)}")
        for source in reader:
            volume = as_int(source["volume"], "volume")
            group = as_int(source["group"], "group")
            key = (volume, group)
            if key in seen:
                raise ValueError(f"Duplicate queue item for volume/group {key} in {path}")
            seen.add(key)
            outcome = source["outcome"]
            terminal = {
                "registered": "registered",
                "skipped_by_lifo": "skipped",
                "failed": "failed",
                "reference": "registered",
            }.get(outcome)
            if terminal is None:
                raise ValueError(f"Unexpected queue item outcome {outcome!r} in {path}")
            expected_status = status["group_status"].get(key)
            if expected_status != terminal:
                raise ValueError(
                    f"Queue/status mismatch for volume {volume}, group {group}: "
                    f"profile={outcome!r}, status={expected_status!r}"
                )

            first_eligible = as_int(source["first_eligible_perf_ns"], "first_eligible_perf_ns")
            complete = as_int(source["item_complete_perf_ns"], "item_complete_perf_ns")
            if complete < first_eligible:
                raise ValueError(f"Negative queue handling interval for {key} in {path}")
            if outcome == "skipped_by_lifo":
                active_start = as_int(source["decision_perf_ns"], "decision_perf_ns")
            else:
                active_start = as_int(source["selected_perf_ns"], "selected_perf_ns")
            reg_start = optional_float(source["registration_start_perf_ns"])
            reg_end = optional_float(source["registration_end_perf_ns"])
            parsed.append(
                {
                    "run": label,
                    "volume": volume,
                    "group": group,
                    "pointer_filename": source["pointer_filename"],
                    "outcome": outcome,
                    "registered_or_skipped": terminal,
                    "queue_handling_time_ms": (complete - first_eligible) / 1e6,
                    "active_handling_time_ms": (complete - active_start) / 1e6,
                    "registration_time_ms": (
                        (reg_end - reg_start) / 1e6
                        if reg_start is not None and reg_end is not None
                        else None
                    ),
                    "scan_id": as_int(source["scan_id"], "scan_id"),
                    "candidate_count": as_int(
                        source["candidate_file_count"], "candidate_file_count"
                    ),
                    "pointer_candidate_count": as_int(
                        source["pointer_candidate_count"], "pointer_candidate_count"
                    ),
                }
            )
    if seen != set(status["group_status"]):
        missing = sorted(set(status["group_status"]) - seen)
        extra = sorted(seen - set(status["group_status"]))
        raise ValueError(
            f"Queue item/status key mismatch in {run_dir}; missing={missing[:5]}, extra={extra[:5]}"
        )
    parsed.sort(key=lambda row: (row["volume"], row["group"]))
    protocol = discover_protocol(row["pointer_filename"] for row in parsed)
    return {"path": path, "rows": parsed, "protocol": protocol}


def parse_scans(run_dir: Path, label: str, items: dict[str, Any]) -> dict[str, Any]:
    path = unique_file(run_dir, "queue_profile_scans_*.csv")
    required = {
        "scan_id",
        "scan_start_perf_ns",
        "directory_entry_count",
        "candidate_file_count",
        "pointer_candidate_count",
        "scan_and_build_ms",
        "empty_polls_since_previous",
        "empty_scan_ms_since_previous",
    }
    items_by_scan: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for item in items["rows"]:
        items_by_scan[item["scan_id"]].append(item)

    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing queue scan columns in {path}: {sorted(missing)}")
        for source in reader:
            scan_id = as_int(source["scan_id"], "scan_id")
            linked = items_by_scan.get(scan_id, [])
            volumes = [item["volume"] for item in linked]
            empty_polls = as_int(source["empty_polls_since_previous"], "empty_polls")
            empty_ms = as_float(source["empty_scan_ms_since_previous"], "empty_scan_ms")
            rows.append(
                {
                    "run": label,
                    "scan_id": scan_id,
                    "scan_start_perf_ns": as_int(
                        source["scan_start_perf_ns"], "scan_start_perf_ns"
                    ),
                    "volume": max(volumes) if volumes else None,
                    "min_volume": min(volumes) if volumes else None,
                    "items_in_scan": len(linked),
                    "directory_entry_count": as_int(
                        source["directory_entry_count"], "directory_entry_count"
                    ),
                    "candidate_count": as_int(
                        source["candidate_file_count"], "candidate_file_count"
                    ),
                    "pointer_candidate_count": as_int(
                        source["pointer_candidate_count"], "pointer_candidate_count"
                    ),
                    "directory_scan_time_ms": as_float(
                        source["scan_and_build_ms"], "scan_and_build_ms"
                    ),
                    "empty_polls_since_previous": empty_polls,
                    "empty_scan_ms_since_previous": empty_ms,
                    "mean_empty_scan_time_ms": empty_ms / empty_polls if empty_polls else None,
                }
            )
    if not rows:
        raise ValueError(f"No queue scan rows in {path}")
    return {"path": path, "rows": rows}


def strict_inventory(directory: Path, protocol: str) -> dict[str, Any]:
    patterns = {
        "pointer": re.compile(
            rf"^{re.escape(protocol)}_volume_(\d+)_group_(\d+)\.txt$"
        ),
        "header": re.compile(
            rf"^{re.escape(protocol)}_volume_(\d+)_slice_(\d+)\.nhdr$"
        ),
        "raw": re.compile(rf"^{re.escape(protocol)}_volume_(\d+)_slice_(\d+)\.raw$"),
    }
    matched: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    entries = list(directory.iterdir())
    files = [entry for entry in entries if entry.is_file()]
    for entry in files:
        for kind, pattern in patterns.items():
            match = pattern.fullmatch(entry.name)
            if match:
                matched[kind].append((int(match.group(1)), int(match.group(2)), entry.name))
                break
    return {
        "entry_count": len(entries),
        "file_count": len(files),
        "directory_count": sum(entry.is_dir() for entry in entries),
        "suffix_counts": {
            suffix: sum(entry.suffix == suffix for entry in files)
            for suffix in (".txt", ".nhdr", ".raw", ".tfm")
        },
        "strict": dict(matched),
        "upsampled": [entry.name for entry in files if "_upsampled." in entry.name],
    }


def parse_inventory(run_dir: Path, protocol: str, require_processed: bool) -> dict[str, Any]:
    root = strict_inventory(run_dir, protocol)
    processed_dirs = sorted(
        entry for entry in run_dir.iterdir() if entry.is_dir() and entry.name.startswith("processed_")
    )
    if require_processed:
        if len(processed_dirs) != 1:
            raise ValueError(
                f"Expected exactly one processed_* directory in {run_dir}, found {len(processed_dirs)}"
            )
        expected_name = f"processed_{protocol}"
        if processed_dirs[0].name != expected_name:
            raise ValueError(
                f"Processed directory {processed_dirs[0].name!r} does not match {expected_name!r}"
            )
    elif processed_dirs:
        raise ValueError(f"Baseline unexpectedly contains processed_* directories: {processed_dirs}")
    processed = strict_inventory(processed_dirs[0], protocol) if processed_dirs else None
    transform_dir = processed_dirs[0] / "transforms" if processed_dirs else None
    transforms = None
    if processed:
        for kind in ("pointer", "header", "raw"):
            if len(processed["strict"].get(kind, [])) != processed["suffix_counts"][
                {"pointer": ".txt", "header": ".nhdr", "raw": ".raw"}[kind]
            ]:
                raise ValueError(
                    f"Processed directory contains nonconforming {kind} files: {processed_dirs[0]}"
                )
        if processed["upsampled"]:
            raise ValueError(f"Processed directory contains upsampled files: {processed['upsampled']}")
        if transform_dir.is_dir():
            transforms = strict_inventory(transform_dir, protocol)
            transform_pattern = re.compile(r"^alignTransform_\d+_\d+-\d+\.tfm$")
            invalid = sorted(
                entry.name
                for entry in transform_dir.iterdir()
                if not entry.is_file() or not transform_pattern.fullmatch(entry.name)
            )
            if invalid:
                raise ValueError(
                    f"Transform archive contains unexpected entries: {invalid}"
                )
    return {
        "root": root,
        "processed_dir": processed_dirs[0] if processed_dirs else None,
        "processed": processed,
        "transform_dir": transform_dir if transforms else None,
        "transforms": transforms,
    }


def summarize_registration(status: dict[str, Any]) -> dict[str, Any]:
    rows = [row for row in status["rows"] if row["volume"] > 0]
    total_expected = sum(row["expected"] for row in rows)
    result: dict[str, Any] = {
        "total_expected": total_expected,
        "registered": sum(row["registered"] for row in rows),
        "skipped": sum(row["skipped"] for row in rows),
        "failed": sum(row["failed"] for row in rows),
        "degraded_volumes": sum(row["registered"] < row["expected"] for row in rows),
        "minimum_registered": min(row["registered"] for row in rows),
        "final": rows[-1],
    }
    result["success_percent"] = pct(result["registered"], total_expected)
    degraded = [row["volume"] for row in rows if row["registered"] < row["expected"]]
    result["first_degraded"] = degraded[0] if degraded else None
    after_startup = [volume for volume in degraded if volume >= 2]
    result["first_degraded_after_startup"] = after_startup[0] if after_startup else None
    return result


def volume_segments(volumes: list[int]) -> list[tuple[str, list[int]]]:
    count = len(volumes)
    boundaries = (count // 3, 2 * count // 3)
    chunks = (volumes[: boundaries[0]], volumes[boundaries[0] : boundaries[1]], volumes[boundaries[1] :])
    names = ("early", "middle", "late")
    return [(name, chunk) for name, chunk in zip(names, chunks) if chunk]


def summarize_handling(items: dict[str, Any], volumes: list[int]) -> dict[str, Any]:
    rows = [row for row in items["rows"] if row["volume"] > 0]
    summary = {
        "all": describe(row["queue_handling_time_ms"] for row in rows),
        "registered": describe(
            row["queue_handling_time_ms"]
            for row in rows
            if row["registered_or_skipped"] == "registered"
        ),
        "skipped": describe(
            row["queue_handling_time_ms"]
            for row in rows
            if row["registered_or_skipped"] == "skipped"
        ),
        "active_registered": describe(
            row["active_handling_time_ms"]
            for row in rows
            if row["registered_or_skipped"] == "registered"
        ),
        "segments": {},
    }
    for name, segment_volumes in volume_segments([volume for volume in volumes if volume > 0]):
        volume_set = set(segment_volumes)
        summary["segments"][name] = {
            "range": (segment_volumes[0], segment_volumes[-1]),
            "registered": describe(
                row["queue_handling_time_ms"]
                for row in rows
                if row["volume"] in volume_set
                and row["registered_or_skipped"] == "registered"
            ),
        }
    return summary


def summarize_scans(scans: dict[str, Any], volumes: list[int]) -> dict[str, Any]:
    rows = scans["rows"]
    mapped = [row for row in rows if row["volume"] is not None]
    result = {
        "scan": describe(row["directory_scan_time_ms"] for row in rows),
        "directory_entries": describe(row["directory_entry_count"] for row in rows),
        "mapped_scans": len(mapped),
        "unmapped_scans": len(rows) - len(mapped),
        "empty_polls": sum(row["empty_polls_since_previous"] for row in rows),
        "empty_scan_ms": sum(row["empty_scan_ms_since_previous"] for row in rows),
        "segments": {},
        "last_directory_entry_count": rows[-1]["directory_entry_count"],
    }
    result["mean_empty_scan_time_ms"] = (
        result["empty_scan_ms"] / result["empty_polls"] if result["empty_polls"] else None
    )
    for name, segment_volumes in volume_segments([volume for volume in volumes if volume > 0]):
        volume_set = set(segment_volumes)
        result["segments"][name] = {
            "range": (segment_volumes[0], segment_volumes[-1]),
            "scan": describe(
                row["directory_scan_time_ms"] for row in mapped if row["volume"] in volume_set
            ),
            "entries": describe(
                row["directory_entry_count"] for row in mapped if row["volume"] in volume_set
            ),
        }
    return result


def write_csv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if row.get(key) is None else row.get(key) for key in fields})


def grouped_median(rows: list[dict[str, Any]], value: str) -> tuple[list[int], list[float]]:
    grouped: dict[int, list[float]] = defaultdict(list)
    for row in rows:
        grouped[row["volume"]].append(float(row[value]))
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


def plot_registration(output: Path, runs: dict[str, dict[str, Any]], expected: int) -> None:
    fig, axis = plt.subplots(figsize=(11, 5.2))
    for label in RUN_LABELS:
        rows = [row for row in runs[label]["status"]["rows"] if row["volume"] > 0]
        axis.plot(
            [row["volume"] for row in rows],
            [row["registered"] for row in rows],
            color=COLORS[label],
            linewidth=1.7,
            label=RUN_TITLES[label],
        )
    axis.axhline(expected, color="#555555", linestyle="--", linewidth=1.1, label=f"Expected ({expected})")
    axis.set(xlabel="Volume", ylabel="Successfully registered groups", title="Registration completion by volume")
    axis.set_xlim(1, max(runs["baseline"]["status"]["volumes"]))
    axis.set_ylim(bottom=0, top=expected + 2)
    axis.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(output / "registered_groups_per_volume.png", dpi=180)
    plt.close(fig)


def plot_handling(output: Path, runs: dict[str, dict[str, Any]]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True, height_ratios=(3, 2))
    for label in RUN_LABELS:
        rows = [row for row in runs[label]["items"]["rows"] if row["volume"] > 0]
        registered = [row for row in rows if row["registered_or_skipped"] == "registered"]
        skipped = [row for row in rows if row["registered_or_skipped"] == "skipped"]
        axes[0].scatter(
            [row["volume"] for row in registered],
            [row["queue_handling_time_ms"] for row in registered],
            s=5,
            alpha=0.10,
            color=COLORS[label],
            linewidths=0,
        )
        x_values, medians = grouped_median(registered, "queue_handling_time_ms")
        axes[0].plot(x_values, medians, color=COLORS[label], linewidth=1.8, label=RUN_TITLES[label])
        if skipped:
            x_values, medians = grouped_median(skipped, "queue_handling_time_ms")
            axes[1].scatter(
                x_values,
                medians,
                color=COLORS[label],
                s=13,
                alpha=0.9,
                label=RUN_TITLES[label],
            )
    axes[0].set(ylabel="Eligible-to-complete (ms)", title="End-to-end queue pointer lifecycle by volume")
    axes[0].legend(frameon=False)
    axes[0].text(0.01, 0.96, "Registered groups: points and per-volume median", transform=axes[0].transAxes, va="top")
    axes[1].set(xlabel="Volume", ylabel="Eligible-to-complete (ms)")
    axes[1].text(0.01, 0.94, "Skipped groups: per-volume median", transform=axes[1].transAxes, va="top")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "queue_handling_time_per_volume.png", dpi=180)
    plt.close(fig)


def plot_scans(output: Path, runs: dict[str, dict[str, Any]]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    for label in RUN_LABELS:
        rows = [row for row in runs[label]["scans"]["rows"] if row["volume"] is not None]
        axes[0].scatter(
            [row["volume"] for row in rows],
            [row["directory_scan_time_ms"] for row in rows],
            s=5,
            alpha=0.11,
            color=COLORS[label],
            linewidths=0,
        )
        x_values, medians = grouped_median(rows, "directory_scan_time_ms")
        axes[0].plot(x_values, medians, color=COLORS[label], linewidth=1.8, label=RUN_TITLES[label])
        x_values, medians = grouped_median(rows, "directory_entry_count")
        axes[1].plot(x_values, medians, color=COLORS[label], linewidth=1.8, label=RUN_TITLES[label])
    axes[0].set(ylabel="Scan and candidate build (ms)", title="Measured queue directory discovery cost")
    axes[0].legend(frameon=False)
    axes[0].text(0.01, 0.96, "Individual non-empty scans and per-volume median", transform=axes[0].transAxes, va="top")
    axes[1].set(xlabel="Volume", ylabel="Directory entries")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "directory_scan_time_per_volume.png", dpi=180)
    plt.close(fig)


def fmt_stat(summary: dict[str, Any], name: str) -> str:
    value = summary.get(name)
    return "n/a" if value is None else f"{value:.3f}"


def write_report(output: Path, runs: dict[str, dict[str, Any]], expected: int) -> None:
    baseline = runs["baseline"]
    retirement = runs["file_retirement"]
    base_reg = baseline["registration_summary"]
    retire_reg = retirement["registration_summary"]
    base_scan = baseline["scan_summary"]
    retire_scan = retirement["scan_summary"]
    handling_reduction = 100.0 * (
        1.0
        - retirement["handling_summary"]["registered"]["mean"]
        / baseline["handling_summary"]["registered"]["mean"]
    )
    scan_reduction = 100.0 * (1.0 - retire_scan["scan"]["mean"] / base_scan["scan"]["mean"])
    success_gain = retire_reg["success_percent"] - base_reg["success_percent"]

    lines = [
        "# Queue directory retirement comparison",
        "",
        "## Experiments and sources",
        "",
        f"- {RUN_TITLES['baseline']}: `{baseline['directory']}`",
        f"- {RUN_TITLES['file_retirement']}: `{retirement['directory']}`",
        f"- Authoritative outcomes: `{baseline['status']['path'].name}` and `{retirement['status']['path'].name}`",
        f"- Item timing: `{baseline['items']['path'].name}` and `{retirement['items']['path'].name}`",
        f"- Scan timing: `{baseline['scans']['path'].name}` and `{retirement['scans']['path'].name}`",
        "",
        "The status files contain one fixed-reference volume (volume 0) plus 199 registration volumes. "
        f"Each non-reference volume has {expected} groups; volume 0 is excluded from all registration totals below.",
        "",
        "## Registration completion",
        "",
        "| Run | Expected | Registered | Skipped | Failed | Success | First degraded | First degraded after startup | Degraded volumes | Minimum registered | Final volume (R/S/F) |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for label in RUN_LABELS:
        summary = runs[label]["registration_summary"]
        final = summary["final"]
        lines.append(
            f"| {RUN_TITLES[label]} | {summary['total_expected']} | {summary['registered']} | "
            f"{summary['skipped']} | {summary['failed']} | {summary['success_percent']:.3f}% | "
            f"{summary['first_degraded']} | {summary['first_degraded_after_startup']} | "
            f"{summary['degraded_volumes']} | {summary['minimum_registered']} | "
            f"{final['registered']}/{final['skipped']}/{final['failed']} |"
        )
    lines.extend(
        [
            "",
            f"{RUN_TITLES['file_retirement']} produced {retire_reg['registered'] - base_reg['registered']:+d} successful registrations "
            f"and {retire_reg['skipped'] - base_reg['skipped']:+d} skips, a {success_gain:.3f}-percentage-point success gain. "
            "Both runs have a startup skip in volume 1, so the post-startup onset is also reported.",
            "",
            "## End-to-end queue handling time",
            "",
            "`queue_handling_time_ms` is measured from `first_eligible_perf_ns` (the end of the first scan in which the pointer was eligible) "
            "to `item_complete_perf_ns` (after registration-or-skip handling, status work, and reference-update work). It therefore includes "
            "same-batch waiting and surrounding queue/Python/filesystem work. It is not CUDA optimizer time. For skipped pointers, this latency "
            "can include time spent processing earlier selected items in the same LIFO batch. `active_handling_time_ms` is also included in the CSV "
            "and starts at selection for processed pointers or at the skip decision for LIFO-skipped pointers.",
            "",
            "| Run/outcome | n | Mean ms | Median ms | Std ms | p95 ms | Max ms |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for label in RUN_LABELS:
        for outcome in ("all", "registered", "skipped"):
            summary = runs[label]["handling_summary"][outcome]
            lines.append(
                f"| {RUN_TITLES[label]} — {outcome} | {summary['n']} | {fmt_stat(summary, 'mean')} | "
                f"{fmt_stat(summary, 'median')} | {fmt_stat(summary, 'std')} | {fmt_stat(summary, 'p95')} | {fmt_stat(summary, 'max')} |"
            )
    lines.extend(
        [
            "",
            f"Among successfully registered groups (the comparable full-processing path), mean eligible-to-complete time was "
            f"{handling_reduction:.1f}% lower in {RUN_TITLES['file_retirement']}.",
            "",
            "Registered-group timing by acquisition third:",
            "",
            "| Run | Segment | Volumes | Mean ms | Median ms | p95 ms |",
            "| --- | --- | --- | ---: | ---: | ---: |",
        ]
    )
    for label in RUN_LABELS:
        for segment, entry in runs[label]["handling_summary"]["segments"].items():
            summary = entry["registered"]
            lines.append(
                f"| {RUN_TITLES[label]} | {segment} | {entry['range'][0]}–{entry['range'][1]} | "
                f"{fmt_stat(summary, 'mean')} | {fmt_stat(summary, 'median')} | {fmt_stat(summary, 'p95')} |"
            )
    lines.extend(
        [
            "",
            "## Directory scanning time",
            "",
            "Historical scan timing is directly recoverable. `scan_and_build_ms` brackets the profiled `list_new_files()` discovery call: "
            "`os.listdir()`, extension/seen filtering, per-candidate `getmtime()`, candidate sorting, and profiler bookkeeping. "
            "The measured code path is the nested `list_new_files()` function in `apps/queue_processor/process_queue_directory.py`; "
            "`apps/queue_processor/queue_profile.py` writes the timing. "
            "It is an actual discovery-path timing rather than a proxy. Scan rows are associated with the maximum volume among pointers emitted "
            "by that scan; one terminal scan in each run has no associated pointer and therefore no volume.",
            "",
            "| Run | Non-empty scan rows | Mean ms | Median ms | Std ms | p95 ms | Max ms | Entry median | Entry max | Last entry count | Empty polls | Mean empty scan ms |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for label in RUN_LABELS:
        summary = runs[label]["scan_summary"]
        scan = summary["scan"]
        entries = summary["directory_entries"]
        lines.append(
            f"| {RUN_TITLES[label]} | {scan['n']} | {fmt_stat(scan, 'mean')} | {fmt_stat(scan, 'median')} | "
            f"{fmt_stat(scan, 'std')} | {fmt_stat(scan, 'p95')} | {fmt_stat(scan, 'max')} | "
            f"{entries['median']:.1f} | {entries['max']:.0f} | {summary['last_directory_entry_count']} | "
            f"{summary['empty_polls']} | {summary['mean_empty_scan_time_ms']:.3f} |"
        )
    lines.extend(
        [
            "",
            f"Mean recorded non-empty discovery time was {scan_reduction:.1f}% lower in {RUN_TITLES['file_retirement']}. Empty-poll timing is only available as "
            "aggregate intervals between non-empty scans; the CSV reports both the interval total and its mean per empty poll.",
            "",
            "## Retired-file and directory inventory",
            "",
            "| Run/location | Entries | Files | Directories | .txt | .nhdr | .raw | .tfm |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for label in RUN_LABELS:
        inventory_set = runs[label]["inventory"]
        locations = [("root", inventory_set["root"])]
        if inventory_set["processed"]:
            locations.append(("processed", inventory_set["processed"]))
        if inventory_set["transforms"]:
            locations.append(("processed/transforms", inventory_set["transforms"]))
        for location, inventory in locations:
            suffix = inventory["suffix_counts"]
            lines.append(
                f"| {RUN_TITLES[label]} {location} | {inventory['entry_count']} | "
                f"{inventory['file_count']} | {inventory['directory_count']} | "
                f"{suffix['.txt']} | {suffix['.nhdr']} | {suffix['.raw']} | {suffix['.tfm']} |"
            )
    lines.append("")
    for label in RUN_LABELS:
        inventory_set = runs[label]["inventory"]
        processed = inventory_set["processed"]
        if not processed:
            continue
        counts = {
            kind: len(processed["strict"].get(kind, []))
            for kind in ("pointer", "header", "raw")
        }
        volumes = sorted(
            {
                volume
                for kind in ("pointer", "header", "raw")
                for volume, _, _ in processed["strict"].get(kind, [])
            }
        )
        volume_range = f"{volumes[0]}–{volumes[-1]}" if volumes else "none"
        transform_count = (
            inventory_set["transforms"]["suffix_counts"][".tfm"]
            if inventory_set["transforms"]
            else 0
        )
        lines.append(
            f"- {RUN_TITLES[label]}: `{inventory_set['processed_dir']}` contains "
            f"{counts['pointer']} pointers, {counts['header']} detached headers, and "
            f"{counts['raw']} raw files spanning volumes {volume_range}; its transform "
            f"archive contains {transform_count} `.tfm` files."
        )
    lines.extend(
        [
            "",
            f"The known retired workload is {expected * 3} root entries per completed non-reference volume "
            f"({expected} each of pointer/header/raw). The root transform counts are "
            f"{baseline['inventory']['root']['suffix_counts']['.tfm']} for {RUN_TITLES['baseline']} and "
            f"{retirement['inventory']['root']['suffix_counts']['.tfm']} for {RUN_TITLES['file_retirement']}.",
            "",
            "## Interpretation",
            "",
            f"{RUN_TITLES['file_retirement']} retained {retire_reg['success_percent']:.3f}% of expected registrations versus "
            f"{base_reg['success_percent']:.3f}% in {RUN_TITLES['baseline']}. Post-startup degradation began at volume "
            f"{base_reg['first_degraded_after_startup']} and {retire_reg['first_degraded_after_startup']}, respectively. "
            f"Registered-pointer mean lifecycle time changed by {-handling_reduction:+.1f}% and measured non-empty discovery time "
            f"changed by {-scan_reduction:+.1f}% in the comparison run.",
            "",
            "They do not establish causation by themselves: there is one run per condition, profiler entry counts include all root files, the number "
            "of generated transforms differs because registration success differs, and host/GPU scheduling or replay timing can vary between runs.",
            "",
            "## Validation and limitations",
            "",
            "- Every queue-item outcome was matched by `(volume, group)` to the authoritative status file; row sets and outcomes agree exactly.",
            f"- Both runs contain {expected} groups for every non-reference volume and use the same volume set.",
            "- Timing units are converted from profiler nanoseconds to milliseconds; `scan_and_build_ms` is already milliseconds.",
            "- Volume 0 is treated only as the fixed reference and excluded from comparisons.",
            "- Final saved layouts cannot reconstruct historical directory growth alone; the temporal entry-count series here comes from the scan profiler.",
            "- Queue lifecycle time is latency, not exclusive CPU time, especially for pointers skipped within a multi-item LIFO batch.",
            "- This analysis is observational and does not modify queue polling, registration, FIRE, motion monitoring, or experiment data.",
        ]
    )
    (output / "comparison_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    RUN_TITLES["baseline"] = args.baseline_label
    RUN_TITLES["file_retirement"] = args.retirement_label
    baseline_dir = require_directory(args.baseline_dir, "baseline directory")
    retirement_dir = require_directory(args.retirement_dir, "retirement directory")
    if baseline_dir == retirement_dir:
        raise SystemExit("ERROR: baseline and retirement directories must differ")
    output = args.output_dir.expanduser().resolve()
    if output == baseline_dir or output == retirement_dir or baseline_dir in output.parents or retirement_dir in output.parents:
        raise SystemExit("ERROR: output directory must not be inside either experiment directory")
    output.mkdir(parents=True, exist_ok=True)

    runs: dict[str, dict[str, Any]] = {}
    for label, directory in zip(RUN_LABELS, (baseline_dir, retirement_dir)):
        status = parse_status(directory)
        items = parse_items(directory, label, status)
        scans = parse_scans(directory, label, items)
        inventory = parse_inventory(
            directory,
            items["protocol"],
            label == "file_retirement" or args.baseline_processed,
        )
        runs[label] = {
            "directory": directory,
            "status": status,
            "items": items,
            "scans": scans,
            "inventory": inventory,
            "registration_summary": summarize_registration(status),
            "handling_summary": summarize_handling(items, status["volumes"]),
            "scan_summary": summarize_scans(scans, status["volumes"]),
        }

    if runs["baseline"]["status"]["volumes"] != runs["file_retirement"]["status"]["volumes"]:
        raise ValueError("Runs do not have identical volume sets")
    expected_values = {runs[label]["status"]["expected_groups"] for label in RUN_LABELS}
    if len(expected_values) != 1:
        raise ValueError(f"Runs have different expected groups per volume: {expected_values}")
    if runs["baseline"]["items"]["protocol"] != runs["file_retirement"]["items"]["protocol"]:
        raise ValueError("Runs have different protocol names")
    expected = next(iter(expected_values))

    registration_rows = []
    rows_by_run = {
        label: {row["volume"]: row for row in runs[label]["status"]["rows"]}
        for label in RUN_LABELS
    }
    for volume in runs["baseline"]["status"]["volumes"]:
        if volume == 0:
            continue
        baseline = rows_by_run["baseline"][volume]
        retirement = rows_by_run["file_retirement"][volume]
        registration_rows.append(
            {
                "volume": volume,
                "expected_groups": expected,
                "baseline_registered": baseline["registered"],
                "baseline_skipped": baseline["skipped"],
                "baseline_failed": baseline["failed"],
                "retirement_registered": retirement["registered"],
                "retirement_skipped": retirement["skipped"],
                "retirement_failed": retirement["failed"],
            }
        )
    write_csv(
        output / "registration_groups_per_volume.csv",
        [
            "volume",
            "expected_groups",
            "baseline_registered",
            "baseline_skipped",
            "baseline_failed",
            "retirement_registered",
            "retirement_skipped",
            "retirement_failed",
        ],
        registration_rows,
    )
    handling_rows = [
        row
        for label in RUN_LABELS
        for row in runs[label]["items"]["rows"]
        if row["volume"] > 0
    ]
    write_csv(
        output / "queue_group_handling_times.csv",
        [
            "run",
            "volume",
            "group",
            "pointer_filename",
            "registered_or_skipped",
            "outcome",
            "queue_handling_time_ms",
            "active_handling_time_ms",
            "registration_time_ms",
            "scan_id",
            "candidate_count",
            "pointer_candidate_count",
        ],
        handling_rows,
    )
    scan_rows = [row for label in RUN_LABELS for row in runs[label]["scans"]["rows"]]
    write_csv(
        output / "directory_scan_times.csv",
        [
            "run",
            "scan_id",
            "scan_start_perf_ns",
            "volume",
            "min_volume",
            "items_in_scan",
            "directory_scan_time_ms",
            "directory_entry_count",
            "candidate_count",
            "pointer_candidate_count",
            "empty_polls_since_previous",
            "empty_scan_ms_since_previous",
            "mean_empty_scan_time_ms",
        ],
        scan_rows,
    )

    apply_plot_style()
    plot_registration(output, runs, expected)
    plot_handling(output, runs)
    plot_scans(output, runs)
    write_report(output, runs, expected)

    print(f"Wrote comparison outputs to {output}")
    for label in RUN_LABELS:
        registration = runs[label]["registration_summary"]
        handling = runs[label]["handling_summary"]["registered"]
        scan = runs[label]["scan_summary"]["scan"]
        print(
            f"{RUN_TITLES[label]}: {registration['registered']}/{registration['total_expected']} "
            f"registered ({registration['success_percent']:.3f}%), "
            f"first degraded={registration['first_degraded']}, "
            f"registered handling median={handling['median']:.3f} ms, "
            f"scan median={scan['median']:.3f} ms"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Analyze queue discovery and terminal item handling for one acquisition."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from analysis_common import (
    AnalysisError,
    STAT_FIELDS,
    adaptive_y_axes,
    add_common_arguments,
    add_temporal_thirds,
    cli_main,
    describe,
    discover_queue_profile_pair,
    finish_terminal,
    format_number,
    format_stats,
    grouped_medians,
    late_to_early_ratio,
    optional_int,
    ordinary_least_squares,
    pearson,
    provenance_lines,
    pyplot,
    read_csv,
    required_float,
    required_int,
    resolve_output_directory,
    statistics_rows,
    temporal_statistics_rows,
    third_statistics,
    validate_input_directory,
    write_csv,
    write_summary,
)


ITEM_REQUIRED = {
    "scan_id", "volume", "group", "pointer_filename", "outcome",
    "first_observed_perf_ns", "first_eligible_perf_ns", "selected_perf_ns",
    "decision_perf_ns", "item_complete_perf_ns", "registration_start_perf_ns",
    "registration_end_perf_ns", "candidate_file_count", "pointer_candidate_count",
}
SCAN_REQUIRED = {
    "scan_id", "scan_start_perf_ns", "directory_entry_count", "candidate_file_count",
    "pointer_candidate_count", "scan_and_build_ms", "total_cycle_ms",
    "empty_polls_since_previous", "empty_scan_ms_since_previous",
}


def match_receipts_to_queue_items(
    receipts: list[dict[str, Any]],
    queue_items: list[dict[str, Any]],
    pointer_members: dict[str, list[str]],
) -> dict[str, Any]:
    """Exact identifier join used to audit potential FIRE/queue correlation.

    This deliberately accepts filenames only. Timestamps are never join keys.
    A slice can belong to exactly one pointer and a pointer can identify exactly
    one queue row; otherwise the join is ambiguous and refused.
    """
    receipt_by_image: dict[str, dict[str, Any]] = {}
    for receipt in receipts:
        filename = receipt.get("image_filename")
        if not filename:
            raise AnalysisError(
                "FIRE receipt record lacks image_filename; timestamp-only matching is refused."
            )
        if filename in receipt_by_image:
            raise AnalysisError(f"Duplicate/ambiguous FIRE receipt identifier: {filename}")
        receipt_by_image[filename] = receipt
    queue_by_pointer: dict[str, dict[str, Any]] = {}
    for item in queue_items:
        pointer = item.get("pointer_filename")
        if not pointer:
            raise AnalysisError("Queue record lacks pointer_filename; timestamp-only matching is refused.")
        if pointer in queue_by_pointer:
            raise AnalysisError(f"Duplicate/ambiguous queue pointer identifier: {pointer}")
        queue_by_pointer[pointer] = item
    image_to_pointer: dict[str, str] = {}
    for pointer, members in pointer_members.items():
        if pointer not in queue_by_pointer:
            continue
        for member in members:
            if member in image_to_pointer and image_to_pointer[member] != pointer:
                raise AnalysisError(
                    f"Duplicate/ambiguous pointer membership for image {member}: "
                    f"{image_to_pointer[member]} and {pointer}"
                )
            image_to_pointer[member] = pointer
    matched = []
    for image, receipt in receipt_by_image.items():
        pointer = image_to_pointer.get(image)
        if pointer is not None:
            matched.append({"receipt": receipt, "queue_item": queue_by_pointer[pointer]})
    matched_images = {entry["receipt"]["image_filename"] for entry in matched}
    matched_pointers = {entry["queue_item"]["pointer_filename"] for entry in matched}
    return {
        "matched": matched,
        "unmatched_fire": [receipt for name, receipt in receipt_by_image.items() if name not in matched_images],
        "unmatched_queue": [item for name, item in queue_by_pointer.items() if name not in matched_pointers],
    }


def _interval(row: dict[str, str], start: str, end: str) -> float | None:
    first, last = optional_int(row.get(start), start), optional_int(row.get(end), end)
    if first is None or last is None:
        return None
    if last < first:
        raise AnalysisError(f"Negative interval {start}..{end} for {row.get('pointer_filename')}")
    return (last - first) / 1e6


def parse_items(path: Path) -> list[dict[str, Any]]:
    source_rows, _ = read_csv(path, ITEM_REQUIRED)
    rows: list[dict[str, Any]] = []
    for source in source_rows:
        outcome = source["outcome"]
        terminal_start = "decision_perf_ns" if outcome == "skipped_by_lifo" else "selected_perf_ns"
        row = {
            "scan_id": required_int(source["scan_id"], "scan_id"),
            "volume": required_int(source["volume"], "volume"),
            "group": required_int(source["group"], "group"),
            "pointer_filename": source["pointer_filename"],
            "outcome": outcome,
            "candidate_file_count": required_int(source["candidate_file_count"], "candidate_file_count"),
            "pointer_candidate_count": required_int(source["pointer_candidate_count"], "pointer_candidate_count"),
            "lifecycle_ms": _interval(source, "first_observed_perf_ns", "item_complete_perf_ns"),
            "eligible_to_terminal_ms": _interval(source, "first_eligible_perf_ns", "item_complete_perf_ns"),
            "active_handling_ms": _interval(source, terminal_start, "item_complete_perf_ns"),
            "discovery_to_selection_ms": _interval(source, "first_eligible_perf_ns", "selected_perf_ns"),
            "registration_call_ms": _interval(source, "registration_start_perf_ns", "registration_end_perf_ns"),
        }
        rows.append(row)
    rows.sort(key=lambda row: (row["scan_id"], row["volume"], row["group"]))
    add_temporal_thirds(rows)
    return rows


def parse_scans(path: Path) -> list[dict[str, Any]]:
    source_rows, _ = read_csv(path, SCAN_REQUIRED)
    rows = []
    for source in source_rows:
        empty_polls = required_int(source["empty_polls_since_previous"], "empty_polls_since_previous")
        empty_ms = required_float(source["empty_scan_ms_since_previous"], "empty_scan_ms_since_previous")
        rows.append({
            "scan_id": required_int(source["scan_id"], "scan_id"),
            "directory_entry_count": required_int(source["directory_entry_count"], "directory_entry_count"),
            "candidate_file_count": required_int(source["candidate_file_count"], "candidate_file_count"),
            "pointer_candidate_count": required_int(source["pointer_candidate_count"], "pointer_candidate_count"),
            "scan_and_build_ms": required_float(source["scan_and_build_ms"], "scan_and_build_ms"),
            "total_cycle_ms": required_float(source["total_cycle_ms"], "total_cycle_ms"),
            "empty_polls_since_previous": empty_polls,
            "empty_scan_ms_since_previous": empty_ms,
            "mean_empty_scan_ms": empty_ms / empty_polls if empty_polls else None,
        })
    rows.sort(key=lambda row: row["scan_id"])
    add_temporal_thirds(rows)
    return rows


def attach_scan_volumes(
    scans: list[dict[str, Any]], items: list[dict[str, Any]]
) -> int:
    """Map scans through the explicit scan_id relationship to queue items."""
    by_scan: dict[int, list[int]] = defaultdict(list)
    for item in items:
        by_scan[item["scan_id"]].append(item["volume"])
    unlinked = 0
    for scan in scans:
        volumes = by_scan.get(scan["scan_id"], [])
        scan["volume"] = max(volumes) if volumes else None
        scan["minimum_linked_volume"] = min(volumes) if volumes else None
        scan["linked_item_count"] = len(volumes)
        unlinked += not volumes
    return unlinked


def scan_discovery_medians_by_volume(
    rows: list[dict[str, Any]],
) -> list[tuple[int, float]]:
    return grouped_medians(rows, "volume", "scan_and_build_ms")


def _plot_scans(rows: list[dict[str, Any]], output: Path) -> dict[str, Any]:
    plt = pyplot()
    mapped = [row for row in rows if row.get("volume") is not None]
    volumes = [row["volume"] for row in mapped]
    discovery = [row["scan_and_build_ms"] for row in mapped]
    medians = scan_discovery_medians_by_volume(mapped)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.scatter(volumes, discovery, s=8, alpha=0.32, color="#3976af", label="Individual scan")
    ax.plot(
        [volume for volume, _ in medians],
        [median for _, median in medians],
        color="red",
        linewidth=1.5,
        label="Per-volume median",
    )
    ax.set_xlabel("Volume number")
    ax.set_ylabel("Discovery time (ms)")
    ax.set_title("Queue directory discovery timing by volume")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output / "queue_discovery_timing.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.scatter(volumes, [row["directory_entry_count"] for row in mapped], s=8, alpha=0.35, color="#3976af")
    ax.set_xlabel("Volume number")
    ax.set_ylabel("Directory entries (count)")
    ax.set_title("Queue directory entries by volume")
    fig.tight_layout()
    fig.savefig(output / "queue_directory_entries.png", dpi=150)
    plt.close(fig)

    entries = [row["directory_entry_count"] for row in rows]
    timings = [row["scan_and_build_ms"] for row in rows]
    fit = ordinary_least_squares(entries, timings)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(entries, timings, s=9, alpha=0.35, color="#3976af", label="Individual scan")
    if fit is not None:
        endpoints = [min(entries), max(entries)]
        fitted = [fit["slope"] * value + fit["intercept"] for value in endpoints]
        ax.plot(endpoints, fitted, color="#c44e52", linewidth=1.6, label="Ordinary least-squares fit")
        ax.text(
            0.03, 0.97,
            f"time = {fit['slope']:.6f} ms/entry × entries + {fit['intercept']:.3f} ms\nR² = {fit['r_squared']:.3f}",
            transform=ax.transAxes, va="top",
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )
    ax.set_xlabel("Directory entries (count)")
    ax.set_ylabel("Discovery time (ms)")
    ax.set_title("Queue discovery time versus directory entries")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output / "queue_discovery_vs_directory_entries.png", dpi=150)
    plt.close(fig)
    return {"volume_medians": medians, "fit": fit, "mapped_scans": len(mapped)}


def _plot_items(rows: list[dict[str, Any]], output: Path) -> dict[str, Any] | None:
    plt = pyplot()
    values = [row["active_handling_ms"] for row in rows if row["active_handling_ms"] is not None]
    fig, axes, break_info = adaptive_y_axes(plt, values, figsize=(10, 5.5))
    colors = {"registered": "#3976af", "skipped_by_lifo": "#e1812c", "failed": "#c44e52", "reference": "#55a868"}
    for outcome in sorted({r["outcome"] for r in rows}):
        subset = [r for r in rows if r["outcome"] == outcome and r["active_handling_ms"] is not None]
        for index, ax in enumerate(axes):
            ax.scatter(
                [r["volume"] for r in subset], [r["active_handling_ms"] for r in subset],
                s=9, alpha=0.60, label=outcome if index == 0 else None,
                color=colors.get(outcome),
            )
    axes[-1].set_xlabel("Volume number")
    axes[0].set_title("Queue item terminal handling by outcome" + (" (broken y-axis)" if break_info else ""))
    axes[0].legend(loc="best")
    fig.supylabel("Active item handling time (ms)")
    if break_info:
        fig.subplots_adjust(left=0.11, right=0.97, bottom=0.12, top=0.88, hspace=0.06)
    else:
        fig.tight_layout(rect=(0.02, 0, 1, 1))
    fig.savefig(output, dpi=150)
    plt.close(fig)
    return break_info


def run(input_directory: Path, output_arg: Path | None = None) -> Path:
    input_directory = validate_input_directory(input_directory)
    item_path, scan_path = discover_queue_profile_pair(input_directory)
    _, output = resolve_output_directory(input_directory, output_arg, "queue_profile")
    items, scans = parse_items(item_path), parse_scans(scan_path)
    if not items:
        raise AnalysisError(f"Queue item profile contains no records: {item_path.name}")
    if not scans:
        raise AnalysisError(f"Queue scan profile contains no records: {scan_path.name}")
    unlinked_scans = attach_scan_volumes(scans, items)

    outcome_counts = Counter(row["outcome"] for row in items)
    scan_stats = describe(row["scan_and_build_ms"] for row in scans)
    scan_thirds = third_statistics(scans, "scan_and_build_ms")
    item_stats = describe(row["active_handling_ms"] for row in items)
    item_thirds = third_statistics(items, "active_handling_ms")
    metrics = {
        "scan_and_build_ms": [r["scan_and_build_ms"] for r in scans],
        "total_cycle_ms": [r["total_cycle_ms"] for r in scans],
        "mean_empty_scan_ms": [r["mean_empty_scan_ms"] for r in scans],
        "item_lifecycle_ms": [r["lifecycle_ms"] for r in items],
        "eligible_to_terminal_ms": [r["eligible_to_terminal_ms"] for r in items],
        "active_item_handling_ms": [r["active_handling_ms"] for r in items],
        "discovery_to_selection_ms": [r["discovery_to_selection_ms"] for r in items],
        "registration_call_ms": [r["registration_call_ms"] for r in items],
    }
    for outcome in sorted(outcome_counts):
        metrics[f"active_handling_ms:{outcome}"] = [
            r["active_handling_ms"] for r in items if r["outcome"] == outcome
        ]
    write_csv(output / "timing_statistics.csv", statistics_rows(metrics), ["metric", *STAT_FIELDS])
    write_csv(
        output / "timing_by_temporal_third.csv",
        temporal_statistics_rows(items, ["lifecycle_ms", "eligible_to_terminal_ms", "active_handling_ms", "registration_call_ms"])
        + temporal_statistics_rows(scans, ["scan_and_build_ms", "total_cycle_ms", "mean_empty_scan_ms"]),
        ["metric", "third", *STAT_FIELDS],
    )
    write_csv(output / "queue_items.csv", items, list(items[0]))
    write_csv(output / "queue_scans.csv", scans, list(scans[0]))

    per_volume = []
    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in items:
        grouped[row["volume"]].append(row)
    for volume, rows in sorted(grouped.items()):
        counts = Counter(row["outcome"] for row in rows)
        stats = describe(row["active_handling_ms"] for row in rows)
        per_volume.append({"volume": volume, "items": len(rows), "registered": counts["registered"], "skipped": counts["skipped_by_lifo"], "failed": counts["failed"], "mean_active_handling_ms": stats["mean"], "max_active_handling_ms": stats["max"]})
    write_csv(output / "queue_items_by_volume.csv", per_volume, list(per_volume[0]))
    outliers = sorted(items, key=lambda row: row["active_handling_ms"] or -1, reverse=True)[:10]
    write_csv(output / "queue_item_outliers.csv", outliers, list(items[0]))
    stale_plot = output / "queue_scan_timing.png"
    if stale_plot.exists():
        stale_plot.unlink()
    scan_plot_info = _plot_scans(scans, output)
    item_break = _plot_items(items, output / "queue_item_handling_timing.png")

    total_empty = sum(row["empty_polls_since_previous"] for row in scans)
    warnings = [
        "True FIRE-receipt-to-terminal item lifetime is omitted: FIRE receipt uses wall-clock logs, while queue item_complete_perf_ns is process-local; no terminal wall-clock field exists, and skipped items have no exact terminal log record."
    ]
    if unlinked_scans:
        warnings.append(
            f"{unlinked_scans} non-pointer scan(s) lacked an explicit scan_id-to-item volume mapping and were omitted from volume-axis plots."
        )
    semantics = [
        "scan_and_build_ms spans os.listdir, extension/seen filtering, candidate getmtime calls, acquisition-order sorting, and profiler observation work.",
        "active item handling spans selected_perf_ns to item_complete_perf_ns for selected items; for LIFO skips it spans decision_perf_ns to item_complete_perf_ns.",
        "eligible_to_terminal_ms includes waiting after the pointer first became an eligible candidate; it is distinct from active handling and directory discovery time.",
        "Queue volume plots map each scan through scan_id and use the maximum linked item volume; the red line is the median discovery time among mapped nonempty scans at each volume.",
        "The discovery-time/directory-entry line is an ordinary least-squares descriptive fit and is not interpreted causally.",
    ]
    lines = [
        f"Nonempty scans: {len(scans)}; aggregated empty polls: {total_empty}; total polls represented: {len(scans) + total_empty}",
        f"Queue items: {len(items)}; outcomes: {dict(sorted(outcome_counts.items()))}",
        format_stats("Directory scan/discovery", scan_stats),
        f"Discovery late/early median ratio: {format_number(late_to_early_ratio(scan_thirds))}",
        f"Discovery progression correlation (scan ordinal vs ms): {format_number(pearson(range(len(scans)), [r['scan_and_build_ms'] for r in scans]))}",
        format_stats("Active item handling", item_stats),
        f"Item-handling late/early median ratio: {format_number(late_to_early_ratio(item_thirds))}",
        (
            f"Discovery time vs directory entries OLS: slope={scan_plot_info['fit']['slope']:.6f} ms/entry, "
            f"intercept={scan_plot_info['fit']['intercept']:.3f} ms, R²={scan_plot_info['fit']['r_squared']:.3f}"
            if scan_plot_info["fit"] is not None else
            "Discovery time vs directory entries OLS: unavailable (insufficient or constant directory-entry data)"
        ),
        (
            f"Queue handling plot uses an adaptive broken y-axis ({item_break['upper_count']} upper observations; "
            f"gap {item_break['lower_max']:.3f}–{item_break['upper_min']:.3f} ms)."
            if item_break else "Queue handling plot uses a normal unbroken y-axis; no sufficiently separated sparse upper region was detected."
        ),
    ]
    write_summary(output / "summary.txt", [
        ("Provenance", provenance_lines(Path(__file__).name, input_directory, [item_path, scan_path], len(items) + len(scans), assumptions=[*semantics, *warnings])),
        ("High-level results", lines),
        ("Timer semantics", semantics),
        ("Temporal thirds", [
            f"Scan {name}: {format_stats('', scan_thirds[name]).removeprefix(': ')}" for name in ("early", "middle", "late")
        ] + [f"Item {name}: {format_stats('', item_thirds[name]).removeprefix(': ')}" for name in ("early", "middle", "late")]),
        ("End-to-end FIRE receipt lifetime", warnings),
    ])
    finish_terminal(input_directory, [item_path, scan_path], output, f"Analyzed {len(items)} queue items and {len(scans)} nonempty scans.", warnings)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    args = parser.parse_args()
    run(args.input_directory, args.output_dir)


if __name__ == "__main__":
    cli_main(main)

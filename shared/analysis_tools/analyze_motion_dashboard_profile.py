#!/usr/bin/env python3
"""Analyze motion-dashboard generation profiling for one acquisition."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import Any

from analysis_common import (
    AnalysisError, STAT_FIELDS, add_common_arguments, add_temporal_thirds, cli_main,
    describe, discover_unique, finish_terminal, format_number, format_stats,
    late_to_early_ratio, optional_float, provenance_lines, pyplot, read_csv,
    required_int, resolve_output_directory, statistics_rows, third_statistics,
    temporal_statistics_rows,
    validate_input_directory, write_csv, write_summary,
)


REQUIRED = {"dashboard_index", "outcome", "total_ms", "figure_create_ms", "savefig_ms", "atomic_replace_ms"}
KNOWN_TIMINGS = (
    "trigger_setup_ms", "monitor_data_prep_ms", "tsnr_prepare_ms", "plot_data_extract_ms",
    "figure_create_ms", "plot_population_ms", "annotation_ms", "tsnr_panel_ms", "layout_ms",
    "savefig_ms", "cleanup_ms", "plot_total_ms", "output_validate_ms", "atomic_replace_ms",
    "published_read_ms", "stream_push_ms", "component_sum_ms", "unattributed_ms", "total_ms",
)
PROFILE_TIME_FIELDS = (
    "profile_start_perf_ns", "generation_start_perf_ns", "timestamp_ns", "wall_time_ns",
)


def parse_profile(path: Path) -> tuple[list[dict[str, Any]], list[str]]:
    source_rows, fields = read_csv(path, REQUIRED)
    timing_fields = [field for field in KNOWN_TIMINGS if field in fields]
    rows = []
    for source in source_rows:
        row: dict[str, Any] = {
            "dashboard_index": required_int(source["dashboard_index"], "dashboard_index"),
            "trigger_transform_index": required_int(source["trigger_transform_index"], "trigger_transform_index") if source.get("trigger_transform_index") else None,
            "volume": required_int(source["volume"], "volume") if source.get("volume") else None,
            "group": required_int(source["group"], "group") if source.get("group") else None,
            "outcome": source["outcome"],
        }
        for field in timing_fields:
            row[field] = optional_float(source.get(field), field)
        for field in PROFILE_TIME_FIELDS:
            if field in fields:
                row[field] = optional_float(source.get(field), field)
        rows.append(row)
    rows.sort(key=lambda row: row["dashboard_index"])
    add_temporal_thirds(rows)
    return rows, timing_fields


def dashboard_x_axis(rows: list[dict[str, Any]]) -> tuple[list[float], str, str]:
    """Use structured volume, or an explicit profile timestamp as elapsed time."""
    if rows and all(row.get("volume") is not None for row in rows):
        return [float(row["volume"]) for row in rows], "Volume number", "volume"
    for field in PROFILE_TIME_FIELDS:
        if rows and all(row.get(field) is not None for row in rows):
            start = min(float(row[field]) for row in rows)
            return (
                [(float(row[field]) - start) / 1e9 for row in rows],
                "Run time (s)",
                f"elapsed time from first {field}",
            )
    raise AnalysisError(
        "Dashboard rows lack complete volume identifiers and an explicit profile timestamp; "
        "a volume or run-time axis cannot be established reliably."
    )


def _plots(rows: list[dict[str, Any]], timing_fields: list[str], output: Path) -> str:
    plt = pyplot()
    x_values, x_label, mapping = dashboard_x_axis(rows)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.scatter(x_values, [r["total_ms"] for r in rows], s=22, alpha=0.75, color="#3976af")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Total generation time (ms)")
    ax.set_title(f"Motion-dashboard generation timing by {x_label.lower()}")
    fig.tight_layout()
    fig.savefig(output / "dashboard_generation_timing.png", dpi=150)
    plt.close(fig)

    components = [field for field in timing_fields if field not in {"total_ms", "plot_total_ms", "component_sum_ms", "unattributed_ms"}]
    means = [(field, describe(r.get(field) for r in rows)["mean"]) for field in components]
    means = [(field, value) for field, value in means if value is not None]
    if means:
        fig, ax = plt.subplots(figsize=(10, max(4, len(means) * 0.28)))
        ax.barh([name.removesuffix("_ms") for name, _ in means], [value for _, value in means])
        ax.set_xlabel("Mean time (ms)")
        ax.set_title("Mean dashboard timing components")
        fig.tight_layout()
        fig.savefig(output / "dashboard_component_means.png", dpi=150)
        plt.close(fig)
    return mapping


def run(input_directory: Path, output_arg: Path | None = None) -> Path:
    input_directory = validate_input_directory(input_directory)
    profile_path = discover_unique(
        input_directory, "motion_dashboard_profile_*.csv", "motion-dashboard profile",
        missing_hint="MOTION_DASHBOARD_PROFILE_FLAG may have been disabled.",
    )
    _, output = resolve_output_directory(input_directory, output_arg, "motion_dashboard_profile")
    rows, timing_fields = parse_profile(profile_path)
    if not rows:
        raise AnalysisError(f"Motion-dashboard profile contains no records: {profile_path.name}")
    metrics = {field: [row.get(field) for row in rows] for field in timing_fields}
    write_csv(output / "dashboard_timing_statistics.csv", statistics_rows(metrics), ["metric", *STAT_FIELDS])
    write_csv(output / "dashboard_timing_by_temporal_third.csv", temporal_statistics_rows(rows, timing_fields), ["metric", "third", *STAT_FIELDS])
    write_csv(output / "dashboard_generations.csv", rows, list(rows[0]))
    x_mapping = _plots(rows, timing_fields, output)
    total_stats = describe(row["total_ms"] for row in rows)
    thirds = third_statistics(rows, "total_ms")
    savefig = describe(row.get("savefig_ms") for row in rows)
    fraction = 100 * float(savefig["mean"]) / float(total_stats["mean"]) if savefig["mean"] is not None and total_stats["mean"] else None
    limitations = [
        "Generation cadence and dashboard duty cycle are unavailable because the profiler records no wall-clock or monotonic generation timestamp.",
        "Backlog/recovery cannot be inferred directly from the current profile schema.",
        "plot_total_ms overlaps detailed plotting phases and is not treated as an additive component.",
        "total_ms spans the complete plot_motion_data attempt; component_sum_ms contains the non-overlapping phase fields defined by dashboard_profile.py.",
    ]
    results = [
        f"Dashboard generation attempts: {len(rows)}; outcomes: {dict(sorted(Counter(r['outcome'] for r in rows).items()))}",
        format_stats("Total generation", total_stats),
        format_stats("savefig", savefig),
        f"Mean savefig share of mean total: {format_number(fraction)}%",
        f"Total-time late/early median ratio: {format_number(late_to_early_ratio(thirds))}",
        f"Dashboard timing x-axis mapping: {x_mapping}",
    ]
    write_summary(output / "summary.txt", [
        ("Provenance", provenance_lines(Path(__file__).name, input_directory, [profile_path], len(rows), assumptions=limitations)),
        ("High-level results", results),
        ("Timer semantics", limitations[2:]),
        ("Temporal thirds", [f"{name}: {format_stats('total', thirds[name])}" for name in ("early", "middle", "late")]),
        ("Unavailable metrics", limitations[:2]),
    ])
    finish_terminal(input_directory, [profile_path], output, f"Analyzed {len(rows)} dashboard generation attempts.", limitations[:2])
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    args = parser.parse_args()
    run(args.input_directory, args.output_dir)


if __name__ == "__main__":
    cli_main(main)

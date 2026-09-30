#!/usr/bin/env python3
"""Analyze FIRE transform discovery, selection, and feedback profiling."""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path
from typing import Any

from analysis_common import (
    AnalysisError, STAT_FIELDS, add_common_arguments, add_temporal_thirds, cli_main,
    describe, discover_unique, finish_terminal, format_number, format_stats,
    grouped_medians, late_to_early_ratio, optional_float, optional_int, pearson, provenance_lines,
    pyplot, read_csv, read_json, required_float, required_int,
    resolve_output_directory, statistics_rows, third_statistics,
    temporal_statistics_rows,
    validate_input_directory, write_csv, write_summary,
)


REQUIRED = {
    "lookup_id", "incoming_image_identifier", "volume", "slice", "group",
    "root_entry_count", "valid_transform_count", "newer_candidate_count",
    "selected_transform_filename", "selected_registration_index",
    "last_consumed_registration_index_before", "last_consumed_registration_index_after",
    "discovery_selection_ms", "transform_read_ms", "conversion_ms", "packaging_ms",
    "conversion_package_ms", "feedback_log_ms", "feedback_send_ms", "outcome",
}
TIMINGS = ("discovery_selection_ms", "transform_read_ms", "conversion_ms", "packaging_ms", "conversion_package_ms", "feedback_log_ms", "feedback_send_ms")
TRANSFORM_RE = re.compile(r"^alignTransform_(\d+)_\d+-\d+(?:_identity)?\.tfm$")


def parse_profile(path: Path) -> list[dict[str, Any]]:
    source_rows, _ = read_csv(path, REQUIRED)
    rows = []
    for source in source_rows:
        row: dict[str, Any] = {
            "lookup_id": required_int(source["lookup_id"], "lookup_id"),
            "incoming_image_identifier": required_int(source["incoming_image_identifier"], "incoming_image_identifier"),
            "volume": required_int(source["volume"], "volume"),
            "slice": required_int(source["slice"], "slice"),
            "group": required_int(source["group"], "group"),
            "root_entry_count": required_int(source["root_entry_count"], "root_entry_count") if source["root_entry_count"] else None,
            "valid_transform_count": required_int(source["valid_transform_count"], "valid_transform_count") if source["valid_transform_count"] else None,
            "newer_candidate_count": required_int(source["newer_candidate_count"], "newer_candidate_count") if source["newer_candidate_count"] else None,
            "selected_transform_filename": source["selected_transform_filename"],
            "selected_registration_index": optional_int(source["selected_registration_index"], "selected_registration_index"),
            "consumed_before": optional_int(source["last_consumed_registration_index_before"], "last_consumed_registration_index_before"),
            "consumed_after": optional_int(source["last_consumed_registration_index_after"], "last_consumed_registration_index_after"),
            "outcome": source["outcome"],
        }
        for field in TIMINGS:
            row[field] = optional_float(source[field], field)
        rows.append(row)
    rows.sort(key=lambda row: row["lookup_id"])
    add_temporal_thirds(rows)
    return rows


def transform_inventory(directory: Path) -> tuple[int | None, list[str]]:
    by_index: dict[int, set[str]] = {}
    for path in directory.glob("*.tfm"):
        match = TRANSFORM_RE.fullmatch(path.name)
        if match:
            by_index.setdefault(int(match.group(1)), set()).add(path.name)
    for path in directory.glob("processed_*/transforms/*.tfm"):
        match = TRANSFORM_RE.fullmatch(path.name)
        if match:
            by_index.setdefault(int(match.group(1)), set()).add(path.name)
    ambiguous = [f"index {index}: {sorted(names)}" for index, names in by_index.items() if len(names) > 1]
    return (max(by_index) if by_index and not ambiguous else None), ambiguous


def discovery_medians_by_volume(rows: list[dict[str, Any]]) -> list[tuple[int, float]]:
    return grouped_medians(rows, "volume", "discovery_selection_ms")


def _plots(rows: list[dict[str, Any]], output: Path) -> dict[str, Any]:
    plt = pyplot()
    stale_plot = output / "fire_moco_lookup_timing.png"
    if stale_plot.exists():
        stale_plot.unlink()
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.scatter(
        [r["volume"] for r in rows],
        [r["discovery_selection_ms"] for r in rows],
        s=8, alpha=0.35, color="blue", label="Individual lookup",
    )
    medians = discovery_medians_by_volume(rows)
    ax.plot(
        [volume for volume, _ in medians],
        [median for _, median in medians],
        color="red", linewidth=1.5, label="Per-volume median",
    )
    ax.set_xlabel("Volume number")
    ax.set_ylabel("Discovery/selection time (ms)")
    ax.set_title("FIRE MoCo discovery/selection timing by volume")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output / "fire_moco_discovery_timing.png", dpi=150)
    plt.close(fig)
    return {"x_field": "volume", "medians": medians}


def run(input_directory: Path, output_arg: Path | None = None) -> Path:
    input_directory = validate_input_directory(input_directory)
    profile_path = discover_unique(input_directory, "fire_moco_profile_*.csv", "FIRE MoCo profile", missing_hint="FIRE MoCo profiling may have been disabled or MoCo may have been off.")
    state_candidates = sorted(input_directory.glob("fire_moco_state.moco"))
    state_path = state_candidates[0] if len(state_candidates) == 1 else None
    _, output = resolve_output_directory(input_directory, output_arg, "fire_moco_profile")
    rows = parse_profile(profile_path)
    if not rows:
        raise AnalysisError(f"FIRE MoCo profile contains no records: {profile_path.name}")
    outcomes = Counter(row["outcome"] for row in rows)
    successful = [row for row in rows if row["outcome"] == "feedback_sent"]
    changes = []
    for row in successful:
        before, after = row["consumed_before"], row["consumed_after"]
        if after is not None:
            changes.append(after - before if before is not None else None)
    single = sum(change == 1 for change in changes)
    jumps = sum(change is not None and change > 1 for change in changes)
    regressions = sum(change is not None and change < 0 for change in changes)
    selected_indices = [row["selected_registration_index"] for row in rows if row["selected_registration_index"] is not None]
    duplicate_selections = len(selected_indices) - len(set(selected_indices))
    count_fields = ("root_entry_count", "valid_transform_count", "newer_candidate_count")
    timings = {field: [row[field] for row in rows] for field in (*TIMINGS, *count_fields)}
    write_csv(output / "fire_moco_timing_statistics.csv", statistics_rows(timings), ["metric", *STAT_FIELDS])
    write_csv(output / "fire_moco_metrics_by_temporal_third.csv", temporal_statistics_rows(rows, [*TIMINGS, *count_fields]), ["metric", "third", *STAT_FIELDS])
    write_csv(output / "fire_moco_lookups.csv", rows, list(rows[0]))
    _plots(rows, output)

    warnings: list[str] = []
    artifacts = [profile_path]
    state = None
    if state_path:
        artifacts.append(state_path)
        try:
            state = read_json(state_path)
            if not isinstance(state, dict) or state.get("schema_version") != 1:
                warnings.append("fire_moco_state.moco has an unsupported schema and was not used.")
                state = None
        except AnalysisError as error:
            warnings.append(f"fire_moco_state.moco could not be used: {error}")
    else:
        warnings.append("fire_moco_state.moco is unavailable; final committed state cannot be corroborated.")
    terminal = successful[-1]["consumed_after"] if successful else None
    disagreements = []
    if state is not None:
        state_index = optional_int(state.get("committed_registration_index"), "committed_registration_index")
        if state_index != terminal:
            disagreements.append(f"Profiler terminal consumed index {terminal} differs from .moco committed index {state_index}.")
    final_generated, inventory_ambiguity = transform_inventory(input_directory)
    if inventory_ambiguity:
        warnings.append("Final generated index is ambiguous: " + "; ".join(inventory_ambiguity))
    lag = final_generated - terminal if final_generated is not None and terminal is not None else None
    discovery = describe(row["discovery_selection_ms"] for row in rows)
    thirds = third_statistics(rows, "discovery_selection_ms")
    count_pairs = [(r["root_entry_count"], r["discovery_selection_ms"]) for r in rows if r["root_entry_count"] is not None]
    results = [
        f"Lookup attempts: {len(rows)}; outcomes: {dict(sorted(outcomes.items()))}",
        f"Successful transform advancements: {len(successful)}; no-new-transform lookups: {outcomes['no_new_transform']}",
        format_stats("Discovery/selection", discovery),
        f"Discovery late/early median ratio: {format_number(late_to_early_ratio(thirds))}",
        f"Root-entry-count/discovery correlation: {format_number(pearson([x for x, _ in count_pairs], [y for _, y in count_pairs]))}",
        f"Single-index advances: {single}; multi-index jumps: {jumps}; duplicate selections: {duplicate_selections}; regressions: {regressions}",
        f"First/final consumed registration index: {successful[0]['consumed_after'] if successful else 'unavailable'} / {terminal if terminal is not None else 'unavailable'}",
        f"Final generated index from unique transform inventory: {final_generated if final_generated is not None else 'unavailable'}; generated-versus-consumed lag: {lag if lag is not None else 'unavailable'}",
    ]
    semantics = [
        "discovery_selection_ms wraps transform-directory enumeration, filename validation, newer-index filtering, and highest-index selection.",
        "transform_read_ms wraps SimpleITK transform reading; conversion_ms and packaging_ms are separately scoped, while conversion_package_ms spans both and therefore overlaps them.",
        "feedback_log_ms covers append-to-log work; feedback_send_ms covers connection.send_feedback until return. Advancement commits only after send succeeds.",
        "Registration indices are compared numerically but are not assumed contiguous.",
    ]
    limitations = semantics + warnings
    write_summary(output / "summary.txt", [
        ("Provenance", provenance_lines(Path(__file__).name, input_directory, artifacts, len(rows), optional_unavailable=warnings, assumptions=limitations)),
        ("High-level results", results),
        ("Timer semantics", semantics),
        ("State disagreements", disagreements or ["none"]),
        ("Temporal thirds", [f"{name}: {format_stats('discovery', thirds[name])}" for name in ("early", "middle", "late")]),
    ])
    finish_terminal(input_directory, artifacts, output, f"Analyzed {len(rows)} FIRE MoCo lookups with {len(successful)} advancements.", warnings + disagreements)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    args = parser.parse_args()
    run(args.input_directory, args.output_dir)


if __name__ == "__main__":
    cli_main(main)

#!/usr/bin/env python3
"""Analyze authoritative registration outcomes from registration_status.json."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import Any

from analysis_common import (
    AnalysisError,
    add_common_arguments,
    cli_main,
    discover_unique,
    finish_terminal,
    format_number,
    provenance_lines,
    pyplot,
    read_json,
    resolve_output_directory,
    validate_input_directory,
    write_csv,
    write_summary,
)


VALID_STATUSES = {"registered", "skipped", "failed"}


def parse_status(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    document = read_json(path)
    if not isinstance(document, dict) or document.get("schema_version") != 1:
        raise AnalysisError(f"Unsupported or missing registration status schema_version in {path.name}")
    volumes = document.get("volumes")
    if not isinstance(volumes, dict):
        raise AnalysisError(f"Required field 'volumes' is missing or is not an object in {path.name}")
    volume_rows, exception_rows = [], []
    for volume_text, volume_data in volumes.items():
        try:
            volume = int(volume_text)
        except (TypeError, ValueError) as error:
            raise AnalysisError(f"Invalid volume key in {path.name}: {volume_text!r}") from error
        if not isinstance(volume_data, dict):
            raise AnalysisError(f"Volume {volume} is not an object in {path.name}")
        finalized = volume_data.get("registration_finalized") is True
        groups = volume_data.get("groups")
        if not isinstance(groups, dict):
            raise AnalysisError(f"Required groups object is missing for volume {volume}")
        counts: Counter[str] = Counter()
        for group_text, record in groups.items():
            try:
                group = int(group_text)
            except (TypeError, ValueError) as error:
                raise AnalysisError(f"Invalid group key {group_text!r} in volume {volume}") from error
            if not isinstance(record, dict) or record.get("status") not in VALID_STATUSES:
                raise AnalysisError(f"Invalid terminal status for volume {volume}, group {group}")
            status = record["status"]
            counts[status] += 1
            if status != "registered":
                exception_rows.append({"volume": volume, "group": group, "status": status, "reason": record.get("reason", "")})
        volume_rows.append({
            "volume": volume, "expected": len(groups), "registered": counts["registered"],
            "skipped": counts["skipped"], "failed": counts["failed"],
            "registration_finalized": finalized,
        })
    volume_rows.sort(key=lambda row: row["volume"])
    exception_rows.sort(key=lambda row: (row["volume"], row["group"]))
    return volume_rows, exception_rows


def constant_expected_count(rows: list[dict[str, Any]]) -> int | None:
    expected = {int(row["expected"]) for row in rows}
    return next(iter(expected)) if len(expected) == 1 else None


def build_completeness_figure(rows: list[dict[str, Any]]):
    plt = pyplot()
    from matplotlib.ticker import MaxNLocator, MultipleLocator

    fig, ax = plt.subplots(figsize=(11, 5))
    volumes = [r["volume"] for r in rows]
    registered = [r["registered"] for r in rows]
    skipped = [r["skipped"] for r in rows]
    failed = [r["failed"] for r in rows]
    ax.bar(volumes, registered, label="registered", color="#3976af", width=1.0)
    ax.bar(volumes, skipped, bottom=registered, label="skipped", color="#e1812c", width=1.0)
    bottom = [a + b for a, b in zip(registered, skipped)]
    ax.bar(volumes, failed, bottom=bottom, label="failed", color="#c44e52", width=1.0)
    ax.set_xlabel("Volume")
    ax.set_ylabel("Registration opportunities")
    ax.set_title("Registration completeness by volume")
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    legend = ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=3,
        borderaxespad=0.0,
    )
    expected = constant_expected_count(rows)
    percent_axis = None
    if expected is not None and expected > 0:
        percent_axis = ax.secondary_yaxis(
            "right",
            functions=(
                lambda count: count * 100.0 / expected,
                lambda percent: percent * expected / 100.0,
            ),
        )
        percent_axis.set_ylabel("Percent of volume (%)")
        percent_axis.yaxis.set_major_locator(MultipleLocator(20))
    fig.subplots_adjust(bottom=0.28, right=0.88 if percent_axis else 0.96)
    return fig, ax, percent_axis, expected, legend


def _plot(rows: list[dict[str, Any]], output: Path) -> dict[str, Any]:
    plt = pyplot()
    fig, _ax, percent_axis, expected, _legend = build_completeness_figure(rows)
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return {"percent_axis": percent_axis is not None, "expected": expected}


def run(input_directory: Path, output_arg: Path | None = None) -> Path:
    input_directory = validate_input_directory(input_directory)
    status_path = discover_unique(input_directory, "registration_status.json", "registration status JSON")
    _, output = resolve_output_directory(input_directory, output_arg, "registration_completeness")
    volumes, exceptions = parse_status(status_path)
    if not volumes:
        raise AnalysisError(f"registration_status.json contains no finalized volume records: {status_path}")
    totals = {key: sum(row[key] for row in volumes) for key in ("expected", "registered", "skipped", "failed")}
    success = 100.0 * totals["registered"] / totals["expected"] if totals["expected"] else None
    degraded = [row for row in volumes if row["registered"] < row["expected"]]
    unfinalized = [row["volume"] for row in volumes if not row["registration_finalized"]]
    write_csv(output / "registration_completeness_by_volume.csv", volumes, list(volumes[0]))
    write_csv(output / "skipped_registrations.csv", [r for r in exceptions if r["status"] == "skipped"], ["volume", "group", "status", "reason"])
    write_csv(output / "failed_registrations.csv", [r for r in exceptions if r["status"] == "failed"], ["volume", "group", "status", "reason"])
    plot_info = _plot(volumes, output / "registration_completeness_by_volume.png")
    limitations = [
        "Schema v1 publishes finalized volumes only; volumes absent from the JSON cannot be distinguished from never-expected or unfinalized volumes.",
        "Schema v1 does not identify a startup/reference region, so no early volumes are excluded from degradation metrics.",
        "Skip/failure reasons are unavailable unless a nonstandard record supplies a 'reason' field.",
    ]
    if not plot_info["percent_axis"]:
        limitations.append(
            "The right-side percent axis was omitted because expected registration opportunities vary by volume; no single count-to-percent conversion is valid."
        )
    final = volumes[-1]
    results = [
        f"Volumes represented: {len(volumes)} ({volumes[0]['volume']} through {volumes[-1]['volume']}; contiguity not assumed)",
        f"Expected opportunities: {totals['expected']}",
        f"Registered: {totals['registered']}; skipped: {totals['skipped']}; failed: {totals['failed']}",
        f"Registration success: {format_number(success)}%",
        f"Degraded volumes: {len(degraded)}; first degraded: {degraded[0]['volume'] if degraded else 'none'}",
        f"Minimum registered in any represented volume: {min(row['registered'] for row in volumes)}",
        f"Final represented volume {final['volume']}: {final['registered']} registered, {final['skipped']} skipped, {final['failed']} failed, finalized={final['registration_finalized']}",
        f"Explicit unfinalized volume records: {unfinalized or 'none'}",
        (
            f"Completeness percent axis: counts divided by constant expected count {plot_info['expected']} (0={0}%, full={100}%)."
            if plot_info["percent_axis"] else
            "Completeness percent axis: omitted because expected opportunities vary by volume."
        ),
    ]
    write_summary(output / "summary.txt", [
        ("Provenance", provenance_lines(Path(__file__).name, input_directory, [status_path], len(volumes), assumptions=limitations)),
        ("High-level results", results),
        ("Exact skipped locations", [f"({r['volume']}, {r['group']})" for r in exceptions if r["status"] == "skipped"] or ["none"]),
        ("Exact failed locations", [f"({r['volume']}, {r['group']})" for r in exceptions if r["status"] == "failed"] or ["none"]),
    ])
    finish_terminal(input_directory, [status_path], output, f"Registration success: {format_number(success)}% across {totals['expected']} opportunities.", limitations)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    args = parser.parse_args()
    run(args.input_directory, args.output_dir)


if __name__ == "__main__":
    cli_main(main)

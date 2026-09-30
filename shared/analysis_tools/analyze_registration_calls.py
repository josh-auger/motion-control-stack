#!/usr/bin/env python3
"""Analyze registration-call timing separately from queue lifecycle timing."""

from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from analysis_common import (
    AnalysisError, STAT_FIELDS, adaptive_y_axes, add_common_arguments, add_temporal_thirds, cli_main,
    describe, finish_terminal, format_number, format_stats,
    late_to_early_ratio, optional_int, provenance_lines, pyplot, read_csv,
    read_json, required_int, resolve_output_directory, statistics_rows, third_statistics,
    temporal_statistics_rows,
    validate_input_directory, write_csv, write_summary,
)


REQUIRED = {
    "volume", "group", "outcome", "output_label", "reg_engine", "cuda_execution_mode",
    "registration_start_perf_ns", "registration_end_perf_ns",
}
PROTOCOL_LINE = re.compile(r"SLIMM REGISTER\s+(\S+)\s+completed.*?SLIMM_PROTOCOL\s+SUCCESS\s+(.*)$")
LOAD_REFERENCE_LINE = re.compile(
    r"SLIMM LOAD_REFERENCE completed.*?SLIMM_PROTOCOL\s+SUCCESS\s+(.*)$"
)
TIMING_LINE = re.compile(r"\[TIMING\]\s+([A-Za-z0-9_]+):\s+([-+0-9.eE]+)")
TRANSFORM_LABEL = re.compile(r"alignTransform_(\d+_\d+-\d+)\.tfm")
RUN_START = re.compile(r"Running registration call\s+(\d+)")
NLopt = re.compile(r"NLopt elapsed time \(sec\)\s*:\s*([-+0-9.eE]+)")
CALL_ELAPSED = re.compile(
    r"(?:CUDA )?Registration call elapsed runtime \(sec\)\s*:\s*([-+0-9.eE]+)",
    re.IGNORECASE,
)


def _protocol_fields(text: str) -> dict[str, str]:
    result = {}
    for token in text.split():
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value.strip('"')
    return result


def parse_log_components(paths: list[Path]) -> tuple[dict[str, dict[str, float]], list[str]]:
    """Parse backend-native components, keyed only by explicit output labels."""
    records: dict[str, list[dict[str, float]]] = defaultdict(list)
    warnings: list[str] = []
    for path in paths:
        current: dict[str, Any] | None = None
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if RUN_START.search(line):
                if current and current.get("label") and current.get("timings"):
                    records[current["label"]].append(current["timings"])
                current = {"label": None, "timings": {}}
            protocol = PROTOCOL_LINE.search(line)
            if protocol:
                label, fields = protocol.groups()
                parsed = {}
                for key, value in _protocol_fields(fields).items():
                    if key.endswith("_ms"):
                        try:
                            parsed[f"persistent_{key}"] = float(value)
                        except ValueError:
                            warnings.append(f"Ignored invalid {key} in {path.name}")
                records[label].append(parsed)
                if current is not None:
                    current["label"] = label
            reference_load = LOAD_REFERENCE_LINE.search(line)
            if reference_load and current is not None:
                fields = _protocol_fields(reference_load.group(1))
                try:
                    current["timings"]["persistent_reference_initialization_ms"] = float(fields["total_ms"])
                except (KeyError, ValueError):
                    warnings.append(f"Ignored invalid LOAD_REFERENCE total_ms in {path.name}")
            if current is None:
                continue
            label_match = TRANSFORM_LABEL.search(line)
            if label_match:
                current["label"] = label_match.group(1)
            timing = TIMING_LINE.search(line)
            if timing:
                key, value = timing.groups()
                try:
                    current["timings"][f"standalone_{key}"] = float(value)
                except ValueError:
                    warnings.append(f"Ignored invalid [TIMING] {key} in {path.name}")
            nlopt = NLopt.search(line)
            if nlopt:
                current["timings"]["sms_mi_reg_optimizer_ms"] = float(nlopt.group(1)) * 1000
            call_elapsed = CALL_ELAPSED.search(line)
            if call_elapsed:
                current["timings"]["log_full_registration_elapsed_ms"] = float(call_elapsed.group(1)) * 1000
        if current and current.get("label") and current.get("timings"):
            records[current["label"]].append(current["timings"])

    unique: dict[str, dict[str, float]] = {}
    for label, candidates in records.items():
        merged: dict[str, float] = {}
        ambiguous = False
        for candidate in candidates:
            for key, value in candidate.items():
                if key in merged and abs(merged[key] - value) > 1e-9:
                    ambiguous = True
                merged[key] = value
        if ambiguous:
            warnings.append(f"Conflicting component timing records for output label {label}; components omitted.")
        else:
            unique[label] = merged
    return unique, warnings


def parse_queue_items(path: Path) -> list[dict[str, Any]]:
    source_rows, _ = read_csv(path, REQUIRED)
    rows = []
    for source in source_rows:
        start = optional_int(source["registration_start_perf_ns"], "registration_start_perf_ns")
        end = optional_int(source["registration_end_perf_ns"], "registration_end_perf_ns")
        if start is None and end is None:
            continue
        if start is None or end is None or end < start:
            raise AnalysisError(f"Invalid registration timing interval for {source.get('pointer_filename', source.get('output_label'))}")
        rows.append({
            "registration_index": int(source["output_label"].split("_", 1)[0]) if re.fullmatch(r"\d+_.*", source["output_label"]) else None,
            "output_label": source["output_label"],
            "volume": required_int(source["volume"], "volume"),
            "group": required_int(source["group"], "group"),
            "outcome": source["outcome"],
            "reg_engine": source["reg_engine"],
            "cuda_execution_mode": source["cuda_execution_mode"],
            "full_registration_elapsed_ms": (end - start) / 1e6,
            "_sort_ns": start,
        })
    rows.sort(key=lambda row: row["_sort_ns"])
    add_temporal_thirds(rows)
    return rows


def log_only_rows(components: dict[str, dict[str, float]]) -> list[dict[str, Any]]:
    rows = []
    label_pattern = re.compile(r"^(\d+)_(\d+)-(\d+)$")
    for label, values in components.items():
        match = label_pattern.fullmatch(label)
        elapsed = values.get("log_full_registration_elapsed_ms")
        if match is None or elapsed is None:
            continue
        index, volume, group = (int(value) for value in match.groups())
        if any(key.startswith("persistent_") for key in values):
            engine, mode = "cuda", "persistent"
        elif any(key.startswith("standalone_") for key in values):
            engine, mode = "cuda", "standalone"
        else:
            engine, mode = "sms-mi-reg", "not_applicable"
        row = {
            "registration_index": index, "output_label": label, "volume": volume,
            "group": group, "outcome": "unavailable", "reg_engine": engine,
            "cuda_execution_mode": mode, "full_registration_elapsed_ms": elapsed,
            "timing_source": "queue_log", "_sort_ns": index,
        }
        row.update(values)
        rows.append(row)
    rows.sort(key=lambda row: row["_sort_ns"])
    add_temporal_thirds(rows)
    return rows


def apply_status_outcomes(rows: list[dict[str, Any]], path: Path) -> None:
    document = read_json(path)
    if not isinstance(document, dict) or document.get("schema_version") != 1:
        raise AnalysisError(f"Unsupported registration status schema in {path.name}")
    volumes = document.get("volumes")
    if not isinstance(volumes, dict):
        raise AnalysisError(f"Required field 'volumes' is missing from {path.name}")
    outcomes = {}
    for volume, record in volumes.items():
        if not isinstance(record, dict) or not isinstance(record.get("groups"), dict):
            continue
        for group, group_record in record["groups"].items():
            if isinstance(group_record, dict) and isinstance(group_record.get("status"), str):
                outcomes[(int(volume), int(group))] = group_record["status"]
    for row in rows:
        if row["outcome"] == "unavailable":
            row["outcome"] = outcomes.get((row["volume"], row["group"]), "unavailable")


def split_reference_initialization(
    rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Exclude only records explicitly carrying a persistent reference load."""
    excluded = [row for row in rows if row.get("persistent_reference_initialization_ms") is not None]
    included = [row for row in rows if row.get("persistent_reference_initialization_ms") is None]
    if not included:
        return rows, []
    return included, excluded


def _plots(rows: list[dict[str, Any]], components: list[str], output: Path) -> dict[str, Any]:
    plt = pyplot()
    plotted, excluded = split_reference_initialization(rows)
    full_values = [row["full_registration_elapsed_ms"] for row in plotted]
    fig, axes, full_break = adaptive_y_axes(plt, full_values, figsize=(10, 5.5))
    for ax in axes:
        ax.scatter(
            range(1, len(plotted) + 1), full_values,
            s=9, alpha=0.50, color="#3976af",
        )
    axes[-1].set_xlabel("Registration opportunity (temporal order)")
    axes[0].set_title(
        "Registration-call elapsed time"
        + (" (explicit reference initialization excluded)" if excluded else "")
        + (" (broken y-axis)" if full_break else "")
    )
    fig.supylabel("Full registration call time (ms)")
    if full_break:
        fig.subplots_adjust(left=0.11, right=0.97, bottom=0.12, top=0.88, hspace=0.06)
    else:
        fig.tight_layout(rect=(0.02, 0, 1, 1))
    fig.savefig(output / "registration_call_timing.png", dpi=150)
    plt.close(fig)
    component_break = None
    component_point_count = 0
    if components:
        component_values = [
            row[field] for row in plotted for field in components
            if row.get(field) is not None
        ]
        fig, axes, component_break = adaptive_y_axes(
            plt, component_values, figsize=(12, 6)
        )
        for field in components:
            points = [(r["volume"], r.get(field)) for r in plotted if r.get(field) is not None]
            if points:
                component_point_count += len(points)
                for index, ax in enumerate(axes):
                    ax.scatter(
                        [x for x, _ in points], [y for _, y in points],
                        s=8, alpha=0.40, label=field.removesuffix("_ms") if index == 0 else None,
                    )
        axes[-1].set_xlabel("Volume number")
        axes[0].set_title(
            "Backend-native registration timing components"
            + (" (explicit reference initialization excluded)" if excluded else "")
            + (" (broken y-axis)" if component_break else "")
        )
        fig.supylabel("Time (ms)")
        handles, labels = axes[0].get_legend_handles_labels()
        axes[-1].legend(
            handles,
            labels,
            fontsize="x-small",
            ncol=len(labels),
            loc="upper center",
            bbox_to_anchor=(0.5, -0.22),
            borderaxespad=0.0,
        )
        if component_break:
            fig.subplots_adjust(left=0.10, right=0.97, bottom=0.25, top=0.88, hspace=0.06)
        else:
            fig.subplots_adjust(left=0.08, right=0.97, bottom=0.25, top=0.90)
        fig.savefig(output / "registration_timing_by_volume.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    return {
        "excluded_initialization": excluded,
        "full_break": full_break,
        "component_break": component_break,
        "component_point_count": component_point_count,
    }


def run(input_directory: Path, output_arg: Path | None = None) -> Path:
    input_directory = validate_input_directory(input_directory)
    item_candidates = sorted(input_directory.glob("queue_profile_items_*.csv"))
    if len(item_candidates) > 1:
        raise AnalysisError(
            "Ambiguous queue item profile: "
            + ", ".join(path.name for path in item_candidates)
        )
    item_path = item_candidates[0] if item_candidates else None
    log_paths = sorted(input_directory.glob("log_local_queue_processor_*.log"))
    status_path = input_directory / "registration_status.json"
    if not status_path.is_file():
        status_path = None
    _, output = resolve_output_directory(input_directory, output_arg, "registration_calls")
    components, warnings = parse_log_components(log_paths) if log_paths else ({}, ["No queue log was found; backend-native timing components are unavailable."])
    if item_path is not None:
        rows = parse_queue_items(item_path)
        for row in rows:
            row["timing_source"] = "queue_profile"
            row.update(components.get(row["output_label"], {}))
    else:
        rows = log_only_rows(components)
        warnings.append(
            "queue_profile_items_*.csv is unavailable; full-call timing uses the queue log's run_MIregistration timer."
        )
    if not rows:
        required = "queue_profile_items_*.csv or an explicitly labeled registration timing block in log_local_queue_processor_*.log"
        raise AnalysisError(f"No registration-call timing records were found; required input is {required}.")
    if status_path is not None:
        apply_status_outcomes(rows, status_path)
    component_fields = sorted({key for row in rows for key in row if key.endswith("_ms")} - {"full_registration_elapsed_ms"})
    clean_rows = [{k: v for k, v in row.items() if k != "_sort_ns"} for row in rows]
    all_fields = ["registration_index", "output_label", "volume", "group", "outcome", "reg_engine", "cuda_execution_mode", "timing_source", "third", "full_registration_elapsed_ms", *component_fields]
    write_csv(output / "registration_call_records.csv", clean_rows, all_fields)
    metrics = {"full_registration_elapsed_ms": [r["full_registration_elapsed_ms"] for r in rows]}
    metrics.update({field: [row.get(field) for row in rows] for field in component_fields})
    write_csv(output / "registration_timing_statistics.csv", statistics_rows(metrics), ["metric", *STAT_FIELDS])
    write_csv(output / "registration_timing_by_temporal_third.csv", temporal_statistics_rows(rows, ["full_registration_elapsed_ms", *component_fields]), ["metric", "third", *STAT_FIELDS])
    outliers = sorted(clean_rows, key=lambda row: row["full_registration_elapsed_ms"], reverse=True)[:10]
    write_csv(output / "registration_call_outliers.csv", outliers, all_fields)
    plot_info = _plots(rows, component_fields, output)
    full = describe(r["full_registration_elapsed_ms"] for r in rows)
    thirds = third_statistics(rows, "full_registration_elapsed_ms")
    modes = Counter(f"{r['reg_engine']}/{r['cuda_execution_mode']}" for r in rows)
    unavailable = []
    expected_groups = {
        "persistent CUDA request/optimizer/target-read/target-add": any(field.startswith("persistent_") for field in component_fields),
        "standalone CUDA internal phases": any(field.startswith("standalone_") for field in component_fields),
        "sms-mi-reg optimizer": "sms_mi_reg_optimizer_ms" in component_fields,
    }
    for label, present in expected_groups.items():
        if not present:
            unavailable.append(f"{label} timing unavailable for this run")
    semantics = [
        "With a queue profile, full_registration_elapsed_ms is registration_end_perf_ns - registration_start_perf_ns, immediately around run_MIregistration in the queue process.",
        "When the queue profile is absent, full_registration_elapsed_ms falls back to the explicitly labeled queue-log run_MIregistration wall-clock timer; timing_source identifies this case.",
        "persistent_request_ms is the backend-reported REGISTER request scope; persistent optimizer/target fields are nested components from the same protocol response.",
        "standalone_total_ms and standalone phase fields are emitted inside the CUDA executable; the queue full-call scope additionally includes subprocess/Python overhead.",
        "No residual overhead is derived: backend component scope/non-overlap cannot be proven from source code present in this repository.",
    ]
    results = [
        f"Registration calls: {len(rows)}; configurations: {dict(sorted(modes.items()))}",
        format_stats("Full registration call", full),
        f"Full-call late/early median ratio: {format_number(late_to_early_ratio(thirds))}",
        f"Backend timing components found: {', '.join(component_fields) if component_fields else 'none'}",
        (
            "Registration timing plots excluded explicitly identified persistent reference initialization: "
            + ", ".join(row["output_label"] for row in plot_info["excluded_initialization"])
            if plot_info["excluded_initialization"] else
            "No distinct reference initialization was explicitly identifiable; registration timing plots retain all calls and use adaptive broken axes only when warranted."
        ),
    ]
    artifacts = [*([item_path] if item_path is not None else []), *log_paths, *([status_path] if status_path is not None else [])]
    write_summary(output / "summary.txt", [
        ("Provenance", provenance_lines(Path(__file__).name, input_directory, artifacts, len(rows), optional_unavailable=unavailable, assumptions=semantics)),
        ("High-level results", results),
        ("Timer semantics", semantics),
        ("Unavailable timing", unavailable or ["none"]),
        ("Temporal thirds", [f"{name}: {format_stats('full call', thirds[name])}" for name in ("early", "middle", "late")]),
    ])
    finish_terminal(input_directory, artifacts, output, f"Analyzed {len(rows)} registration calls.", warnings + unavailable)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    args = parser.parse_args()
    run(args.input_directory, args.output_dir)


if __name__ == "__main__":
    cli_main(main)

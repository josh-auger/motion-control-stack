#!/usr/bin/env python3
"""Shared, read-only helpers for acquisition post-run analysis."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import sys
import warnings
from collections.abc import Iterable, Sequence
from datetime import datetime
from pathlib import Path
from typing import Any


class AnalysisError(RuntimeError):
    """An actionable input/artifact error suitable for CLI display."""


STAT_FIELDS = ("count", "mean", "median", "std", "p95", "p99", "min", "max")
THIRD_NAMES = ("early", "middle", "late")


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input_directory", type=Path, help="completed savedData_* directory")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="analysis root; relative paths are resolved below the input directory",
    )


def validate_input_directory(path: Path) -> Path:
    path = path.expanduser().resolve()
    if not path.exists():
        raise AnalysisError(f"Input directory does not exist: {path}")
    if not path.is_dir():
        raise AnalysisError(f"Input path is not a directory: {path}")
    return path


def resolve_output_directory(
    input_directory: Path,
    requested: Path | None,
    analysis_name: str,
    *,
    now: datetime | None = None,
) -> tuple[Path, Path]:
    if requested is None:
        stamp = (now or datetime.now()).strftime("%Y%m%d_%H%M%S")
        root = input_directory / f"analysis_{stamp}"
    elif requested.is_absolute():
        root = requested.expanduser()
    else:
        root = input_directory / requested
    root = root.resolve()
    output = root / analysis_name
    output.mkdir(parents=True, exist_ok=True)
    return root, output


def matching_files(directory: Path, pattern: str) -> list[Path]:
    return sorted(path for path in directory.glob(pattern) if path.is_file())


def discover_unique(
    directory: Path,
    pattern: str,
    description: str,
    *,
    missing_hint: str | None = None,
) -> Path:
    matches = matching_files(directory, pattern)
    if not matches:
        hint = f" {missing_hint}" if missing_hint else ""
        raise AnalysisError(
            f"Required {description} was not found ({pattern}) in {directory}.{hint}"
        )
    if len(matches) != 1:
        candidates = ", ".join(path.name for path in matches)
        raise AnalysisError(
            f"Ambiguous {description}: {len(matches)} files match {pattern!r}: "
            f"{candidates}. Their contents do not provide an acquisition identifier "
            "that safely selects one; isolate the acquisition artifacts and retry."
        )
    return matches[0]


def discover_queue_profile_pair(directory: Path) -> tuple[Path, Path]:
    items = matching_files(directory, "queue_profile_items_*.csv")
    scans = matching_files(directory, "queue_profile_scans_*.csv")
    if not items:
        raise AnalysisError(
            f"Required queue item profile was not found "
            f"(queue_profile_items_*.csv) in {directory}. QUEUE_PROFILE_FLAG may have been disabled."
        )
    if not scans:
        raise AnalysisError(
            f"Required queue scan profile was not found "
            f"(queue_profile_scans_*.csv) in {directory}. QUEUE_PROFILE_FLAG may have been disabled."
        )
    item_by_suffix = {p.name.removeprefix("queue_profile_items_"): p for p in items}
    scan_by_suffix = {p.name.removeprefix("queue_profile_scans_"): p for p in scans}
    shared = sorted(set(item_by_suffix) & set(scan_by_suffix))
    if len(shared) == 1 and len(items) == 1 and len(scans) == 1:
        suffix = shared[0]
        return item_by_suffix[suffix], scan_by_suffix[suffix]
    if len(shared) == 1:
        # A shared writer suffix is a deterministic pair only when every other
        # candidate is unpaired residue.
        suffix = shared[0]
        return item_by_suffix[suffix], scan_by_suffix[suffix]
    candidates = [p.name for p in items + scans]
    raise AnalysisError(
        "Ambiguous queue profiles: could not identify exactly one item/scan pair "
        f"from shared writer suffixes. Candidates: {', '.join(candidates)}"
    )


def optional_files(directory: Path, pattern: str) -> list[Path]:
    return matching_files(directory, pattern)


def read_csv(path: Path, required_columns: Iterable[str]) -> tuple[list[dict[str, str]], list[str]]:
    try:
        with path.open("r", newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream, strict=True)
            fields = list(reader.fieldnames or [])
            missing = sorted(set(required_columns) - set(fields))
            if missing:
                raise AnalysisError(
                    f"Required CSV column(s) missing from {path.name}: {', '.join(missing)}"
                )
            rows = []
            for line_number, row in enumerate(reader, 2):
                if None in row:
                    raise AnalysisError(
                        f"Malformed CSV row in {path.name} at line {line_number}: too many fields"
                    )
                rows.append(dict(row))
    except AnalysisError:
        raise
    except (OSError, UnicodeError, csv.Error) as error:
        raise AnalysisError(f"Could not read CSV {path}: {error}") from error
    return rows, fields


def read_json(path: Path) -> Any:
    try:
        with path.open("r", encoding="utf-8") as stream:
            return json.load(stream)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise AnalysisError(f"Could not read JSON {path}: {error}") from error


def required_int(value: Any, field: str, *, context: str = "") -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        suffix = f" ({context})" if context else ""
        raise AnalysisError(f"Invalid integer in {field}{suffix}: {value!r}") from error


def optional_int(value: Any, field: str = "value") -> int | None:
    return None if value in (None, "") else required_int(value, field)


def required_float(value: Any, field: str, *, context: str = "") -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        suffix = f" ({context})" if context else ""
        raise AnalysisError(f"Invalid number in {field}{suffix}: {value!r}") from error
    if not math.isfinite(result):
        raise AnalysisError(f"Non-finite number in {field}: {value!r}")
    return result


def optional_float(value: Any, field: str = "value") -> float | None:
    return None if value in (None, "") else required_float(value, field)


def percentile(sorted_values: Sequence[float], percent: float) -> float:
    if not sorted_values:
        raise ValueError("percentile requires at least one value")
    position = (len(sorted_values) - 1) * percent / 100.0
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(sorted_values[lower])
    fraction = position - lower
    return float(sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction)


def describe(values: Iterable[float | int | None]) -> dict[str, float | int | None]:
    cleaned = sorted(
        float(value)
        for value in values
        if value is not None and math.isfinite(float(value))
    )
    if not cleaned:
        return {field: (0 if field == "count" else None) for field in STAT_FIELDS}
    return {
        "count": len(cleaned),
        "mean": statistics.fmean(cleaned),
        "median": statistics.median(cleaned),
        "std": statistics.pstdev(cleaned),
        "p95": percentile(cleaned, 95),
        "p99": percentile(cleaned, 99),
        "min": cleaned[0],
        "max": cleaned[-1],
    }


def temporal_third(index: int, count: int) -> str:
    if count <= 0 or index < 0 or index >= count:
        raise ValueError("temporal-third index must be inside a nonempty sequence")
    # Equivalent to three ordered chunks whose sizes differ by at most one,
    # with any remainder assigned to earlier chunks (numpy.array_split).
    base, remainder = divmod(count, 3)
    boundary_1 = base + (1 if remainder > 0 else 0)
    boundary_2 = boundary_1 + base + (1 if remainder > 1 else 0)
    return "early" if index < boundary_1 else "middle" if index < boundary_2 else "late"


def add_temporal_thirds(rows: list[dict[str, Any]], field: str = "third") -> None:
    for index, row in enumerate(rows):
        row[field] = temporal_third(index, len(rows))


def third_statistics(rows: Sequence[dict[str, Any]], value_field: str) -> dict[str, dict[str, Any]]:
    return {
        name: describe(row.get(value_field) for row in rows if row.get("third") == name)
        for name in THIRD_NAMES
    }


def late_to_early_ratio(stats: dict[str, dict[str, Any]], key: str = "median") -> float | None:
    early = stats.get("early", {}).get(key)
    late = stats.get("late", {}).get(key)
    if early in (None, 0) or late is None:
        return None
    return float(late) / float(early)


def pearson(xs: Iterable[float], ys: Iterable[float]) -> float | None:
    pairs = [(float(x), float(y)) for x, y in zip(xs, ys)]
    if len(pairs) < 2:
        return None
    x_values, y_values = zip(*pairs)
    x_mean, y_mean = statistics.fmean(x_values), statistics.fmean(y_values)
    x_deviation = [value - x_mean for value in x_values]
    y_deviation = [value - y_mean for value in y_values]
    denominator = math.sqrt(
        sum(value * value for value in x_deviation)
        * sum(value * value for value in y_deviation)
    )
    if denominator == 0:
        return None
    return sum(x * y for x, y in zip(x_deviation, y_deviation)) / denominator


def grouped_medians(
    rows: Sequence[dict[str, Any]], group_field: str, value_field: str
) -> list[tuple[Any, float]]:
    """Return one median per sorted group, ignoring missing values."""
    grouped: dict[Any, list[float]] = {}
    for row in rows:
        group, value = row.get(group_field), row.get(value_field)
        if group is not None and value is not None:
            grouped.setdefault(group, []).append(float(value))
    return [
        (group, float(statistics.median(values)))
        for group, values in sorted(grouped.items())
    ]


def running_median(
    values: Sequence[float], window: int | None = None
) -> tuple[list[float], int]:
    """Return a centered running median and the odd window used.

    The default window is the nearest odd integer at or above sqrt(n), bounded
    to 5..101 records when possible and never larger than the dataset. Edge
    windows are truncated rather than padded. This follows acquisition order
    without assuming an acquisition duration or number of volumes.
    """
    cleaned = [float(value) for value in values]
    if not cleaned:
        return [], 0
    if window is None:
        window = max(5, min(101, math.ceil(math.sqrt(len(cleaned)))))
    if window <= 0:
        raise ValueError("running-median window must be positive")
    if window % 2 == 0:
        window += 1
    if window > len(cleaned):
        window = len(cleaned) if len(cleaned) % 2 else max(1, len(cleaned) - 1)
    half = window // 2
    medians = [
        float(statistics.median(cleaned[max(0, index - half):min(len(cleaned), index + half + 1)]))
        for index in range(len(cleaned))
    ]
    return medians, window


def ordinary_least_squares(
    xs: Sequence[float], ys: Sequence[float]
) -> dict[str, float | int] | None:
    """Fit y = slope*x + intercept and return R-squared."""
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    x_values, y_values = [float(x) for x in xs], [float(y) for y in ys]
    x_mean, y_mean = statistics.fmean(x_values), statistics.fmean(y_values)
    denominator = sum((x - x_mean) ** 2 for x in x_values)
    if denominator == 0:
        return None
    slope = sum(
        (x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values)
    ) / denominator
    intercept = y_mean - slope * x_mean
    predictions = [slope * x + intercept for x in x_values]
    residual = sum((y - predicted) ** 2 for y, predicted in zip(y_values, predictions))
    total = sum((y - y_mean) ** 2 for y in y_values)
    r_squared = 1.0 - residual / total if total else 1.0
    return {
        "count": len(x_values),
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
    }


def adaptive_break_limits(
    values: Iterable[float | int | None],
) -> dict[str, tuple[float, float] | float | int] | None:
    """Find a clearly separated, sparse upper outlier region.

    A break is used only for at least eight observations. The largest gap among
    splits that leave no more than 10% of observations in the upper region must
    begin above the outer Tukey fence (Q3 + 3*IQR), and must be at least both
    2*IQR and half the full lower-cluster span. This deliberately avoids breaking
    ordinary continuous or broadly multimodal distributions while still finding
    one extreme point beyond several merely high observations.
    """
    cleaned = sorted(
        float(value)
        for value in values
        if value is not None and math.isfinite(float(value))
    )
    if len(cleaned) < 8:
        return None
    q1, q3 = percentile(cleaned, 25), percentile(cleaned, 75)
    iqr = q3 - q1
    if iqr <= 0:
        return None
    fence = q3 + 3.0 * iqr
    maximum_upper = max(1, math.ceil(len(cleaned) * 0.10))
    candidates = []
    for upper_count in range(1, maximum_upper + 1):
        split = len(cleaned) - upper_count
        gap = cleaned[split] - cleaned[split - 1]
        candidates.append((gap, split))
    _candidate_gap, split = max(candidates)
    lower, upper = cleaned[:split], cleaned[split:]
    if upper[0] <= fence:
        return None
    gap = upper[0] - lower[-1]
    lower_span = lower[-1] - lower[0]
    if gap < max(2.0 * iqr, 0.5 * lower_span):
        return None
    lower_padding = max(0.04 * max(lower_span, iqr), 1e-9)
    upper_span = upper[-1] - upper[0]
    upper_padding = max(0.08 * max(upper_span, iqr), 1e-9)
    return {
        "lower_ylim": (lower[0] - lower_padding, lower[-1] + lower_padding),
        "upper_ylim": (upper[0] - upper_padding, upper[-1] + upper_padding),
        "lower_max": lower[-1],
        "upper_min": upper[0],
        "upper_count": len(upper),
    }


def adaptive_y_axes(plt, values, *, figsize=(10, 5)):
    """Create one axis or aligned broken-y axes using adaptive_break_limits."""
    break_info = adaptive_break_limits(values)
    if break_info is None:
        fig, axis = plt.subplots(figsize=figsize)
        return fig, [axis], None
    fig, (upper, lower) = plt.subplots(
        2, 1, sharex=True, figsize=figsize,
        gridspec_kw={"height_ratios": [1, 3], "hspace": 0.06},
    )
    upper.set_ylim(*break_info["upper_ylim"])
    lower.set_ylim(*break_info["lower_ylim"])
    upper.spines["bottom"].set_visible(False)
    lower.spines["top"].set_visible(False)
    upper.tick_params(labelbottom=False, bottom=False)
    lower.xaxis.tick_bottom()
    marker = 0.012
    kwargs = dict(color="black", clip_on=False, linewidth=1.0)
    upper.plot((-marker, +marker), (-marker, +marker), transform=upper.transAxes, **kwargs)
    upper.plot((1 - marker, 1 + marker), (-marker, +marker), transform=upper.transAxes, **kwargs)
    lower.plot((-marker, +marker), (1 - marker, 1 + marker), transform=lower.transAxes, **kwargs)
    lower.plot((1 - marker, 1 + marker), (1 - marker, 1 + marker), transform=lower.transAxes, **kwargs)
    return fig, [upper, lower], break_info


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field) for field in fields})


def statistics_rows(metrics: dict[str, Iterable[float | int | None]]) -> list[dict[str, Any]]:
    return [{"metric": name, **describe(values)} for name, values in metrics.items()]


def temporal_statistics_rows(
    rows: Sequence[dict[str, Any]], fields: Sequence[str]
) -> list[dict[str, Any]]:
    output = []
    for field in fields:
        for third in THIRD_NAMES:
            output.append({
                "metric": field,
                "third": third,
                **describe(row.get(field) for row in rows if row.get("third") == third),
            })
    return output


def format_number(value: Any, digits: int = 3) -> str:
    if value is None:
        return "unavailable"
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.{digits}f}"


def format_stats(label: str, stats: dict[str, Any], unit: str = "ms") -> str:
    return (
        f"{label}: n={stats['count']}, mean={format_number(stats['mean'])} {unit}, "
        f"median={format_number(stats['median'])} {unit}, std={format_number(stats['std'])} {unit}, "
        f"p95={format_number(stats['p95'])} {unit}, p99={format_number(stats['p99'])} {unit}, "
        f"min={format_number(stats['min'])} {unit}, max={format_number(stats['max'])} {unit}"
    )


def provenance_lines(
    script_name: str,
    input_directory: Path,
    artifacts: Sequence[Path],
    record_count: int,
    *,
    optional_unavailable: Sequence[str] = (),
    assumptions: Sequence[str] = (),
) -> list[str]:
    lines = [
        f"Analysis script: {script_name}",
        f"Analysis timestamp: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"Input directory: {input_directory.resolve()}",
        "Source artifacts: " + (", ".join(path.name for path in artifacts) or "none"),
        f"Records analyzed: {record_count}",
    ]
    if optional_unavailable:
        lines.append("Optional artifacts/data unavailable: " + "; ".join(optional_unavailable))
    if assumptions:
        lines.append("Assumptions/limitations:")
        lines.extend(f"- {item}" for item in assumptions)
    return lines


def write_summary(path: Path, sections: Sequence[tuple[str, Sequence[str]]]) -> None:
    output: list[str] = []
    for title, lines in sections:
        if output:
            output.append("")
        output.append(title)
        output.append("=" * len(title))
        output.extend(str(line) for line in lines)
    path.write_text("\n".join(output) + "\n", encoding="utf-8")


def pyplot():
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    warnings.filterwarnings("ignore", message="Unable to import Axes3D")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise AnalysisError(
            "Plotting requires matplotlib, which is listed in the project requirements."
        ) from error
    return plt


def finish_terminal(
    input_directory: Path,
    artifacts: Sequence[Path],
    output_directory: Path,
    headline: str,
    warnings: Sequence[str] = (),
) -> None:
    print(f"Input directory: {input_directory}")
    print("Artifacts used: " + ", ".join(path.name for path in artifacts))
    print(f"Output directory: {output_directory}")
    for warning in warnings:
        print(f"WARNING: {warning}")
    print(headline)


def cli_main(function) -> None:
    try:
        function()
    except AnalysisError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2) from error

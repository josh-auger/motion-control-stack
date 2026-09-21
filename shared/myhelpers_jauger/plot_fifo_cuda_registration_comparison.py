"""Plot the FIFO-on CUDA A/B comparison from saved data and analysis tables.

Run from the repository root:
    python3 shared/myhelpers_jauger/plot_fifo_cuda_registration_comparison.py \
        --standalone-dir data/savedData_20260921_fetal1_SMS1_TR2.5_fifo_standalone \
        --persistent-dir data/savedData_20260921_fetal1_SMS1_TR2.5_fifo_persistent \
        --analysis-dir data/analysis_20260921_fifo_cuda_comparison \
        --logs-dir data/logs_motioncontrolstack

Reads the comparison CSV/JSON, original motion CSVs, and entrypoint logs.
Writes figures only to --analysis-dir; does not overwrite analysis tables.
"""

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


HERE = None
RUNS = {}
SUMMARY = None
BLUE = "#1768AC"
ORANGE = "#D05A2D"
NEUTRAL = "#52606D"
AXIS_COLORS = ("#315B8A", "#A04766", "#16866A")
DPI = 180

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "figure.titlesize": 15,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "savefig.facecolor": "white",
})


def read_csv(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--standalone-dir", type=Path, required=True,
                        help="Saved acquisition directory for standalone CUDA")
    parser.add_argument("--persistent-dir", type=Path, required=True,
                        help="Saved acquisition directory for persistent CUDA")
    parser.add_argument("--analysis-dir", type=Path, required=True,
                        help="Directory containing comparison tables; figures are saved here")
    parser.add_argument("--logs-dir", type=Path, required=True,
                        help="Directory containing the matching entrypoint logs")
    args = parser.parse_args(argv)
    for name in ("standalone_dir", "persistent_dir", "analysis_dir", "logs_dir"):
        path = getattr(args, name).expanduser().resolve()
        if not path.is_dir():
            parser.error(f"--{name.replace('_', '-')} is not a directory: {path}")
        setattr(args, name, path)
    if args.standalone_dir == args.persistent_dir:
        parser.error("--standalone-dir and --persistent-dir must differ")
    if any(args.analysis_dir.is_relative_to(path)
           for path in (args.standalone_dir, args.persistent_dir)):
        parser.error("--analysis-dir must be outside both acquisition directories")
    return args


def configured_motion_threshold(logs_dir):
    matched = {"standalone": [], "persistent": []}
    for path in sorted(logs_dir.glob("log_entrypoint*")):
        content = path.read_text()
        if not re.search(r"^\s*REG_ENGINE=cuda\s*$", content, re.MULTILINE):
            continue
        if not re.search(r"^\s*FIFO_FLAG=on\s*$", content, re.MULTILINE):
            continue
        for mode in matched:
            if re.search(rf"^\s*CUDA_EXECUTION_MODE={mode}\s*$", content, re.MULTILINE):
                matched[mode].append((path, content))
    values = []
    for mode, entries in matched.items():
        if len(entries) != 1:
            raise ValueError(f"Expected one FIFO-on CUDA {mode} entrypoint log in {logs_dir}; found {len(entries)}")
        path, content = entries[0]
        for key, expected in (("REG_ENGINE", "cuda"), ("CUDA_EXECUTION_MODE", mode), ("FIFO_FLAG", "on")):
            observed = re.search(rf"^\s*{key}=(\S+)\s*$", content, re.MULTILINE)
            assert observed and observed.group(1) == expected, (path, key)
        threshold = re.search(r"^\s*MOTION_THRESH=(\S+)\s*$", content, re.MULTILINE)
        assert threshold, path
        values.append(float(threshold.group(1)))
    assert values[0] == values[1] == 0.3
    return values[0]


def same(actual, expected, name):
    assert np.isclose(actual, expected, rtol=1e-10, atol=1e-10), (name, actual, expected)


def save(fig, filename):
    path = HERE / filename
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(path)


def mode_legend(fig, *, y=0.95):
    handles = [
        Line2D([0], [0], color=BLUE, lw=2, label="Standalone CUDA"),
        Line2D([0], [0], color=ORANGE, lw=2, label="Persistent CUDA"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, y),
               frameon=False, ncol=2)


def load_and_validate(threshold):
    timing = read_csv(HERE / "registration_timing_comparison.csv")
    transform = read_csv(HERE / "transform_difference_summary.csv")
    motion_difference = read_csv(HERE / "motion_trace_difference_summary.csv")
    motion = {}
    for mode in ("standalone", "persistent"):
        paths = list(RUNS[mode].glob("motionMonitor_data*.csv"))
        if len(paths) != 1:
            raise ValueError(f"Expected one motion CSV in {RUNS[mode]}; found {len(paths)}")
        rows = read_csv(paths[0])
        keyed = {(int(row["Volume_index"]), int(row["Slice_group_index"])): row for row in rows}
        assert len(keyed) == len(rows) == SUMMARY["runs"][mode]["motion_rows"] == 8787
        motion[mode] = keyed

    index = np.asarray([int(row["registration_index"]) for row in timing])
    assert len(index) == 8786 and np.array_equal(index, np.arange(1, 8787))
    times = {
        "standalone": np.asarray([float(row["standalone_ms"]) for row in timing]),
        "persistent": np.asarray([float(row["persistent_ms"]) for row in timing]),
    }
    for mode, values in times.items():
        expected = SUMMARY["runs"][mode]["queue_timing_ms"]
        assert len(values) == expected["n"]
        for field, actual in (("mean", np.mean(values)), ("median", np.median(values)),
                              ("min", np.min(values)), ("max", np.max(values)),
                              ("p95", np.percentile(values, 95))):
            same(actual, expected[field], f"{mode} timing {field}")

    transform_index = np.asarray([int(row["index"]) for row in transform])
    assert len(transform) == SUMMARY["comparison"]["transforms"]["matched_cuda_registrations"] == 8786
    assert np.array_equal(transform_index, index)
    rotations = np.abs(np.asarray([[float(row[f"d_{axis}"]) for axis in ("rx", "ry", "rz")] for row in transform]))
    translations = np.abs(np.asarray([[float(row[f"d_{axis}"]) for axis in ("tx", "ty", "tz")] for row in transform]))
    same(np.max(rotations), SUMMARY["comparison"]["transforms"]["rotation_rad"]["max_abs"], "max rotation difference")
    same(np.mean(rotations), SUMMARY["comparison"]["transforms"]["rotation_rad"]["mean_abs"], "mean rotation difference")
    same(np.max(translations), SUMMARY["comparison"]["transforms"]["translation_mm"]["max_abs"], "max translation difference")
    same(np.mean(translations), SUMMARY["comparison"]["transforms"]["translation_mm"]["mean_abs"], "mean translation difference")

    difference_keys = [(int(row["volume"]), int(row["group"])) for row in motion_difference]
    assert len(difference_keys) == len(set(difference_keys)) == 8787
    assert set(difference_keys) == set(motion["standalone"]) == set(motion["persistent"])
    ordered_keys = sorted(difference_keys, key=lambda key: int(motion["standalone"][key]["reg_index"]))
    motion_index = np.asarray([int(motion["standalone"][key]["reg_index"]) for key in ordered_keys])
    assert np.array_equal(motion_index, np.arange(1, 8788))
    for key in ordered_keys:
        assert motion["standalone"][key]["reg_index"] == motion["persistent"][key]["reg_index"]
    difference_by_key = dict(zip(difference_keys, motion_difference))

    motion_columns = (
        "X_rotation(rad)", "Y_rotation(rad)", "Z_rotation(rad)",
        "X_translation(mm)", "Y_translation(mm)", "Z_translation(mm)",
        "Displacement(mm)",
    )
    trace = {}
    for mode in ("standalone", "persistent"):
        trace[mode] = {column: np.asarray([float(motion[mode][key][column]) for key in ordered_keys])
                       for column in motion_columns}
        fd = trace[mode]["Displacement(mm)"]
        expected = SUMMARY["runs"][mode]
        same(np.mean(fd), expected["fd_mm"]["mean"], f"{mode} mean FD")
        same(np.max(fd), expected["fd_mm"]["max"], f"{mode} maximum FD")
        assert np.count_nonzero(fd > threshold) == expected["fd_above_0_3_count"]

    for column in motion_columns:
        previous_delta = np.asarray([float(difference_by_key[key][f"delta_{column}"]) for key in ordered_keys])
        current_delta = trace["persistent"][column] - trace["standalone"][column]
        assert np.allclose(previous_delta, current_delta, rtol=0, atol=1e-12), column
    return index, times, transform_index, rotations, translations, motion_index, trace


def registration_timing(index, times):
    a, b = times["standalone"], times["persistent"]
    mean_a, mean_b = float(np.mean(a)), float(np.mean(b))
    load_ms = float(SUMMARY["runs"]["persistent"]["persistent_load_details"][0]["total_ms"])
    fig, (top, detail) = plt.subplots(2, 1, figsize=(13.5, 7.1), sharex=True,
                                      gridspec_kw={"height_ratios": [2.5, 1]},
                                      constrained_layout=True)
    fig.suptitle("Registration wall time across the full FIFO-on acquisition")
    top.plot(index, a, color=BLUE, lw=0.55, alpha=0.72, label="Standalone CUDA")
    top.plot(index, b, color=ORANGE, lw=0.75, alpha=0.82, label="Persistent CUDA")
    top.axhline(mean_a, color=BLUE, ls="--", lw=1.5, label=f"Standalone mean: {mean_a:.2f} ms")
    top.axhline(mean_b, color=ORANGE, ls="--", lw=1.5, label=f"Persistent mean: {mean_b:.2f} ms")
    top.scatter([index[0]], [b[0]], color=ORANGE, edgecolor="white", s=50, zorder=5)
    top.text(0.70, 0.87, f"First persistent call: {b[0]:.2f} ms\nincludes {load_ms:.2f} ms reference load",
             transform=top.transAxes, ha="left", va="top", fontsize=9,
             bbox={"facecolor": "white", "edgecolor": "#D7DEE5", "boxstyle": "round,pad=0.4"})
    top.set(ylabel="Queue registration time (ms)", ylim=(0, max(a.max(), b.max()) * 1.06))
    top.legend(loc="upper right", bbox_to_anchor=(1, 0.76), frameon=True,
               facecolor="white", edgecolor="#D7DEE5", ncol=2, fontsize=8)
    top.grid(axis="y", color="#DCE3E9", lw=0.7)

    detail.plot(index, b, color=ORANGE, lw=0.6, alpha=0.8)
    detail.axhline(mean_b, color=ORANGE, ls="--", lw=1.2)
    detail.set(xlabel="Registration index", ylabel="Persistent detail (ms)", ylim=(0, 15.5), xlim=(1, len(index)))
    detail.text(0.99, 0.90, "0–15 ms detail; first call remains in full-scale panel",
                transform=detail.transAxes, ha="right", va="top", fontsize=9, color=NEUTRAL)
    detail.grid(axis="y", color="#DCE3E9", lw=0.7)
    save(fig, "registration_timing_vs_index.png")


def timing_distribution(times):
    a, b = times["standalone"], times["persistent"]
    fig, ax = plt.subplots(figsize=(8.5, 6.5), constrained_layout=True)
    fig.suptitle("Registration timing distributions retain every outlier")
    box = ax.boxplot([a, b], widths=0.48, patch_artist=True, showfliers=True,
                     showmeans=True, whis=1.5,
                     medianprops={"color": "#202B35", "linewidth": 1.8},
                     meanprops={"marker": "D", "markerfacecolor": "#202B35", "markeredgecolor": "white", "markersize": 6},
                     flierprops={"marker": "o", "markersize": 2.5, "markerfacecolor": NEUTRAL,
                                  "markeredgecolor": "none", "alpha": 0.6})
    for patch, color in zip(box["boxes"], (BLUE, ORANGE)):
        patch.set(facecolor=color, alpha=0.66, edgecolor=color, linewidth=1.4)
    for whisker, color in zip(box["whiskers"], (BLUE, BLUE, ORANGE, ORANGE)):
        whisker.set(color=color, linewidth=1.4)
    for cap, color in zip(box["caps"], (BLUE, BLUE, ORANGE, ORANGE)):
        cap.set(color=color, linewidth=1.4)
    ax.set_xticks((1, 2), (f"Standalone CUDA\nn={len(a):,}; median={np.median(a):.2f} ms",
                           f"Persistent CUDA\nn={len(b):,}; median={np.median(b):.2f} ms"))
    ax.set_yscale("log")
    ax.set_ylim(1, 500)
    ax.set_yticks((1, 2, 5, 10, 20, 50, 100, 200, 500))
    ax.set_yticklabels(("1", "2", "5", "10", "20", "50", "100", "200", "500"))
    ax.set_ylabel("Queue registration time (ms; logarithmic scale)")
    ax.grid(axis="y", which="both", color="#DCE3E9", lw=0.7)
    ax.legend(handles=[Line2D([0], [0], marker="D", linestyle="None", color="#202B35",
                              markersize=6, label="Mean")], loc="upper right", frameon=False)
    ax.text(0.5, 0.98, f"Means: {np.mean(a):.2f} ms standalone  •  {np.mean(b):.2f} ms persistent",
            transform=ax.transAxes, ha="center", va="top", fontsize=9, color=NEUTRAL)
    save(fig, "registration_timing_distribution.png")


def motion_parameters(index, trace):
    fig, axes = plt.subplots(3, 2, figsize=(15, 9.2), sharex=True)
    fig.suptitle("Aligned rigid-body motion parameters across the acquisition", y=0.995)
    mode_legend(fig, y=0.967)
    panels = [
        ("X_rotation(rad)", "rx rotation (rad)"), ("X_translation(mm)", "tx translation (mm)"),
        ("Y_rotation(rad)", "ry rotation (rad)"), ("Y_translation(mm)", "ty translation (mm)"),
        ("Z_rotation(rad)", "rz rotation (rad)"), ("Z_translation(mm)", "tz translation (mm)"),
    ]
    for ax, (column, label) in zip(axes.flat, panels):
        a, b = trace["standalone"][column], trace["persistent"][column]
        ax.plot(index, a, color=BLUE, lw=0.55, alpha=0.70)
        ax.plot(index, b, color=ORANGE, lw=0.55, alpha=0.70)
        ax.set_ylabel(label)
        low = min(0, float(a.min()), float(b.min()))
        high = max(0, float(a.max()), float(b.max()))
        padding = (high - low) * 0.05 or 0.01
        ax.set_ylim(low - padding, high + padding)
        ax.set_xlim(1, len(index))
        ax.axhline(0, color="#84919C", lw=0.65)
        ax.grid(axis="y", color="#E4E9ED", lw=0.65)
    for ax in axes[-1]:
        ax.set_xlabel("Aligned motion-sample index (1 includes reference identity)")
    fig.text(0.5, 0.015, "Each parameter uses its native scale; both modes share the same scale within each panel.",
             ha="center", fontsize=9, color=NEUTRAL)
    fig.tight_layout(rect=(0, 0.035, 1, 0.935), h_pad=1.1, w_pad=1.6)
    save(fig, "motion_parameter_trace_comparison.png")


def framewise_displacement(index, trace, threshold):
    a = trace["standalone"]["Displacement(mm)"]
    b = trace["persistent"]["Displacement(mm)"]
    difference = np.abs(b - a)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(13.5, 7), sharex=True,
                                      gridspec_kw={"height_ratios": [2.3, 1]},
                                      constrained_layout=True)
    fig.suptitle("Framewise displacement: overall trajectory and pointwise disagreement")
    top.plot(index, a, color=BLUE, lw=0.65, alpha=0.72, label="Standalone CUDA")
    top.plot(index, b, color=ORANGE, lw=0.65, alpha=0.72, label="Persistent CUDA")
    top.axhline(threshold, color=NEUTRAL, lw=1.2, ls="--", label=f"Motion threshold: {threshold:g} mm")
    top.set(ylabel="Framewise displacement (mm)", ylim=(0, max(a.max(), b.max()) * 1.07))
    top.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#D7DEE5", fontsize=9)
    top.text(0.01, 0.98,
             f"Standalone: mean {a.mean():.3f}, max {a.max():.3f} mm; above threshold {(a > threshold).sum():,}\n"
             f"Persistent: mean {b.mean():.3f}, max {b.max():.3f} mm; above threshold {(b > threshold).sum():,}",
             transform=top.transAxes, ha="left", va="top", fontsize=9,
             bbox={"facecolor": "white", "edgecolor": "#D7DEE5", "boxstyle": "round,pad=0.4"})
    top.grid(axis="y", color="#DCE3E9", lw=0.7)
    bottom.plot(index, difference, color="#62458A", lw=0.65, alpha=0.8)
    bottom.set(xlabel="Aligned motion-sample index (1 includes reference identity)",
               ylabel="Absolute A/B FD difference (mm)", ylim=(0, difference.max() * 1.08),
               xlim=(1, len(index)))
    maximum = int(np.argmax(difference))
    bottom.scatter([index[maximum]], [difference[maximum]], color="#62458A", s=30, zorder=4)
    bottom.text(0.99, 0.93, f"Maximum difference: {difference[maximum]:.3f} mm at sample {index[maximum]:,}",
                transform=bottom.transAxes, ha="right", va="top", fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "#D7DEE5", "boxstyle": "round,pad=0.3"})
    bottom.grid(axis="y", color="#DCE3E9", lw=0.7)
    save(fig, "framewise_displacement_comparison.png")


def transform_differences(index, rotations, translations):
    fig, axes = plt.subplots(2, 1, figsize=(13.5, 7.2), sharex=True, constrained_layout=True)
    fig.suptitle("Absolute differences between matched CUDA transform parameters")
    for ax, values, names, unit, maximum in (
        (axes[0], rotations, ("rx", "ry", "rz"), "rad", SUMMARY["comparison"]["transforms"]["largest_rotation"]),
        (axes[1], translations, ("tx", "ty", "tz"), "mm", SUMMARY["comparison"]["transforms"]["largest_translation"]),
    ):
        for column, name, color in zip(values.T, names, AXIS_COLORS):
            ax.plot(index, column, color=color, lw=0.58, alpha=0.72, label=f"|Δ{name}|")
        ax.set(ylabel=f"Absolute parameter difference ({unit})", ylim=(0, values.max() * 1.09), xlim=(1, len(index)))
        ax.grid(axis="y", color="#DCE3E9", lw=0.7)
        ax.legend(loc="upper right", bbox_to_anchor=(1, 0.78), frameon=True,
                  facecolor="white", edgecolor="#D7DEE5", ncol=3, fontsize=9)
        ax.text(0.99, 0.96,
                f"Maximum: {values.max():.4f} {unit}  •  index {maximum['index']:,}, volume {maximum['volume']}, group {maximum['group']}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "#D7DEE5", "boxstyle": "round,pad=0.3"})
    axes[1].set_xlabel("Matched registration index")
    save(fig, "transform_parameter_differences.png")


def acquisition_summary():
    a, b = SUMMARY["runs"]["standalone"], SUMMARY["runs"]["persistent"]
    standalone = np.asarray([a["queue_elapsed_logged_s"], a["drain_after_stream_close_s"]]) / 60
    persistent = np.asarray([b["queue_elapsed_logged_s"], b["drain_after_stream_close_s"]]) / 60
    fig, ax = plt.subplots(figsize=(11, 5.2))
    fig.suptitle("Queue processing and post-transmission drain", y=0.97)
    ypos = np.arange(2)
    ax.barh(ypos - 0.19, standalone, height=0.32, color=BLUE, label="Standalone CUDA")
    ax.barh(ypos + 0.19, persistent, height=0.32, color=ORANGE, label="Persistent CUDA")
    for i, value in enumerate(standalone):
        ax.text(value + 0.45, ypos[i] - 0.19, f"{value:.2f} min", va="center", fontsize=10, color=BLUE)
    ax.text(persistent[0] + 0.45, ypos[0] + 0.19, f"{persistent[0]:.2f} min", va="center", fontsize=10, color=ORANGE)
    ax.text(0.45, ypos[1] + 0.19, f"{b['drain_after_stream_close_s']:.3f} s", va="center", fontsize=10, color=ORANGE)
    ax.set_yticks(ypos, ("First queue event → final registration", "Stream close → final registration"))
    ax.invert_yaxis()
    ax.set(xlabel="Elapsed time (minutes)", xlim=(0, max(standalone) * 1.18))
    ax.grid(axis="x", color="#DCE3E9", lw=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", frameon=False)
    mean_speedup = SUMMARY["comparison"]["queue_timing"]["mean_speedup"]
    queue_speedup = SUMMARY["comparison"]["queue_timing"]["total_queue_speedup"]
    fig.text(0.5, 0.89, f"Mean per-registration speedup: {mean_speedup:.2f}×     •     Whole-queue speedup: {queue_speedup:.2f}×",
             ha="center", fontsize=10, color=NEUTRAL)
    fig.subplots_adjust(left=0.31, right=0.95, top=0.83, bottom=0.16)
    save(fig, "acquisition_processing_time_summary.png")


def main(argv=None):
    global HERE, RUNS, SUMMARY
    args = parse_args(argv)
    HERE = args.analysis_dir
    RUNS = {"standalone": args.standalone_dir, "persistent": args.persistent_dir}
    SUMMARY = json.loads((HERE / "comparison_summary.json").read_text())
    threshold = configured_motion_threshold(args.logs_dir)
    timing_index, times, transform_index, rotations, translations, motion_index, trace = load_and_validate(threshold)
    registration_timing(timing_index, times)
    timing_distribution(times)
    motion_parameters(motion_index, trace)
    framewise_displacement(motion_index, trace, threshold)
    transform_differences(transform_index, rotations, translations)
    acquisition_summary()
    print("Validated: 8,786 timing pairs, 8,786 matched transform pairs, 8,787 aligned motion samples; summary statistics unchanged.")


if __name__ == "__main__":
    main()

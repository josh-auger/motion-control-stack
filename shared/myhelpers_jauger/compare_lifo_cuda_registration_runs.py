r"""Compare standalone and persistent CUDA runs with FIFO_FLAG=off (LIFO).

Run from the repository root:
    python3 shared/myhelpers_jauger/compare_lifo_cuda_registration_runs.py \
        --standalone-dir data/savedData_20260921_fetal1_SMS1_TR2.5_lifo_standalone \
        --persistent-dir data/savedData_20260921_fetal1_SMS1_TR2.5_lifo_persistent \
        --analysis-dir data/analysis_20260921_lifo_cuda_comparison \
        --logs-dir data/logs_motioncontrolstack

If --logs-dir contains multiple FIFO-off logs for either execution mode, pass
--standalone-entrypoint-log and --persistent-entrypoint-log to select the pair.
Inputs are read-only; all generated tables, JSON, and figures go to --analysis-dir.

The two modes register different groups. Transform and motion comparisons use
(volume, slice group), never registration row number. FD comparisons with the
same predecessor are reported separately because FD depends on the sequence of
registered groups.
"""

import argparse
import csv
import filecmp
import hashlib
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OUT = None
LOGS = None
RUNS = {}
ENTRYPOINT = {}
IDENTITY_KEY = None
GROUP_STRIDE = None
MODES = ("standalone", "persistent")
COLORS = {"standalone": "#1768AC", "persistent": "#D05A2D"}
STAMP = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}) - ")
TFM = re.compile(r"alignTransform_(\d+)_(\d+)-(\d+)(?:_identity)?\.tfm$")
POINTER = re.compile(r"volume_(\d+)_group_(\d+)\.txt$")
METRIC = re.compile(r"\[TIMING\] (\w+): ([+-]?[\d.]+)")
REF_METRIC = re.compile(r"\[REF_SETUP_TIMING\] (\w+): ([+-]?[\d.]+)")
MOTION_COLUMNS = ("X_rotation(rad)", "Y_rotation(rad)", "Z_rotation(rad)",
                  "X_translation(mm)", "Y_translation(mm)", "Z_translation(mm)")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.titlesize": 11, "figure.titlesize": 14,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.facecolor": "white"})


def timestamp(line):
    match = STAMP.match(line)
    return datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S,%f") if match else None


def select_entrypoint(parser, args, mode):
    explicit = getattr(args, f"{mode}_entrypoint_log")
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.is_file():
            parser.error(f"--{mode}-entrypoint-log is not a file: {path}")
        return path
    candidates = []
    for path in args.logs_dir.glob("log_entrypoint*.log"):
        content = path.read_text(errors="replace")
        if (re.search(r"^\s*REG_ENGINE=cuda\s*$", content, re.MULTILINE)
                and re.search(r"^\s*FIFO_FLAG=off\s*$", content, re.MULTILINE)
                and re.search(rf"^\s*CUDA_EXECUTION_MODE={mode}\s*$", content, re.MULTILINE)):
            candidates.append(path)
    if len(candidates) != 1:
        parser.error(f"Found {len(candidates)} FIFO-off CUDA {mode} entrypoint logs in "
                     f"{args.logs_dir}; pass --{mode}-entrypoint-log explicitly")
    return candidates[0]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--standalone-dir", type=Path, required=True,
                        help="Saved acquisition directory for standalone CUDA")
    parser.add_argument("--persistent-dir", type=Path, required=True,
                        help="Saved acquisition directory for persistent CUDA")
    parser.add_argument("--analysis-dir", type=Path, required=True,
                        help="Directory for comparison CSV, JSON, and figure outputs")
    parser.add_argument("--logs-dir", type=Path, required=True,
                        help="Directory containing entrypoint and queue-processor stdout logs")
    for mode in MODES:
        parser.add_argument(f"--{mode}-entrypoint-log", type=Path,
                            help=f"Specific {mode} entrypoint log if --logs-dir contains multiple runs")
    args = parser.parse_args(argv)
    for name in ("standalone_dir", "persistent_dir", "logs_dir"):
        path = getattr(args, name).expanduser().resolve()
        if not path.is_dir():
            parser.error(f"--{name.replace('_', '-')} is not a directory: {path}")
        setattr(args, name, path)
    args.analysis_dir = args.analysis_dir.expanduser().resolve()
    if args.standalone_dir == args.persistent_dir:
        parser.error("--standalone-dir and --persistent-dir must differ")
    if any(args.analysis_dir.is_relative_to(path)
           for path in (args.standalone_dir, args.persistent_dir)):
        parser.error("--analysis-dir must be outside both acquisition directories")
    args.entrypoint = {mode: select_entrypoint(parser, args, mode) for mode in MODES}
    return args


def acquisition_position(volume, group):
    return volume * GROUP_STRIDE + group


def one_log(directory, pattern, marker):
    paths = [p for p in directory.glob(pattern) if marker in p.read_text(errors="replace")]
    if len(paths) != 1:
        raise ValueError(f"Expected one {pattern} containing {marker!r} in {directory}; found {paths}")
    return paths[0]


def protocol_fields(line):
    return dict(part.split("=", 1) for part in
                line.split("SLIMM_PROTOCOL\tSUCCESS\t", 1)[1].strip().split("\t")
                if "=" in part)


def parse_queue(mode):
    directory = RUNS[mode]
    path = one_log(directory, "log_local_queue_processor*", "Running registration call")
    records, skipped, skip_times, reference_events, loads, finalizations, failures = [], [], [], [], [], [], []
    start = reset = last_processed = None
    current = None
    last_elapsed_s = None
    for line in path.open(encoding="utf-8", errors="replace"):
        ts = timestamp(line)
        if "Started monitoring at :" in line:
            start = ts
        if "Provisional reference volume set to :" in line:
            reference_events.append({"kind": "initial", "time": str(ts),
                                     "path": line.rsplit(" : ", 1)[1].strip()})
        if "Assigning new provisional reference volume :" in line:
            reference_events.append({"kind": "update", "time": str(ts),
                                     "path": line.rsplit(" : ", 1)[1].strip()})
        if "Reference volume calibration successful" in line:
            reference_events.append({"kind": "calibration_accepted", "time": str(ts)})
        if "Transform file not found: /data/alignTransform_0001.tfm" in line:
            reference_events.append({"kind": "calibration_missing_transform", "time": str(ts)})
        if "Old file found. Skipping" in line:
            match = POINTER.search(line.strip())
            if not match:
                raise ValueError(f"Cannot identify skipped group: {line}")
            skipped.append((int(match.group(1)), int(match.group(2))))
            skip_times.append(ts)
        if "Registration FAILED" in line or "CUDA registration FAILED" in line:
            failures.append(line.strip())
        match = re.search(r"Running registration call (\d+)", line)
        if match:
            if current is not None:
                raise ValueError(f"Registration {current['index']} has no elapsed timing")
            current = {"index": int(match.group(1)), "started": ts,
                       "timing": {}, "ref_setup": {}}
            records.append(current)
        if current is not None:
            match = METRIC.search(line)
            if match:
                current["timing"][match.group(1)] = float(match.group(2))
            match = REF_METRIC.search(line)
            if match:
                current["ref_setup"][match.group(1)] = float(match.group(2))
            if "SLIMM REGISTER " in line and "completed in" in line:
                current["protocol"] = protocol_fields(line)
            if "Alignment transform saved as /data/alignTransform_" in line:
                match = TFM.search(line.strip())
                if not match:
                    raise ValueError(f"Cannot identify transform: {line}")
                current["transform_filename"] = match.group(0)
                current["key"] = (int(match.group(2)), int(match.group(3)))
            if "CUDA registration call elapsed runtime (sec) :" in line:
                current["elapsed_ms"] = float(line.rsplit(" : ", 1)[1]) * 1000
                current["finished"] = ts
                current = None
        if "SLIMM LOAD_REFERENCE completed" in line:
            fields = protocol_fields(line)
            fields["logged_at"] = str(ts)
            loads.append(fields)
        match = re.search(r"REG STATUS: volume (\d+) finalized \((\d+) registered, (\d+) skipped, (\d+) failed\)", line)
        if match:
            finalizations.append({"volume": int(match.group(1)), "time": ts,
                                  "registered": int(match.group(2)),
                                  "skipped": int(match.group(3)),
                                  "failed": int(match.group(4))})
        if re.search(r"Processed item \d+ in", line):
            last_processed = ts
        if "Total elapsed time (sec) :" in line:
            last_elapsed_s = float(line.rsplit(" : ", 1)[1])
        if "Reset trigger detected" in line:
            reset = ts
    if current is not None or not start or not reset or not finalizations:
        raise ValueError(f"Incomplete queue log: {path}")
    if [r["index"] for r in records] != list(range(1, len(records) + 1)):
        raise ValueError(f"Nonsequential registration calls: {path}")
    if any("elapsed_ms" not in r or "key" not in r for r in records):
        raise ValueError(f"Missing registration result in {path}")
    return {"path": path, "records": records, "skipped": skipped,
            "skip_times": skip_times,
            "reference_events": reference_events, "loads": loads,
            "finalizations": finalizations, "failures": failures,
            "start": start, "reset": reset, "last_processed": last_processed,
            "last_elapsed_s": last_elapsed_s}


def parse_fire(mode):
    path = one_log(RUNS[mode], "log_python-fire-server*", "Received MRD_MESSAGE_CLOSE")
    result = {"path": path, "timeout": False}
    for line in path.open(encoding="utf-8", errors="replace"):
        if "Received MRD_MESSAGE_CLOSE" in line:
            result["stream_close"] = timestamp(line)
        if "Receive elapsed time (sec)" in line:
            result["receive_elapsed_s"] = float(line.rsplit(" : ", 1)[1])
        if "Timed out waiting for close acknowledgement" in line:
            result["timeout"] = True
        if "Moved output files into subdirectory" in line:
            result["consolidated_at"] = timestamp(line)
    return result


def parse_status(mode):
    path = RUNS[mode] / "registration_status.json"
    data = json.loads(path.read_text())
    entries = {(int(v), int(g)): entry
               for v, item in data["volumes"].items()
               for g, entry in item["groups"].items()}
    return {"path": path, "schema_version": data["schema_version"],
            "revision": data["revision"], "entries": entries,
            "finalized_volumes": sum(x["registration_finalized"] for x in data["volumes"].values()),
            "statuses": Counter(x["status"] for x in entries.values())}


def parse_transforms(mode):
    result = {}
    for path in RUNS[mode].glob("*.tfm"):
        match = TFM.fullmatch(path.name)
        if not match:
            raise ValueError(f"Unexpected transform filename: {path}")
        data = {}
        for line in path.open():
            if line.startswith("Transform:"):
                data["type"] = line.split(":", 1)[1].strip()
            elif line.startswith("Parameters:"):
                data["parameters"] = np.asarray([float(x) for x in line.split(":", 1)[1].split()])
            elif line.startswith("FixedParameters:"):
                data["fixed_parameters"] = np.asarray([float(x) for x in line.split(":", 1)[1].split()])
        key = (int(match.group(2)), int(match.group(3)))
        if key in result:
            raise ValueError(f"Duplicate transform for {key}: {path}")
        data.update({"path": path, "filename": path.name, "index": int(match.group(1))})
        result[key] = data
    return result


def parse_motion(mode):
    paths = list(RUNS[mode].glob("motionMonitor_data*.csv"))
    if len(paths) != 1:
        raise ValueError(f"Expected one motion CSV in {RUNS[mode]}: {paths}")
    with paths[0].open(newline="") as source:
        reader = csv.DictReader(source)
        rows = list(reader)
        columns = reader.fieldnames
    by_key = {(int(r["Volume_index"]), int(r["Slice_group_index"])): r for r in rows}
    if len(by_key) != len(rows):
        raise ValueError(f"Duplicate motion key in {paths[0]}")
    if [int(r["reg_index"]) for r in rows] != list(range(1, len(rows) + 1)):
        raise ValueError(f"Motion rows are not in registration order: {paths[0]}")
    predecessor = {key: prev for key, prev in zip(by_key, [None, *list(by_key)[:-1]])}
    return {"path": paths[0], "rows": rows, "by_key": by_key,
            "predecessor": predecessor, "columns": columns}


def parse_config(mode):
    path = ENTRYPOINT[mode]
    content = path.read_text()
    keys = ("REG_ENGINE", "CUDA_EXECUTION_MODE", "FIFO_FLAG", "REG_TYPE",
            "MOCO_FLAG", "STREAM_FLAG", "HEAD_RADIUS", "MOTION_THRESH",
            "SEND_DASHBOARD_FLAG")
    config = {}
    for key in keys:
        match = re.search(rf"^\s*{key}=(\S+)\s*$", content, re.MULTILINE)
        if not match:
            raise ValueError(f"Missing {key} in {path}")
        config[key] = match.group(1)
    if (config["REG_ENGINE"], config["CUDA_EXECUTION_MODE"], config["FIFO_FLAG"]) != ("cuda", mode, "off"):
        raise ValueError(f"Unexpected configuration in {path}: {config}")
    return config


def persistent_lifecycle():
    match = re.fullmatch(r"log_entrypoint_(\d{8}_\d{6})\.log", ENTRYPOINT["persistent"].name)
    path = LOGS / f"STDOUT_queue_processor_{match.group(1)}.log" if match else None
    if path is None or not path.is_file():
        return {"stdout_log": None, "start_pids": [], "stop_pids": [],
                "ready_times_s": [], "same_pid_on_shutdown": None}
    content = path.read_text()
    started = re.findall(r"Persistent CUDA registration process started \(PID=(\d+)\)", content)
    stopped = re.findall(r"Persistent CUDA registration process stopped \(PID=(\d+)\)", content)
    ready = re.findall(r"SLIMM persistent registration ready \(([\d.]+) s\)", content)
    return {"stdout_log": str(path), "start_pids": started, "stop_pids": stopped,
            "ready_times_s": [float(x) for x in ready],
            "same_pid_on_shutdown": bool(started) and started == stopped}


def stats(values):
    a = np.asarray(values, dtype=float)
    if not len(a):
        return None
    return {"n": len(a), "min": float(np.min(a)), "max": float(np.max(a)),
            "mean": float(np.mean(a)), "median": float(np.median(a)),
            "std_population": float(np.std(a)), "p5": float(np.percentile(a, 5)),
            "p95": float(np.percentile(a, 95))}


def diff_stats(values):
    a = np.asarray(values, dtype=float)
    if not len(a):
        return None
    return {"n": len(a), "mean_abs": float(np.mean(np.abs(a))),
            "max_abs": float(np.max(np.abs(a))),
            "rms": float(np.sqrt(np.mean(a * a))),
            "p95_abs": float(np.percentile(np.abs(a), 95))}


def write_csv(name, rows, columns):
    with (OUT / name).open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def plot_and_save(fig, name):
    fig.savefig(OUT / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_figures(queue, status, motion, transform_rows, summary):
    # The common x coordinate is acquisition order. LIFO registers sparse,
    # different subsets; registration call number would suggest false pairing.
    fig, (ax, detail) = plt.subplots(2, 1, figsize=(13, 7), sharex=True,
                                      gridspec_kw={"height_ratios": [2, 1]})
    for mode in MODES:
        records = queue[mode]["records"]
        x = np.asarray([acquisition_position(*r["key"]) for r in records])
        y = np.asarray([r["elapsed_ms"] for r in records])
        ax.plot(x, y, ".", ms=2.3, alpha=0.65, color=COLORS[mode],
                label=f"{mode.title()} CUDA (n={len(y):,})")
        ax.axhline(y.mean(), color=COLORS[mode], ls="--", lw=1.4,
                   label=f"{mode.title()} mean {y.mean():.1f} ms")
        if mode == "persistent":
            detail.plot(x, y, ".", ms=2.3, alpha=0.7, color=COLORS[mode])
            detail.axhline(y.mean(), color=COLORS[mode], ls="--", lw=1.3)
            first_load = any(r["started"] <= datetime.fromisoformat(load["logged_at"]) <= r["finished"]
                             for r in records[:1] for load in queue[mode]["loads"])
            note = "\nincludes reference load" if first_load else ""
            ax.annotate(f"First persistent call {y[0]:.1f} ms{note}",
                        xy=(x[0], y[0]), xytext=(x.max() * 0.13, 145),
                        arrowprops={"arrowstyle": "->", "color": COLORS[mode]}, fontsize=9)
    ax.set(ylabel="Queue registration wall time (ms)", ylim=(0, None))
    persistent_values = [r["elapsed_ms"] for r in queue["persistent"]["records"]]
    steady_values = persistent_values[1:] if len(persistent_values) > 1 else persistent_values
    detail_limit = max(9, float(np.percentile(steady_values, 99.5)) * 1.2)
    detail.set(xlabel=f"Acquisition group position (volume × {GROUP_STRIDE} + group)",
               ylabel="Persistent detail (ms)", ylim=(0, detail_limit))
    fig.suptitle("Registration time across the FIFO-off acquisition")
    ax.grid(axis="y", color="#DFE5EA")
    detail.grid(axis="y", color="#DFE5EA")
    ax.legend(frameon=False, ncol=2, fontsize=9)
    fig.tight_layout()
    plot_and_save(fig, "registration_timing_vs_index.png")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.8))
    vals = [[r["elapsed_ms"] for r in queue[m]["records"]] for m in MODES]
    for ax, mode, values, label in zip(axes,
                                       ("standalone", "persistent", "persistent"),
                                       (vals[0], vals[1], vals[1]),
                                       ("Standalone", "Persistent: full", "Persistent: detail")):
        box = ax.boxplot([values], tick_labels=[f"{label}\nn={len(values):,}"],
                         showfliers=True, showmeans=True, patch_artist=True,
                         flierprops={"marker": ".", "markersize": 2, "alpha": 0.3})
        box["boxes"][0].set(facecolor=COLORS[mode], alpha=0.42)
        ax.set(ylabel="Queue registration wall time (ms)", ylim=(0, None))
        ax.grid(axis="y", color="#DFE5EA")
        ax.text(0.58, 0.88, f"Mean {np.mean(values):.2f} ms\nMedian {np.median(values):.2f} ms",
                transform=ax.transAxes, fontsize=9)
    axes[2].set_ylim(0, detail_limit)
    axes[2].set_title(f"Detail: values above {detail_limit:.1f} ms outside view")
    fig.suptitle("Registration time distributions (full ranges and persistent detail)")
    fig.tight_layout()
    plot_and_save(fig, "registration_timing_distribution.png")

    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
    for ax, column in zip(axes.flat, MOTION_COLUMNS):
        for mode in MODES:
            rows = motion[mode]["rows"]
            x = [acquisition_position(int(r["Volume_index"]), int(r["Slice_group_index"])) for r in rows]
            y = [float(r[column]) for r in rows]
            ax.plot(x, y, ".", ms=2.1, color=COLORS[mode], alpha=0.65,
                    label=mode.title())
        ax.set_title(column)
        ax.grid(axis="y", color="#E4E9ED")
    for ax in axes[1]:
        ax.set_xlabel("Acquisition group position")
    axes[0, 0].set_ylabel("Rotation (rad)")
    axes[1, 0].set_ylabel("Translation (mm)")
    axes[0, 0].legend(frameon=False, markerscale=2)
    fig.suptitle("Rigid-body motion traces over acquired groups")
    fig.tight_layout()
    plot_and_save(fig, "motion_parameter_trace_comparison.png")

    fig, ax = plt.subplots(figsize=(13, 5))
    for mode in MODES:
        rows = motion[mode]["rows"]
        x = [acquisition_position(int(r["Volume_index"]), int(r["Slice_group_index"])) for r in rows]
        y = [float(r["Displacement(mm)"]) for r in rows]
        ax.plot(x, y, ".", ms=2.2, color=COLORS[mode], alpha=0.6,
                label=f"{mode.title()}: mean {np.mean(y):.3f}, max {np.max(y):.3f} mm")
    ax.axhline(summary["threshold_mm"], color="#52606D", ls="--", lw=1.2,
               label=f"Motion threshold {summary['threshold_mm']:.1f} mm")
    ax.set(xlabel="Acquisition group position", ylabel="Framewise displacement (mm)",
           title="Framewise displacement on each mode's registered sequence", ylim=(0, None))
    ax.grid(axis="y", color="#E4E9ED")
    ax.legend(frameon=False)
    fig.tight_layout()
    plot_and_save(fig, "framewise_displacement_comparison.png")

    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
    x = [acquisition_position(r["volume"], r["group"]) for r in transform_rows]
    for i, (ax, names, unit) in enumerate(((axes[0], ("rx", "ry", "rz"), "rad"),
                                            (axes[1], ("tx", "ty", "tz"), "mm"))):
        for name, color in zip(names, ("#315B8A", "#A04766", "#16866A")):
            ax.plot(x, [abs(r[f"d_{name}"]) for r in transform_rows], ".", ms=3,
                    color=color, alpha=0.72, label=name)
        ax.set_ylabel(f"Absolute difference ({unit})")
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="#E4E9ED")
        ax.legend(frameon=False, ncol=3)
    axes[1].set_xlabel("Acquisition group position")
    fig.suptitle(f"Transform agreement for {len(transform_rows):,} shared CUDA-registered groups")
    fig.tight_layout()
    plot_and_save(fig, "transform_parameter_differences.png")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3))
    labels = ["Queue monitoring to final group status", "Stream close to final group status"]
    a = [summary["runs"]["standalone"]["queue_to_final_status_s"],
         summary["runs"]["standalone"]["drain_to_final_status_s"]]
    b = [summary["runs"]["persistent"]["queue_to_final_status_s"],
         summary["runs"]["persistent"]["drain_to_final_status_s"]]
    for ax, label, pair, scale, unit in zip(axes, labels, ((a[0], b[0]), (a[1], b[1])),
                                            (1, 1000), ("s", "ms")):
        for i, (mode, value) in enumerate(zip(MODES, pair)):
            ax.barh(i, value * scale, color=COLORS[mode], height=0.55)
            ax.text(value * scale * 1.01, i, f"{value * scale:.1f} {unit}", va="center", fontsize=9)
        ax.set(yticks=(0, 1), yticklabels=("Standalone", "Persistent"),
               xlabel=f"Elapsed time ({unit})", title=label,
               xlim=(0, max(pair) * scale * 1.18))
        ax.invert_yaxis()
        ax.grid(axis="x", color="#E4E9ED")
        ax.set_axisbelow(True)
    fig.suptitle("Queue completion and post-stream drain")
    fig.tight_layout()
    plot_and_save(fig, "acquisition_processing_time_summary.png")

    fig, ax = plt.subplots(figsize=(13, 4.5))
    volume_numbers = sorted({v for v, g in status["standalone"]["entries"]
                             if (v, g) != IDENTITY_KEY})
    for mode in MODES:
        counts = Counter(v for (v, g), item in status[mode]["entries"].items()
                         if (v, g) != IDENTITY_KEY and item["status"] == "registered")
        ax.plot(volume_numbers, [counts[v] for v in volume_numbers],
                color=COLORS[mode], lw=1.3, alpha=0.82, label=mode.title())
    ax.set(xlabel="Volume", ylabel=f"Registered slice groups (of {GROUP_STRIDE})",
           title="Registration completeness by volume", ylim=(0, GROUP_STRIDE))
    ax.grid(axis="y", color="#E4E9ED")
    ax.legend(frameon=False)
    fig.tight_layout()
    plot_and_save(fig, "registration_completeness_by_volume.png")


def main(argv=None):
    global OUT, LOGS, RUNS, ENTRYPOINT, IDENTITY_KEY, GROUP_STRIDE
    args = parse_args(argv)
    OUT = args.analysis_dir
    LOGS = args.logs_dir
    RUNS = {"standalone": args.standalone_dir, "persistent": args.persistent_dir}
    ENTRYPOINT = args.entrypoint
    OUT.mkdir(parents=True, exist_ok=True)
    config = {m: parse_config(m) for m in MODES}
    config_diff = {k: [config[m][k] for m in MODES]
                   for k in config["standalone"]
                   if config["standalone"][k] != config["persistent"][k]}
    if config_diff != {"CUDA_EXECUTION_MODE": ["standalone", "persistent"]}:
        raise ValueError(f"Unexpected configuration differences: {config_diff}")
    threshold = float(config["standalone"]["MOTION_THRESH"])
    queue = {m: parse_queue(m) for m in MODES}
    fire = {m: parse_fire(m) for m in MODES}
    status = {m: parse_status(m) for m in MODES}
    transforms = {m: parse_transforms(m) for m in MODES}
    motion = {m: parse_motion(m) for m in MODES}

    keys = set(status["standalone"]["entries"])
    if keys != set(status["persistent"]["entries"]):
        raise ValueError("The two status files have different acquisition group keys")
    identities = {mode: [key for key, item in status[mode]["entries"].items()
                         if item.get("transform_filename", "").endswith("_identity.tfm")]
                  for mode in MODES}
    if any(len(identities[mode]) != 1 for mode in MODES) or identities["standalone"] != identities["persistent"]:
        raise ValueError(f"Expected the same single identity reference group in both runs: {identities}")
    IDENTITY_KEY = identities["standalone"][0]
    if len(keys) == 1:
        raise ValueError("No CUDA registration opportunities in status files")
    GROUP_STRIDE = max(key[1] for key in keys if key != IDENTITY_KEY) + 1
    expected_groups = len(keys) - 1
    runs = {}
    timing_rows = []
    status_rows = []
    for key in sorted(keys):
        status_rows.append({"volume": key[0], "group": key[1],
                            "standalone_status": status["standalone"]["entries"][key]["status"],
                            "persistent_status": status["persistent"]["entries"][key]["status"]})
    for mode in MODES:
        q, f, s, t, mo = queue[mode], fire[mode], status[mode], transforms[mode], motion[mode]
        registered = {key for key, item in s["entries"].items() if item["status"] == "registered"}
        skipped = {key for key, item in s["entries"].items() if item["status"] == "skipped"}
        if len(q["records"]) != len(registered) - 1 or len(q["skipped"]) != len(skipped):
            raise ValueError(f"Queue/status count discrepancy in {mode}")
        if set(q["skipped"]) != skipped or {r["key"] for r in q["records"]} != registered - {IDENTITY_KEY}:
            raise ValueError(f"Queue/status group discrepancy in {mode}")
        if set(t) != registered or set(mo["by_key"]) != registered:
            raise ValueError(f"Transform or motion/status group discrepancy in {mode}")
        if any(s["entries"][key]["transform_filename"] != t[key]["filename"]
               for key in registered):
            raise ValueError(f"Status/transform filename discrepancy in {mode}")
        if any(r["transform_filename"] != t[r["key"]]["filename"] for r in q["records"]):
            raise ValueError(f"Queue/transform filename discrepancy in {mode}")
        pointer_keys = {(int(match.group(1)), int(match.group(2)))
                        for path in RUNS[mode].glob("*group_*.txt")
                        if (match := POINTER.search(path.name))}
        if pointer_keys != keys:
            raise ValueError(f"Pointer/status group discrepancy in {mode}")
        elapsed = [r["elapsed_ms"] for r in q["records"]]
        fd = [float(r["Displacement(mm)"]) for r in mo["rows"]]
        for r in q["records"]:
            timing_rows.append({"mode": mode, "registration_index": r["index"],
                                "volume": r["key"][0], "group": r["key"][1],
                                "acquisition_position": acquisition_position(*r["key"]),
                                "queue_elapsed_ms": r["elapsed_ms"],
                                "slimm_total_ms": r["timing"].get("total_ms", ""),
                                "slimm_optimizer_ms": (r["protocol"].get("optimizer_ms", "")
                                                       if "protocol" in r else r["timing"].get("optimizer_ms", "")),
                                "persistent_register_ms": (r["protocol"].get("request_ms", "")
                                                           if "protocol" in r else "")})
        last_final = q["finalizations"][-1]["time"]
        runs[mode] = {
            "registration_opportunities": expected_groups,
            "attempts": len(q["records"]), "completed_cuda": len(registered) - 1,
            "completed_including_identity": len(registered),
            "skipped": len(skipped), "failed": s["statuses"].get("failed", 0),
            "transforms": len(t), "motion_rows": len(mo["rows"]),
            "pointer_files": len(pointer_keys), "raw_files": len(list(RUNS[mode].glob("*.raw"))),
            "nhdr_files": len(list(RUNS[mode].glob("*.nhdr"))),
            "finalized_volumes": s["finalized_volumes"], "status_revision": s["revision"],
            "completion_rate_pct": 100 * (len(registered) - 1) / expected_groups,
            "queue_timing_ms": stats(elapsed),
            "summed_registration_wall_s": sum(elapsed) / 1000,
            "slimm_total_ms": stats([r["timing"]["total_ms"] for r in q["records"]
                                     if "total_ms" in r["timing"]]),
            "optimizer_ms": stats([float(r["protocol"]["optimizer_ms"]) if "protocol" in r
                                   else r["timing"]["optimizer_ms"] for r in q["records"]]),
            "persistent_register_ms": stats([float(r["protocol"]["request_ms"])
                                             for r in q["records"] if "protocol" in r]),
            "standalone_reference_setup_ms": stats([r["ref_setup"]["total_ms"]
                                                     for r in q["records"] if "total_ms" in r["ref_setup"]]),
            "standalone_target_read_ms": stats([r["timing"]["target_read_ms"]
                                                for r in q["records"] if "target_read_ms" in r["timing"]]),
            "standalone_binary_startup_ms": stats([r["timing"]["startup_ms"]
                                                    for r in q["records"] if "startup_ms" in r["timing"]]),
            "persistent_target_read_ms": stats([float(r["protocol"]["target_read_ms"])
                                                for r in q["records"] if "protocol" in r]),
            "first_registration_ms": elapsed[0],
            "persistent_after_first_ms": stats(elapsed[1:]) if mode == "persistent" else None,
            "queue_monitor_start": str(q["start"]), "last_registration_finished": str(q["records"][-1]["finished"]),
            "final_group_status": str(last_final), "queue_reset": str(q["reset"]),
            "stream_close": str(f["stream_close"]),
            "stream_duration_logged_s": f["receive_elapsed_s"],
            "queue_to_last_registration_s": (q["records"][-1]["finished"] - q["start"]).total_seconds(),
            "queue_to_final_status_s": (last_final - q["start"]).total_seconds(),
            "drain_to_last_registration_s": (q["records"][-1]["finished"] - f["stream_close"]).total_seconds(),
            "drain_to_final_status_s": (last_final - f["stream_close"]).total_seconds(),
            "completed_by_stream_close": sum(r["finished"] <= f["stream_close"] for r in q["records"]),
            "attempted_by_stream_close": sum(r["started"] <= f["stream_close"] for r in q["records"]),
            "skipped_by_stream_close": sum(ts <= f["stream_close"] for ts in q["skip_times"]),
            "fire_timeout": f["timeout"], "consolidated_at": str(f.get("consolidated_at")),
            "reference_events": q["reference_events"], "reference_updates": sum(
                e["kind"] == "update" for e in q["reference_events"]),
            "persistent_load_reference_count": len(q["loads"]),
            "persistent_load_details": q["loads"],
            "queue_failure_lines": q["failures"],
            "fd_mm": stats(fd), "fd_above_threshold": sum(x > threshold for x in fd),
            "motion_flag_count": sum(int(r["Motion_flag"]) for r in mo["rows"]),
        }
        runs[mode]["unresolved_at_stream_close"] = (expected_groups -
            runs[mode]["completed_by_stream_close"] - runs[mode]["skipped_by_stream_close"])

    shared = {m: {key for key, entry in status[m]["entries"].items()
                  if entry["status"] == "registered"} for m in MODES}
    both = shared["standalone"] & shared["persistent"]
    matched_cuda = sorted(both - {IDENTITY_KEY})
    transform_rows = []
    for key in matched_cuda:
        a, b = transforms["standalone"][key], transforms["persistent"][key]
        if a["type"] != b["type"] or a["type"] != "Euler3DTransform_double_3_3":
            raise ValueError(f"Transform representation mismatch at {key}")
        if not np.array_equal(a["fixed_parameters"], b["fixed_parameters"]):
            raise ValueError(f"Fixed transform parameters mismatch at {key}")
        delta = b["parameters"] - a["parameters"]
        transform_rows.append({"volume": key[0], "group": key[1],
                               "standalone_index": a["index"], "persistent_index": b["index"],
                               **{f"d_{name}": float(delta[i]) for i, name in
                                  enumerate(("rx", "ry", "rz", "tx", "ty", "tz"))}})
    rotation = np.asarray([[r[f"d_{name}"] for name in ("rx", "ry", "rz")]
                           for r in transform_rows])
    translation = np.asarray([[r[f"d_{name}"] for name in ("tx", "ty", "tz")]
                              for r in transform_rows])
    max_rotation_row = transform_rows[int(np.argmax(np.max(np.abs(rotation), axis=1)))]
    max_translation_row = transform_rows[int(np.argmax(np.max(np.abs(translation), axis=1)))]
    timing_by_key = {m: {r["key"]: r["elapsed_ms"] for r in queue[m]["records"]}
                     for m in MODES}
    matched_timing = {m: stats([timing_by_key[m][key] for key in matched_cuda])
                      for m in MODES}

    motion_rows = []
    comparable_fd_rows = []
    for key in sorted(both):
        a, b = motion["standalone"]["by_key"][key], motion["persistent"]["by_key"][key]
        row = {"volume": key[0], "group": key[1],
               "standalone_reg_index": int(a["reg_index"]),
               "persistent_reg_index": int(b["reg_index"]),
               "same_predecessor": motion["standalone"]["predecessor"][key] ==
                                   motion["persistent"]["predecessor"][key]}
        for column in (*MOTION_COLUMNS, "Displacement(mm)", "Cumulative_displacement(mm)"):
            row[f"delta_{column}"] = float(b[column]) - float(a[column])
        row["standalone_motion_flag"] = a["Motion_flag"]
        row["persistent_motion_flag"] = b["Motion_flag"]
        motion_rows.append(row)
        if row["same_predecessor"] and key != IDENTITY_KEY:
            comparable_fd_rows.append(row)

    raw_names = {m: {p.name: p for p in RUNS[m].glob("*.raw")} for m in MODES}
    shared_raw_names = set(raw_names["standalone"]) & set(raw_names["persistent"])
    differing_raw = [name for name in sorted(shared_raw_names)
                     if not filecmp.cmp(raw_names["standalone"][name],
                                        raw_names["persistent"][name], shallow=False)]
    upsampled_hashes = {}
    for mode in MODES:
        paths = list(RUNS[mode].glob("*_upsampled.raw"))
        if len(paths) > 1:
            raise ValueError(f"Multiple upsampled reference raw files in {RUNS[mode]}: {paths}")
        upsampled_hashes[mode] = sha256(paths[0]) if paths else None
    speedup = runs["standalone"]["queue_timing_ms"]["mean"] / runs["persistent"]["queue_timing_ms"]["mean"]
    queue_a, queue_b = [runs[m]["queue_to_final_status_s"] for m in MODES]
    fd_delta = [r["delta_Displacement(mm)"] for r in comparable_fd_rows]
    max_fd_row = max(motion_rows, key=lambda r: abs(r["delta_Displacement(mm)"]))
    first_status_difference = next((key for key in sorted(keys)
        if status["standalone"]["entries"][key]["status"] !=
           status["persistent"]["entries"][key]["status"]), None)
    summary = {
        "sources": {m: {"directory": str(RUNS[m]), "entrypoint": str(ENTRYPOINT[m]),
                        "queue_log": str(queue[m]["path"]), "fire_log": str(fire[m]["path"]),
                        "status": str(status[m]["path"]), "motion": str(motion[m]["path"])}
                    for m in MODES},
        "configuration": config, "configuration_differences": config_diff,
        "threshold_mm": threshold, "runs": runs,
        "persistent_lifecycle": persistent_lifecycle(),
        "raw_inputs": {"shared_raw_names": len(shared_raw_names),
                       "standalone_only": len(set(raw_names["standalone"]) - shared_raw_names),
                       "persistent_only": len(set(raw_names["persistent"]) - shared_raw_names),
                       "differing_shared_raw_content": len(differing_raw),
                       "first_differences": differing_raw[:10],
                       "upsampled_reference_sha256": upsampled_hashes},
        "comparison": {
            "mean_registration_speedup": speedup,
            "mean_registration_reduction_ms": runs["standalone"]["queue_timing_ms"]["mean"] - runs["persistent"]["queue_timing_ms"]["mean"],
            "mean_registration_reduction_pct": 100 * (1 - 1 / speedup),
            "matched_target_timing_ms": matched_timing,
            "matched_target_mean_speedup": (matched_timing["standalone"]["mean"] /
                                            matched_timing["persistent"]["mean"]),
            "queue_to_final_status_saved_s": queue_a - queue_b,
            "queue_to_final_status_reduction_pct": 100 * (1 - queue_b / queue_a),
            "queue_to_final_status_speedup": queue_a / queue_b,
            "shared_registered_including_identity": len(both),
            "shared_cuda_registered": len(matched_cuda),
            "standalone_only_registered": len(shared["standalone"] - both),
            "persistent_only_registered": len(shared["persistent"] - both),
            "status_same_keys": keys == set(status["persistent"]["entries"]),
            "status_difference_count": sum(status["standalone"]["entries"][key]["status"] !=
                                           status["persistent"]["entries"][key]["status"] for key in keys),
            "first_status_difference": {"volume": first_status_difference[0],
                                        "group": first_status_difference[1]} if first_status_difference else None,
            "transform": {"representation": "Euler3DTransform_double_3_3",
                          "matched_cuda": len(transform_rows),
                          "rotation_rad": diff_stats(rotation.ravel()),
                          "translation_mm": diff_stats(translation.ravel()),
                          "largest_rotation_group": {"volume": max_rotation_row["volume"],
                                                     "group": max_rotation_row["group"],
                                                     "max_abs_rad": float(np.max(np.abs([max_rotation_row[f"d_{x}"] for x in ("rx", "ry", "rz")]))),
                                                     "standalone_index": max_rotation_row["standalone_index"],
                                                     "persistent_index": max_rotation_row["persistent_index"]},
                          "largest_translation_group": {"volume": max_translation_row["volume"],
                                                        "group": max_translation_row["group"],
                                                        "max_abs_mm": float(np.max(np.abs([max_translation_row[f"d_{x}"] for x in ("tx", "ty", "tz")]))),
                                                        "standalone_index": max_translation_row["standalone_index"],
                                                        "persistent_index": max_translation_row["persistent_index"]}},
            "motion": {"shared_samples_including_identity": len(motion_rows),
                       "parameter_differences": {col: diff_stats([r[f"delta_{col}"] for r in motion_rows])
                                                 for col in MOTION_COLUMNS},
                       "fd_direct_matched_key_difference_mm": diff_stats(
                           [r["delta_Displacement(mm)"] for r in motion_rows]),
                       "fd_same_predecessor_samples": len(comparable_fd_rows),
                       "fd_same_predecessor_difference_mm": diff_stats(fd_delta),
                       "largest_direct_fd_difference": {"volume": max_fd_row["volume"],
                                                        "group": max_fd_row["group"],
                                                        "difference_mm": max_fd_row["delta_Displacement(mm)"],
                                                        "same_predecessor": max_fd_row["same_predecessor"]},
                       "motion_flag_disagreements_on_shared_keys": sum(
                           r["standalone_motion_flag"] != r["persistent_motion_flag"] for r in motion_rows)}}}

    write_csv("registration_timing_comparison.csv", timing_rows,
              ["mode", "registration_index", "volume", "group", "acquisition_position",
               "queue_elapsed_ms", "slimm_total_ms", "slimm_optimizer_ms", "persistent_register_ms"])
    write_csv("registration_status_comparison.csv", status_rows,
              ["volume", "group", "standalone_status", "persistent_status"])
    write_csv("transform_difference_summary.csv", transform_rows,
              ["volume", "group", "standalone_index", "persistent_index",
               "d_rx", "d_ry", "d_rz", "d_tx", "d_ty", "d_tz"])
    write_csv("motion_trace_difference_summary.csv", motion_rows,
              ["volume", "group", "standalone_reg_index", "persistent_reg_index",
               "same_predecessor", *[f"delta_{c}" for c in (*MOTION_COLUMNS,
               "Displacement(mm)", "Cumulative_displacement(mm)")],
               "standalone_motion_flag", "persistent_motion_flag"])
    table_metrics = (
        ("Registration opportunities", "registration_opportunities"),
        ("Attempted CUDA registrations", "attempts"),
        ("Completed CUDA registrations", "completed_cuda"),
        ("Skipped groups", "skipped"), ("Failed groups", "failed"),
        ("Transform files including identity", "transforms"),
        ("CUDA completion rate (%)", "completion_rate_pct"),
        ("Mean registration time (ms)", ("queue_timing_ms", "mean")),
        ("Median registration time (ms)", ("queue_timing_ms", "median")),
        ("P95 registration time (ms)", ("queue_timing_ms", "p95")),
        ("Max registration time (ms)", ("queue_timing_ms", "max")),
        ("Queue monitoring to final status (s)", "queue_to_final_status_s"),
        ("Stream close to final status (s)", "drain_to_final_status_s"),
        ("Unresolved at stream close", "unresolved_at_stream_close"),
        ("Reference updates", "reference_updates"),
        ("Mean FD (mm)", ("fd_mm", "mean")),
        ("Max FD (mm)", ("fd_mm", "max")),
        (f"FD over {threshold:g} mm", "fd_above_threshold"),
    )
    table = []
    for label, field in table_metrics:
        table.append({"metric": label, **{mode: runs[mode][field] if isinstance(field, str)
                      else runs[mode][field[0]][field[1]] for mode in MODES}})
    table.extend((
        {"metric": "Shared CUDA transforms compared", "standalone": len(transform_rows),
         "persistent": len(transform_rows)},
        {"metric": "Maximum matched rotation difference (rad)", "standalone": "",
         "persistent": summary["comparison"]["transform"]["rotation_rad"]["max_abs"]},
        {"metric": "Maximum matched translation difference (mm)", "standalone": "",
         "persistent": summary["comparison"]["transform"]["translation_mm"]["max_abs"]},
    ))
    write_csv("lifo_cuda_comparison_summary.csv", table, ["metric", *MODES])
    (OUT / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_figures(queue, status, motion, transform_rows, summary)
    print(json.dumps({"runs": runs, "comparison": summary["comparison"],
                      "raw_inputs": summary["raw_inputs"]}, indent=2))


if __name__ == "__main__":
    main()

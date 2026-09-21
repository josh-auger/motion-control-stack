"""Compare the saved FIFO-on standalone and persistent CUDA acquisitions.

Run from the repository root:
    python3 shared/myhelpers_jauger/compare_fifo_cuda_registration_runs.py \
        --standalone-dir data/savedData_20260921_fetal1_SMS1_TR2.5_fifo_standalone \
        --persistent-dir data/savedData_20260921_fetal1_SMS1_TR2.5_fifo_persistent \
        --analysis-dir data/analysis_20260921_fifo_cuda_comparison

The acquisition directories are read-only inputs. Results go to --analysis-dir.
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

import numpy as np


OUT = None
RUNS = {}
TIMESTAMP = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}) - ")
POINTER = re.compile(r"volume_(\d{4})_group_(\d{4})\.txt$")
TRANSFORM = re.compile(r"alignTransform_(\d{4})_(\d{4})-(\d{4})(?:_identity)?\.tfm$")
METRIC = re.compile(r"\[TIMING\] (\w+): ([+-]?[\d.]+)")
REF_METRIC = re.compile(r"\[REF_SETUP_TIMING\] (\w+): ([+-]?[\d.]+)")


def stamp(line):
    match = TIMESTAMP.match(line)
    return datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S,%f") if match else None


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--standalone-dir", type=Path, required=True,
                        help="Saved acquisition directory for standalone CUDA")
    parser.add_argument("--persistent-dir", type=Path, required=True,
                        help="Saved acquisition directory for persistent CUDA")
    parser.add_argument("--analysis-dir", type=Path, required=True,
                        help="Directory for comparison CSV, JSON, and figure outputs")
    args = parser.parse_args(argv)
    for mode in ("standalone", "persistent"):
        path = getattr(args, f"{mode}_dir").expanduser().resolve()
        if not path.is_dir():
            parser.error(f"--{mode}-dir is not a directory: {path}")
        setattr(args, f"{mode}_dir", path)
    args.analysis_dir = args.analysis_dir.expanduser().resolve()
    if args.standalone_dir == args.persistent_dir:
        parser.error("--standalone-dir and --persistent-dir must differ")
    if any(args.analysis_dir.is_relative_to(path)
           for path in (args.standalone_dir, args.persistent_dir)):
        parser.error("--analysis-dir must be outside both acquisition directories")
    return args


def stats(values):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return None
    return {
        "n": len(values), "min": float(np.min(values)), "max": float(np.max(values)),
        "mean": float(np.mean(values)), "median": float(np.median(values)),
        "std_population": float(np.std(values, ddof=0)),
        "p5": float(np.percentile(values, 5)), "p95": float(np.percentile(values, 95)),
    }


def differences(values):
    values = np.asarray(values, dtype=float)
    return {
        "n": len(values), "mean_abs": float(np.mean(np.abs(values))),
        "max_abs": float(np.max(np.abs(values))),
        "rms": float(np.sqrt(np.mean(values ** 2))),
        "p95_abs": float(np.percentile(np.abs(values), 95)),
    }


def protocol_fields(line):
    protocol = line.split("SLIMM_PROTOCOL\tSUCCESS\t", 1)[1].strip()
    return dict(field.split("=", 1) for field in protocol.split("\t") if "=" in field)


def read_queue(mode):
    path = next(RUNS[mode].glob("log_local_queue_processor*"))
    records = []
    loads = []
    current = None
    started = finalized = final_processed = reset = None
    elapsed_summary = None
    reference_events = []
    failures = []
    for line in path.open(encoding="utf-8", errors="replace"):
        ts = stamp(line)
        if "Started monitoring at :" in line:
            started = ts
        elif "Provisional reference volume set to :" in line:
            reference_events.append({"time": str(ts), "kind": "initial", "path": line.split(" : ", 2)[-1].strip()})
        elif "Assigning new provisional reference volume :" in line:
            reference_events.append({"time": str(ts), "kind": "update", "path": line.split(" : ", 2)[-1].strip()})
        elif "Reference volume calibration successful" in line:
            reference_events.append({"time": str(ts), "kind": "accepted"})
        elif "Transform file not found: /data/alignTransform_0001.tfm" in line:
            reference_events.append({"time": str(ts), "kind": "calibration_missing_transform"})
        if "Registration FAILED" in line or "CUDA registration FAILED" in line:
            failures.append(line.strip())

        match = re.search(r"Running registration call (\d+)", line)
        if match:
            current = {"index": int(match.group(1)), "started": ts, "timing": {}, "ref_setup": {}}
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
            if "CUDA registration call elapsed runtime (sec) :" in line:
                current["elapsed_ms"] = float(line.rsplit(" : ", 1)[1]) * 1000
                current["finished"] = ts
                current = None
        if "SLIMM LOAD_REFERENCE completed" in line:
            fields = protocol_fields(line)
            fields["logged_at"] = str(ts)
            loads.append(fields)
        if "REG STATUS: volume 191 finalized" in line:
            finalized = ts
        if "Processed item 8788 in" in line:
            final_processed = ts
        if "Total elapsed time (sec) :" in line:
            elapsed_summary = float(line.rsplit(" : ", 1)[1])
        if "Reset trigger detected" in line:
            reset = ts
    return {
        "path": str(path), "records": records, "loads": loads,
        "started": started, "finalized": finalized, "final_processed": final_processed,
        "reset": reset, "elapsed_summary_s": elapsed_summary,
        "reference_events": reference_events, "failures": failures,
    }


def read_fire(mode):
    path = next(RUNS[mode].glob("log_python-fire-server*"))
    result = {"path": str(path), "timeout": False}
    for line in path.open(encoding="utf-8", errors="replace"):
        if "Received MRD_MESSAGE_CLOSE" in line:
            result["stream_close"] = stamp(line)
        if "Receive elapsed time (sec)" in line:
            result["receive_elapsed_s"] = float(line.rsplit(" : ", 1)[1])
        if "Timed out waiting for close acknowledgement" in line:
            result["timeout"] = True
        if "Moved output files into subdirectory" in line:
            result["consolidated_at"] = stamp(line)
    return result


def read_status(mode):
    path = RUNS[mode] / "registration_status.json"
    status = json.loads(path.read_text())
    entries = {}
    for volume, volume_data in status["volumes"].items():
        for group, entry in volume_data["groups"].items():
            entries[(int(volume), int(group))] = entry
    return {"path": str(path), "schema_version": status["schema_version"],
            "revision": status["revision"], "volumes": status["volumes"], "entries": entries,
            "statuses": dict(Counter(entry["status"] for entry in entries.values())),
            "finalized_volumes": sum(bool(v["registration_finalized"]) for v in status["volumes"].values())}


def read_transforms(mode):
    transforms = {}
    for path in RUNS[mode].glob("*.tfm"):
        match = TRANSFORM.match(path.name)
        assert match, path
        fields = {}
        for line in path.open():
            if line.startswith("Transform:"):
                fields["type"] = line.split(":", 1)[1].strip()
            elif line.startswith("Parameters:"):
                fields["parameters"] = [float(x) for x in line.split(":", 1)[1].split()]
            elif line.startswith("FixedParameters:"):
                fields["fixed_parameters"] = [float(x) for x in line.split(":", 1)[1].split()]
        fields["volume"] = int(match.group(2))
        fields["group"] = int(match.group(3))
        fields["index"] = int(match.group(1))
        transforms[path.name] = fields
    return transforms


def read_motion(mode):
    path = next(RUNS[mode].glob("motionMonitor_data*.csv"))
    with path.open(newline="") as source:
        reader = csv.DictReader(source)
        rows = list(reader)
        columns = reader.fieldnames
    keyed = {(int(row["Volume_index"]), int(row["Slice_group_index"])): row for row in rows}
    assert len(keyed) == len(rows)
    return {"path": str(path), "columns": columns, "rows": rows, "keyed": keyed}


def write_csv(name, rows, fields):
    with (OUT / name).open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv=None):
    global OUT, RUNS
    args = parse_args(argv)
    OUT = args.analysis_dir
    RUNS = {"standalone": args.standalone_dir, "persistent": args.persistent_dir}
    OUT.mkdir(parents=True, exist_ok=True)
    queue = {mode: read_queue(mode) for mode in RUNS}
    fire = {mode: read_fire(mode) for mode in RUNS}
    status = {mode: read_status(mode) for mode in RUNS}
    transforms = {mode: read_transforms(mode) for mode in RUNS}
    motion = {mode: read_motion(mode) for mode in RUNS}

    timings = []
    for index in range(1, 8787):
        a = queue["standalone"]["records"][index - 1]
        b = queue["persistent"]["records"][index - 1]
        assert a["index"] == b["index"] == index
        timings.append({"registration_index": index, "standalone_ms": a["elapsed_ms"],
                        "persistent_ms": b["elapsed_ms"],
                        "persistent_register_ms": float(b["protocol"]["request_ms"]),
                        "standalone_slimm_total_ms": a["timing"].get("total_ms"),
                        "standalone_optimizer_ms": a["timing"].get("optimizer_ms"),
                        "persistent_optimizer_ms": float(b["protocol"]["optimizer_ms"])})
    write_csv("registration_timing_comparison.csv", timings, list(timings[0]))

    names_a = set(transforms["standalone"])
    names_b = set(transforms["persistent"])
    matched = sorted(names_a & names_b)
    transform_rows = []
    for name in matched:
        a = transforms["standalone"][name]
        b = transforms["persistent"][name]
        assert a["type"] == b["type"]
        assert len(a["parameters"]) == len(b["parameters"]) == 6
        if "identity" in name:
            continue
        diffs = np.asarray(b["parameters"]) - np.asarray(a["parameters"])
        transform_rows.append({"transform": name, "index": a["index"],
                               "volume": a["volume"], "group": a["group"],
                               **{f"d_{axis}": float(diffs[i]) for i, axis in enumerate(("rx", "ry", "rz", "tx", "ty", "tz"))},
                               "max_abs_rotation_rad": float(np.max(np.abs(diffs[:3]))),
                               "max_abs_translation_mm": float(np.max(np.abs(diffs[3:])))})
    write_csv("transform_difference_summary.csv", transform_rows, list(transform_rows[0]))

    keys_a = set(motion["standalone"]["keyed"])
    keys_b = set(motion["persistent"]["keyed"])
    motion_columns = ["X_rotation(rad)", "Y_rotation(rad)", "Z_rotation(rad)",
                      "X_translation(mm)", "Y_translation(mm)", "Z_translation(mm)",
                      "Displacement(mm)", "Cumulative_displacement(mm)"]
    motion_rows = []
    for volume, group in sorted(keys_a & keys_b):
        a = motion["standalone"]["keyed"][(volume, group)]
        b = motion["persistent"]["keyed"][(volume, group)]
        assert int(a["reg_index"]) == int(b["reg_index"])
        motion_rows.append({"reg_index": int(a["reg_index"]), "volume": volume, "group": group,
                            **{f"delta_{column}": float(b[column]) - float(a[column]) for column in motion_columns}})
    write_csv("motion_trace_difference_summary.csv", motion_rows, list(motion_rows[0]))

    raw_a = {p.name: p for p in RUNS["standalone"].glob("*.raw")}
    raw_b = {p.name: p for p in RUNS["persistent"].glob("*.raw")}
    raw_matched = set(raw_a) & set(raw_b)
    raw_differing = [name for name in sorted(raw_matched) if not filecmp.cmp(raw_a[name], raw_b[name], shallow=False)]
    reference_hashes = {
        mode: sha256(next(RUNS[mode].glob("*volume_0000*_upsampled.raw")))
        for mode in RUNS
    }

    result = {
        "sources": {mode: {"queue": queue[mode]["path"], "fire": fire[mode]["path"],
                           "status": status[mode]["path"], "motion": motion[mode]["path"]} for mode in RUNS},
        "raw_inputs": {"matched_names": len(raw_matched), "unmatched_standalone": len(set(raw_a) - set(raw_b)),
                       "unmatched_persistent": len(set(raw_b) - set(raw_a)),
                       "differing_content": len(raw_differing), "first_differing": raw_differing[:10],
                       "upsampled_reference_sha256": reference_hashes},
        "runs": {},
        "comparison": {},
    }
    for mode in RUNS:
        q, f, s, t, m = queue[mode], fire[mode], status[mode], transforms[mode], motion[mode]
        records = q["records"]
        elapsed = [r["elapsed_ms"] for r in records]
        completed_at_close = sum(r["finished"] <= f["stream_close"] for r in records)
        started_at_close = sum(r["started"] <= f["stream_close"] for r in records)
        pointers = {POINTER.search(p.name).groups() for p in RUNS[mode].glob("*group_*.txt") if POINTER.search(p.name)}
        result["runs"][mode] = {
            "queue_start": str(q["started"]), "queue_finalized": str(q["finalized"]),
            "queue_final_processed": str(q["final_processed"]), "queue_reset": str(q["reset"]),
            "queue_elapsed_logged_s": q["elapsed_summary_s"],
            "queue_elapsed_timestamp_s": (q["finalized"] - q["started"]).total_seconds(),
            "stream_close": str(f["stream_close"]), "stream_duration_logged_s": f["receive_elapsed_s"],
            "drain_after_stream_close_s": (q["finalized"] - f["stream_close"]).total_seconds(),
            "fire_timeout": f["timeout"], "consolidated_at": str(f.get("consolidated_at")),
            "pointer_count": len(pointers), "opportunities": len(pointers) - 1,
            "attempts": len(records), "attempt_indices_unique": len({r["index"] for r in records}),
            "timed_registrations": sum("elapsed_ms" in r for r in records),
            "attempted_by_stream_close": started_at_close,
            "completed_by_stream_close": completed_at_close,
            "remaining_at_stream_close": len(records) - completed_at_close,
            "queue_timing_ms": stats(elapsed),
            "slimm_total_ms": stats([r["timing"]["total_ms"] for r in records if "total_ms" in r["timing"]]),
            "slimm_optimizer_ms": stats([r["timing"]["optimizer_ms"] for r in records if "optimizer_ms" in r["timing"]] if mode == "standalone" else [float(r["protocol"]["optimizer_ms"]) for r in records]),
            "target_read_ms": stats([r["timing"]["target_read_ms"] for r in records if "target_read_ms" in r["timing"]] if mode == "standalone" else [float(r["protocol"]["target_read_ms"]) for r in records]),
            "reference_setup_ms": stats([r["ref_setup"]["total_ms"] for r in records if "total_ms" in r["ref_setup"]]),
            "reference_allocation_ms": stats([r["ref_setup"]["moving_device_allocation_ms"] for r in records if "moving_device_allocation_ms" in r["ref_setup"]]),
            "reference_allocation_resized_count": sum(r["ref_setup"].get("moving_device_allocation_resized") == 1 for r in records),
            "persistent_register_request_ms": stats([float(r["protocol"]["request_ms"]) for r in records if "protocol" in r]),
            "persistent_loads": len(q["loads"]), "persistent_load_details": q["loads"],
            "first_registration_ms": elapsed[0],
            "next_100_registration_ms": stats(elapsed[1:101]),
            "steady_after_first_100_ms": stats(elapsed[100:]),
            "reference_events": q["reference_events"], "queue_failure_log_lines": q["failures"][:10],
            "status_schema_version": s["schema_version"], "status_revision": s["revision"],
            "status_counts_including_identity": s["statuses"], "finalized_volumes": s["finalized_volumes"],
            "status_entry_count": len(s["entries"]),
            "transform_count_including_identity": len(t),
            "transform_types": dict(Counter(x["type"] for x in t.values())),
            "motion_rows": len(m["rows"]), "motion_columns": m["columns"],
            "fd_mm": stats([float(row["Displacement(mm)"]) for row in m["rows"]]),
            "fd_above_0_3_count": sum(float(row["Displacement(mm)"]) > 0.3 for row in m["rows"]),
            "motion_flag_count": sum(int(row["Motion_flag"]) for row in m["rows"]),
        }

    qa = result["runs"]["standalone"]
    qb = result["runs"]["persistent"]
    result["comparison"]["queue_timing"] = {
        "mean_speedup": qa["queue_timing_ms"]["mean"] / qb["queue_timing_ms"]["mean"],
        "mean_reduction_ms": qa["queue_timing_ms"]["mean"] - qb["queue_timing_ms"]["mean"],
        "mean_reduction_pct": 100 * (1 - qb["queue_timing_ms"]["mean"] / qa["queue_timing_ms"]["mean"]),
        "total_queue_saved_s": qa["queue_elapsed_logged_s"] - qb["queue_elapsed_logged_s"],
        "total_queue_speedup": qa["queue_elapsed_logged_s"] / qb["queue_elapsed_logged_s"],
        "total_queue_reduction_pct": 100 * (1 - qb["queue_elapsed_logged_s"] / qa["queue_elapsed_logged_s"]),
    }
    result["comparison"]["transforms"] = {
        "matched_including_identity": len(matched),
        "matched_cuda_registrations": len(transform_rows),
        "unmatched_standalone": sorted(names_a - names_b), "unmatched_persistent": sorted(names_b - names_a),
        "fixed_parameter_mismatches": sum(transforms["standalone"][name]["fixed_parameters"] != transforms["persistent"][name]["fixed_parameters"] for name in matched),
        "rotation_rad": differences([row[f"d_{axis}"] for row in transform_rows for axis in ("rx", "ry", "rz")]),
        "translation_mm": differences([row[f"d_{axis}"] for row in transform_rows for axis in ("tx", "ty", "tz")]),
        "per_axis": {axis: differences([row[f"d_{axis}"] for row in transform_rows]) for axis in ("rx", "ry", "rz", "tx", "ty", "tz")},
        "largest_rotation": max(transform_rows, key=lambda row: row["max_abs_rotation_rad"]),
        "largest_translation": max(transform_rows, key=lambda row: row["max_abs_translation_mm"]),
        "over_0_001_rad": sum(row["max_abs_rotation_rad"] > 0.001 for row in transform_rows),
        "over_0_1_mm": sum(row["max_abs_translation_mm"] > 0.1 for row in transform_rows),
        "first_over_0_001_rad": next(row["transform"] for row in transform_rows if row["max_abs_rotation_rad"] > 0.001),
        "first_over_0_1_mm": next(row["transform"] for row in transform_rows if row["max_abs_translation_mm"] > 0.1),
    }
    result["comparison"]["motion"] = {
        "matched_rows": len(motion_rows), "unmatched_standalone": len(keys_a - keys_b),
        "unmatched_persistent": len(keys_b - keys_a),
        "differences": {column: differences([row[f"delta_{column}"] for row in motion_rows]) for column in motion_columns},
        "largest_fd_difference": max(motion_rows, key=lambda row: abs(row["delta_Displacement(mm)"])),
        "motion_flag_disagreements": sum(motion["standalone"]["keyed"][key]["Motion_flag"] != motion["persistent"]["keyed"][key]["Motion_flag"] for key in keys_a & keys_b),
        "fd_pearson_correlation": float(np.corrcoef(
            [float(motion["standalone"]["keyed"][key]["Displacement(mm)"]) for key in sorted(keys_a & keys_b)],
            [float(motion["persistent"]["keyed"][key]["Displacement(mm)"]) for key in sorted(keys_a & keys_b)]
        )[0, 1]),
        "fd_difference_over_0_3_mm": sum(abs(row["delta_Displacement(mm)"]) > 0.3 for row in motion_rows),
    }
    result["comparison"]["status"] = {
        "same_keys": set(status["standalone"]["entries"]) == set(status["persistent"]["entries"]),
        "different_entries": [f"{v}:{g}" for (v, g) in sorted(set(status["standalone"]["entries"]) & set(status["persistent"]["entries"])) if status["standalone"]["entries"][(v, g)] != status["persistent"]["entries"][(v, g)]],
        "missing_transform_references": {mode: [f"{v}:{g}" for (v, g), entry in status[mode]["entries"].items() if entry["status"] == "registered" and entry["transform_filename"] not in transforms[mode]] for mode in RUNS},
    }

    (OUT / "comparison_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"runs": result["runs"], "comparison": result["comparison"], "raw_inputs": result["raw_inputs"]}, indent=2))

    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(11, 4))
        x = [row["registration_index"] for row in timings]
        ax.plot(x, [row["standalone_ms"] for row in timings], lw=0.65, alpha=0.75, label="standalone")
        ax.plot(x, [row["persistent_ms"] for row in timings], lw=0.65, alpha=0.75, label="persistent")
        ax.set(xlabel="Registration index", ylabel="Queue registration time (ms)", title="FIFO CUDA A/B: per-registration wall time")
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUT / "registration_timing_vs_index.png", dpi=160)
        plt.close(fig)
    except ImportError:
        pass


if __name__ == "__main__":
    main()

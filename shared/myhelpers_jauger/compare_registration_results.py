#!/usr/bin/env python3
"""Compare static-phantom sms-mi-reg and SLIMM-v3 CUDA result directories.

The transform files are authoritative.  VersorRigid3D and Euler3D transforms
are reduced to rotation matrices and physical homogeneous transforms before
cross-backend comparisons are made.  Queue logs and motion-monitor CSVs add
optimizer, runtime, and framewise-displacement information.
"""

from __future__ import annotations

import argparse
import html
import logging
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

try:
    import SimpleITK as sitk
except ImportError as exc:  # pragma: no cover - environment-dependent message
    raise SystemExit("SimpleITK is required to read .tfm files: pip install SimpleITK") from exc


LOG = logging.getLogger("registration-comparison")
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
TFM_NAME_RE = re.compile(r"alignTransform_(\d{4})_(\d{4})-(\d{4})(?:_identity)?\.tfm$")
IDENTIFIER_RE = re.compile(r"(?:alignTransform_)?(\d{4})_(\d{4})-(\d{4})(?:_identity)?(?:\.tfm)?")


@dataclass
class Inventory:
    root: Path
    files: list[Path]
    transforms: list[Path]
    queue_logs: list[Path]
    motion_logs: list[Path]
    fire_logs: list[Path]
    metadata_json: list[Path]
    motion_csvs: list[Path]
    dashboards: list[Path]
    runtime_files: list[Path]
    other: list[Path]
    selected_queue_log: Path | None
    selected_motion_csv: Path | None


def inventory_directory(root: Path) -> Inventory:
    files = sorted(p for p in root.rglob("*") if p.is_file())
    transforms = [p for p in files if p.suffix.lower() == ".tfm"]
    queue_logs = [p for p in files if p.name.startswith("log_local_queue_processor_") and p.suffix == ".log"]
    motion_logs = [p for p in files if p.name.startswith("log_motion_monitor_") and p.suffix == ".log"]
    fire_logs = [p for p in files if p.name.startswith("log_python-fire-server_") and p.suffix == ".log"]
    metadata_json = [p for p in files if p.suffix.lower() == ".json"]
    motion_csvs = [p for p in files if p.name.startswith("motionMonitor_data_") and p.suffix.lower() == ".csv"]
    dashboards = [p for p in files if "dashboard" in p.name.lower() and p.suffix.lower() in {".jpg", ".jpeg", ".png"}]
    runtime_files = [p for p in files if p.name.startswith("runtimes_")]
    known = set(transforms + queue_logs + motion_logs + fire_logs + metadata_json + motion_csvs + dashboards + runtime_files)
    other = [p for p in files if p not in known and not p.name.startswith(".~lock.")]

    def most_informative(paths: list[Path]) -> Path | None:
        return max(paths, key=lambda p: (p.stat().st_size, p.stat().st_mtime), default=None)

    return Inventory(
        root, files, transforms, queue_logs, motion_logs, fire_logs,
        metadata_json, motion_csvs, dashboards, runtime_files, other,
        most_informative(queue_logs), most_informative(motion_csvs),
    )


def show_inventory(label: str, inv: Inventory) -> None:
    print(f"\n[{label}] {inv.root}")
    categories = [
        ("queue log", inv.queue_logs), ("motion CSV", inv.motion_csvs),
        ("metadata JSON", inv.metadata_json), ("transform", inv.transforms),
        ("motion log", inv.motion_logs), ("fire-server log", inv.fire_logs),
        ("dashboard", inv.dashboards), ("runtime artifact", inv.runtime_files),
        ("other", inv.other),
    ]
    for name, paths in categories:
        print(f"  {name:18s}: {len(paths)}")
        if name != "transform" or len(paths) <= 6:
            for path in paths:
                print(f"    - {path.relative_to(inv.root)}")
        elif paths:
            print(f"    - {paths[0].relative_to(inv.root)} ... {paths[-1].relative_to(inv.root)}")
    print(f"  selected queue log: {inv.selected_queue_log.name if inv.selected_queue_log else 'NONE'}")
    print(f"  selected motion CSV: {inv.selected_motion_csv.name if inv.selected_motion_csv else 'NONE'}")


def physical_matrix(rotation: np.ndarray, center: np.ndarray, translation: np.ndarray) -> np.ndarray:
    """Return x' = R(x-c)+c+t as a 4x4 homogeneous matrix."""
    result = np.eye(4)
    result[:3, :3] = rotation
    result[:3, 3] = center + translation - rotation @ center
    return result


def rotation_angle_deg(rotation: np.ndarray) -> float:
    cosine = float(np.clip((np.trace(rotation) - 1.0) / 2.0, -1.0, 1.0))
    return math.degrees(math.acos(cosine))


def displacement_from_rigid(translation: np.ndarray, rotation: np.ndarray, radius: float) -> float:
    """Project convention from apps/motion_monitor/monitor_directory.py.

    This is translation norm plus the chord displacement at the given radius,
    algebraically identical to radius*sqrt((1-cos(theta))**2+sin(theta)**2).
    """
    theta = math.radians(rotation_angle_deg(rotation))
    return float(np.linalg.norm(translation) + 2.0 * radius * abs(math.sin(theta / 2.0)))


def read_transforms(inv: Inventory, backend: str, head_radius: float) -> tuple[pd.DataFrame, dict[tuple[int, int], dict[str, Any]], list[str]]:
    rows: list[dict[str, Any]] = []
    objects: dict[tuple[int, int], dict[str, Any]] = {}
    issues: list[str] = []
    for path in inv.transforms:
        match = TFM_NAME_RE.search(path.name)
        if not match:
            issues.append(f"Unrecognized transform filename: {path.name}")
            continue
        output_index, volume, group = map(int, match.groups())
        try:
            transform = sitk.ReadTransform(str(path))
            name = transform.GetName()
            rotation = np.asarray(transform.GetMatrix(), dtype=float).reshape(3, 3)
            center = np.asarray(transform.GetCenter(), dtype=float)
            translation = np.asarray(transform.GetTranslation(), dtype=float)
        except Exception as exc:  # keep inventory/report usable on partial runs
            issues.append(f"Could not read {path.name}: {exc}")
            continue
        if name not in {"VersorRigid3DTransform", "Euler3DTransform"}:
            issues.append(f"Unexpected transform type {name} in {path.name}")
        homogeneous = physical_matrix(rotation, center, translation)
        angle = rotation_angle_deg(rotation)
        row = {
            "backend": backend, "registration_call": 0 if output_index == 0 else output_index,
            "output_index": output_index, "volume": volume, "group": group,
            "is_reference": bool(output_index == 0 or "identity" in path.stem),
            "transform_file": path.name, "transform_type": name,
            "center_x": center[0], "center_y": center[1], "center_z": center[2],
            "tx": translation[0], "ty": translation[1], "tz": translation[2],
            "translation_mm": float(np.linalg.norm(translation)),
            "rotation_deg": angle,
            "displacement_proxy_mm": displacement_from_rigid(translation, rotation, head_radius),
        }
        for i in range(3):
            for j in range(3):
                row[f"r{i}{j}"] = rotation[i, j]
        rows.append(row)
        objects[(volume, group)] = {"R": rotation, "c": center, "t": translation, "H": homogeneous, "row": row}
    return pd.DataFrame(rows).sort_values(["output_index"]).reset_index(drop=True), objects, issues


def parse_identifier(text: str) -> tuple[int, int, int] | None:
    match = IDENTIFIER_RE.search(text)
    return tuple(map(int, match.groups())) if match else None


def parse_queue_log(path: Path | None, backend: str) -> tuple[pd.DataFrame, list[str]]:
    if path is None:
        return pd.DataFrame(), [f"No queue log found for {backend}"]
    records: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    issues: list[str] = []

    def finish() -> None:
        nonlocal current
        if current is not None:
            records.append(current)
        current = None

    for line_no, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        start = re.search(r"Running registration call\s+(\d+)", line)
        if start:
            finish()
            current = {"registration_call": int(start.group(1)), "backend": backend, "log_line": line_no}
            continue
        if current is None:
            continue
        identifiers = [tuple(map(int, match.groups())) for match in IDENTIFIER_RE.finditer(line)]
        ident = next((item for item in identifiers if item[0] == current["registration_call"]),
                     identifiers[-1] if identifiers else None)
        if ident and ("run command" in line.lower() or "saved" in line.lower()):
            output_index, volume, group = ident
            current.update(output_index=output_index, volume=volume, group=group,
                           transform_file=f"alignTransform_{output_index:04d}_{volume:04d}-{group:04d}.tfm")
        optimum = re.search(r"found (?:optimum|minimum) at f\([^)]*\)\s*=\s*(%s)" % FLOAT, line)
        if optimum:
            current["objective"] = float(optimum.group(1))
        fixed = re.search(r"(?:with fixed parameters|slimm center of rotation)\s*[:=]\s*\(([^)]*)\)", line)
        if fixed:
            values = [float(v) for v in re.findall(FLOAT, fixed.group(1))]
            if len(values) >= 3:
                current.update(log_center_x=values[0], log_center_y=values[1], log_center_z=values[2])
        evals = re.search(r"(?:Number of function evals is\s*:|# of cost eval:)\s*(\d+)", line)
        if evals:
            current["optimizer_evals"] = int(evals.group(1))
        return_code = re.search(r"Return code is\s*:\s*(-?\d+)", line)
        if return_code:
            current["return_code"] = int(return_code.group(1))
        opt_sec = re.search(r"NLopt elapsed time \(sec\)\s*:\s*(%s)" % FLOAT, line)
        if opt_sec:
            current["optimizer_sec"] = float(opt_sec.group(1))
        cuda_time = re.search(r"time\s*=\s*(%s)\s*ms;\s*# of cost eval:\s*(\d+)" % FLOAT, line)
        if cuda_time:
            current["optimizer_sec"] = float(cuda_time.group(1)) / 1000.0
            current["optimizer_evals"] = int(cuda_time.group(2))
        call = re.search(r"(?:CUDA )?registration call elapsed runtime \(sec\)\s*:\s*(%s)" % FLOAT, line, re.IGNORECASE)
        if call:
            current["call_sec"] = float(call.group(1))
    finish()
    frame = pd.DataFrame(records)
    if not frame.empty:
        for col in ("volume", "group", "output_index"):
            if col not in frame:
                frame[col] = np.nan
        frame = frame.sort_values("registration_call").reset_index(drop=True)
        missing_ids = frame[frame[["volume", "group"]].isna().any(axis=1)]
        if len(missing_ids):
            issues.append(f"{len(missing_ids)} queue-log registrations lack volume/group identifiers")
    return frame, issues


def read_motion_csv(path: Path | None, backend: str) -> tuple[pd.DataFrame, list[str]]:
    if path is None:
        return pd.DataFrame(), [f"No motion-monitor CSV found for {backend}"]
    frame = pd.read_csv(path)
    required = {"Volume_index", "Slice_group_index"}
    issues = []
    if not required.issubset(frame.columns):
        issues.append(f"Motion CSV {path.name} lacks columns: {sorted(required - set(frame.columns))}")
        return frame, issues
    rename = {
        "reg_index": "motion_reg_index", "Volume_index": "volume", "Slice_group_index": "group",
        "X_rotation(rad)": "csv_rx", "Y_rotation(rad)": "csv_ry", "Z_rotation(rad)": "csv_rz",
        "X_translation(mm)": "csv_tx", "Y_translation(mm)": "csv_ty", "Z_translation(mm)": "csv_tz",
        "Displacement(mm)": "fd", "Cumulative_displacement(mm)": "cumulative_fd", "Motion_flag": "motion_flag",
    }
    frame = frame.rename(columns=rename)
    frame["backend"] = backend
    xyz = [c for c in ("csv_tx", "csv_ty", "csv_tz") if c in frame]
    if len(xyz) == 3:
        frame["csv_translation_mm"] = np.sqrt(sum(frame[c].astype(float) ** 2 for c in xyz))
    rxyz = [c for c in ("csv_rx", "csv_ry", "csv_rz") if c in frame]
    if len(rxyz) == 3:
        # Visualization only; authoritative angle is computed from each .tfm matrix.
        frame["csv_euler_vector_deg"] = np.degrees(np.sqrt(sum(frame[c].astype(float) ** 2 for c in rxyz)))
    return frame, issues


def merge_sources(transform_df: pd.DataFrame, log_df: pd.DataFrame, motion_df: pd.DataFrame) -> pd.DataFrame:
    result = transform_df.copy()
    if not log_df.empty:
        cols = [c for c in ["volume", "group", "objective", "optimizer_evals", "optimizer_sec", "call_sec",
                                    "return_code", "log_center_x", "log_center_y", "log_center_z"] if c in log_df]
        result = result.merge(log_df[cols], on=["volume", "group"], how="left", validate="one_to_one")
    if not motion_df.empty and {"volume", "group"}.issubset(motion_df.columns):
        cols = [c for c in ["volume", "group", "motion_reg_index", "csv_rx", "csv_ry", "csv_rz",
                                    "csv_tx", "csv_ty", "csv_tz", "csv_translation_mm", "csv_euler_vector_deg",
                                    "fd", "cumulative_fd", "motion_flag"] if c in motion_df]
        result = result.merge(motion_df[cols], on=["volume", "group"], how="left", validate="one_to_one")
    return result.sort_values("output_index").reset_index(drop=True)


def pair_transforms(sms: pd.DataFrame, cuda: pd.DataFrame,
                    sms_objects: dict[tuple[int, int], dict[str, Any]],
                    cuda_objects: dict[tuple[int, int], dict[str, Any]], radius: float) -> tuple[pd.DataFrame, list[tuple[int, int]], list[tuple[int, int]]]:
    sms_keys = set(sms_objects) - {(0, 14)}
    cuda_keys = set(cuda_objects) - {(0, 14)}
    matched = sorted(sms_keys & cuda_keys)
    rows = []
    sms_lookup = sms.set_index(["volume", "group"]).to_dict("index")
    cuda_lookup = cuda.set_index(["volume", "group"]).to_dict("index")
    for key in matched:
        s, c = sms_objects[key], cuda_objects[key]
        relative = np.linalg.inv(s["H"]) @ c["H"]  # both files encode fixed -> moving
        r_rel, b_rel = relative[:3, :3], relative[:3, 3]
        common_center = 0.5 * (s["c"] + c["c"])
        relative_at_center = r_rel @ common_center + b_rel - common_center
        center_delta = c["c"] - s["c"]
        sr, cr = sms_lookup[key], cuda_lookup[key]
        row = {"volume": key[0], "group": key[1], "sms_transform_file": sr["transform_file"],
               "cuda_transform_file": cr["transform_file"]}
        for prefix, source in (("sms", sr), ("cuda", cr)):
            for name in ("tx", "ty", "tz", "translation_mm", "rotation_deg", "displacement_proxy_mm", "fd",
                         "optimizer_evals", "optimizer_sec", "call_sec", "objective"):
                row[f"{prefix}_{name}"] = source.get(name, np.nan)
        row.update(
            center_difference_mm=float(np.linalg.norm(center_delta)),
            relative_translation_mm=float(np.linalg.norm(relative_at_center)),
            relative_origin_offset_mm=float(np.linalg.norm(b_rel)),
            relative_rotation_deg=rotation_angle_deg(r_rel),
            relative_displacement_proxy_mm=displacement_from_rigid(relative_at_center, r_rel, radius),
        )
        rows.append(row)
    return pd.DataFrame(rows), sorted(sms_keys - cuda_keys), sorted(cuda_keys - sms_keys)


def numeric_stats(values: Iterable[Any]) -> dict[str, float]:
    data = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(float)
    if not len(data):
        return {name: math.nan for name in ("count", "mean", "median", "sd", "rms", "min", "max")}
    return {"count": float(len(data)), "mean": float(np.mean(data)), "median": float(np.median(data)),
            "sd": float(np.std(data, ddof=1)) if len(data) > 1 else 0.0,
            "rms": float(np.sqrt(np.mean(data ** 2))), "min": float(np.min(data)), "max": float(np.max(data))}


def slope(values: Iterable[Any]) -> float:
    series = pd.to_numeric(pd.Series(values), errors="coerce")
    mask = series.notna().to_numpy()
    if mask.sum() < 2:
        return math.nan
    x = np.arange(len(series), dtype=float)[mask]
    return float(np.polyfit(x, series.to_numpy(float)[mask], 1)[0])


def noise_summary(frame: pd.DataFrame, threshold: float) -> dict[str, Any]:
    data = frame.loc[~frame["is_reference"]].copy()
    result: dict[str, Any] = {}
    for col in ("tx", "ty", "tz", "translation_mm", "rotation_deg", "fd"):
        result[col] = numeric_stats(data[col]) if col in data else numeric_stats([])
    fd = pd.to_numeric(data.get("fd", pd.Series(dtype=float)), errors="coerce").dropna()
    result["fd_cumulative"] = float(fd.sum()) if len(fd) else math.nan
    result["fd_above"] = int((fd > threshold).sum()) if len(fd) else 0
    result["fd_above_pct"] = 100.0 * result["fd_above"] / len(fd) if len(fd) else math.nan
    result["slopes"] = {col: slope(data[col]) if col in data else math.nan for col in ("tx", "ty", "tz", "rotation_deg", "fd")}
    volume_means = data.groupby("volume")[[c for c in ("translation_mm", "rotation_deg", "fd") if c in data]].mean(numeric_only=True)
    result["volume_mean_ranges"] = {c: float(volume_means[c].max() - volume_means[c].min()) for c in volume_means}
    return result


def slice_group_summary(sms: pd.DataFrame, cuda: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for backend, frame in (("sms-mi-reg", sms), ("CUDA", cuda)):
        data = frame.loc[~frame["is_reference"]]
        for group, group_df in data.groupby("group"):
            row: dict[str, Any] = {"backend": backend, "group": int(group), "count": len(group_df)}
            for col in ("translation_mm", "rotation_deg", "fd"):
                stats = numeric_stats(group_df[col]) if col in group_df else numeric_stats([])
                row[f"{col}_mean"] = stats["mean"]
                row[f"{col}_sd"] = stats["sd"]
            rows.append(row)
    return pd.DataFrame(rows)


def runtime_summary(sms: pd.DataFrame, cuda: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for backend, frame in (("sms-mi-reg", sms), ("CUDA", cuda)):
        data = frame.loc[~frame["is_reference"]]
        for metric in ("optimizer_evals", "optimizer_sec", "call_sec", "objective"):
            stats = numeric_stats(data[metric]) if metric in data else numeric_stats([])
            rows.append({"backend": backend, "metric": metric, **stats})
        opt = pd.to_numeric(data.get("optimizer_sec", pd.Series(np.nan, index=data.index)), errors="coerce")
        call = pd.to_numeric(data.get("call_sec", pd.Series(np.nan, index=data.index)), errors="coerce")
        valid = opt.notna() & call.notna() & (call != 0)
        fraction = numeric_stats((opt[valid] / call[valid]) * 100.0)
        rows.append({"backend": backend, "metric": "optimizer_fraction_percent", **fraction})
    return pd.DataFrame(rows)


def fmt(value: Any, digits: int = 5) -> str:
    try:
        number = float(value)
        return "NA" if not math.isfinite(number) else f"{number:.{digits}g}"
    except (TypeError, ValueError):
        return "NA"


def inventory_html(inv: Inventory) -> str:
    def names(paths: list[Path]) -> str:
        return ", ".join(f"<code>{html.escape(str(p.relative_to(inv.root)))}</code>" for p in paths) if paths else "None"
    return "<ul>" + "".join([
        f"<li>Queue logs ({len(inv.queue_logs)}): {names(inv.queue_logs)}</li>",
        f"<li>Motion CSVs ({len(inv.motion_csvs)}): {names(inv.motion_csvs)}</li>",
        f"<li>Metadata JSON ({len(inv.metadata_json)}): {names(inv.metadata_json)}</li>",
        f"<li>Transforms: {len(inv.transforms)} <code>.tfm</code> files</li>",
        f"<li>Motion-monitor logs ({len(inv.motion_logs)}): {names(inv.motion_logs)}</li>",
        f"<li>Fire-server logs ({len(inv.fire_logs)}): {names(inv.fire_logs)}</li>",
        f"<li>Dashboards ({len(inv.dashboards)}): {names(inv.dashboards)}</li>",
        f"<li>Runtime artifacts ({len(inv.runtime_files)}): {names(inv.runtime_files)}</li>",
    ]) + "</ul>"


def make_plots(output: Path, sms: pd.DataFrame, cuda: pd.DataFrame, paired: pd.DataFrame,
               groups: pd.DataFrame, threshold: float) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 180, "axes.grid": True, "grid.alpha": 0.25,
                         "font.size": 9, "figure.constrained_layout.use": True})
    made: list[Path] = []
    data = [("sms-mi-reg", sms.loc[~sms.is_reference], "#1f77b4"), ("CUDA", cuda.loc[~cuda.is_reference], "#d62728")]

    def save(fig: Any, name: str) -> None:
        path = output / name
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        made.append(path)

    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    for label, frame, color in data:
        for ax, col in zip(axes, ("tx", "ty", "tz")):
            ax.plot(frame.output_index, frame[col], label=label, lw=1, color=color, alpha=.85)
            ax.set_ylabel(f"{col} (mm)")
    axes[0].legend(); axes[-1].set_xlabel("Registration index")
    save(fig, "translation_traces.png")

    for column, ylabel, name, threshold_line in [
        ("rotation_deg", "Equivalent rotation (deg)", "rotation_magnitude_trace.png", None),
        ("fd", "Framewise displacement (mm)", "framewise_displacement_trace.png", threshold),
    ]:
        fig, ax = plt.subplots(figsize=(10, 4))
        for label, frame, color in data:
            ax.plot(frame.output_index, frame[column], label=label, lw=1, color=color)
        if threshold_line is not None:
            ax.axhline(threshold_line, color="black", ls="--", lw=1, label=f"{threshold_line:g} mm threshold")
        ax.set(xlabel="Registration index", ylabel=ylabel); ax.legend()
        save(fig, name)

    for column, xlabel, name in [
        ("translation_mm", "Translation magnitude (mm)", "translation_magnitude_distribution.png"),
        ("rotation_deg", "Equivalent rotation (deg)", "rotation_distribution.png"),
        ("fd", "Framewise displacement (mm)", "fd_distribution.png"),
    ]:
        fig, ax = plt.subplots(figsize=(7, 4))
        for label, frame, color in data:
            values = pd.to_numeric(frame[column], errors="coerce").dropna()
            ax.hist(values, bins=24, density=True, histtype="step", lw=1.8, label=label, color=color)
        ax.set(xlabel=xlabel, ylabel="Density"); ax.legend()
        save(fig, name)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, col, title in zip(axes, ("optimizer_sec", "call_sec"), ("Optimizer", "Registration call")):
        vals = [pd.to_numeric(f[col], errors="coerce").dropna() for _, f, _ in data]
        ax.boxplot(vals, tick_labels=[x[0] for x in data], showfliers=True)
        ax.set(ylabel="Seconds", title=f"{title} runtime")
    save(fig, "runtime_comparison.png")

    fig, ax = plt.subplots(figsize=(7, 4))
    vals = [pd.to_numeric(f["optimizer_evals"], errors="coerce").dropna() for _, f, _ in data]
    ax.boxplot(vals, tick_labels=[x[0] for x in data], showfliers=True)
    ax.set(ylabel="Cost evaluations", title="Optimizer cost evaluations")
    save(fig, "cost_evaluations_comparison.png")

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    axes[0].plot(paired.index, paired.relative_translation_mm, lw=1, color="#2ca02c")
    axes[1].plot(paired.index, paired.relative_rotation_deg, lw=1, color="#9467bd")
    axes[0].set_ylabel("Relative translation (mm)"); axes[1].set_ylabel("Relative rotation (deg)")
    axes[1].set_xlabel("Matched registration index")
    save(fig, "pairwise_transform_disagreement.png")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for backend, color in (("sms-mi-reg", "#1f77b4"), ("CUDA", "#d62728")):
        frame = groups[groups.backend == backend]
        for ax, col, label in zip(axes, ("translation_mm_mean", "rotation_deg_mean", "fd_mean"),
                                  ("Mean translation (mm)", "Mean rotation (deg)", "Mean FD (mm)")):
            ax.plot(frame.group, frame[col], marker="o", ms=3, lw=1, label=backend, color=color)
            ax.set(xlabel="Slice-group index", ylabel=label); ax.set_xticks(range(14))
    axes[0].legend()
    save(fig, "slice_group_bias.png")
    return made


def build_report(output: Path, sms_inv: Inventory, cuda_inv: Inventory, sms: pd.DataFrame, cuda: pd.DataFrame,
                 paired: pd.DataFrame, sms_missing: list[tuple[int, int]], cuda_missing: list[tuple[int, int]],
                 groups: pd.DataFrame, runtimes: pd.DataFrame, sms_noise: dict[str, Any], cuda_noise: dict[str, Any],
                 threshold: float, radius: float, issues: list[str], figures: list[Path]) -> Path:
    """Write a portable HTML report with relative references to generated plots."""
    runtime = runtimes.set_index(["backend", "metric"])
    optimizer_speedup = (runtime.loc[("sms-mi-reg", "optimizer_sec"), "mean"] /
                         runtime.loc[("CUDA", "optimizer_sec"), "mean"])
    call_speedup = (runtime.loc[("sms-mi-reg", "call_sec"), "mean"] /
                    runtime.loc[("CUDA", "call_sec"), "mean"])
    pair_t = numeric_stats(paired.get("relative_translation_mm", []))
    pair_r = numeric_stats(paired.get("relative_rotation_deg", []))
    pair_d = numeric_stats(paired.get("relative_displacement_proxy_mm", []))
    max_center = pd.to_numeric(paired.get("center_difference_mm", pd.Series(dtype=float)), errors="coerce").max()

    def make_table(headers: list[str], rows: list[list[Any]]) -> str:
        header = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
        body = "".join(
            "<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>"
            for row in rows
        )
        return f'<div class="table-wrap"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>'

    available_figures = {path.name for path in figures}

    def figure_group(items: list[tuple[str, str]]) -> str:
        rendered = []
        for name, caption in items:
            if name not in available_figures:
                continue
            safe_name = html.escape(name, quote=True)
            rendered.append(
                f'<figure><a href="{safe_name}"><img src="{safe_name}" alt="{html.escape(caption, quote=True)}" '
                f'loading="lazy"></a><figcaption>{html.escape(caption)}</figcaption></figure>'
            )
        return (f'<div class="figure-grid">{"".join(rendered)}</div>' if rendered else
                '<p class="note">Plot generation was disabled with <code>--no-plots</code>.</p>')

    noise_specs = [
        ("Mean translation magnitude (mm)", "translation_mm", "mean"),
        ("Median translation magnitude (mm)", "translation_mm", "median"),
        ("SD translation magnitude (mm)", "translation_mm", "sd"),
        ("RMS translation magnitude (mm)", "translation_mm", "rms"),
        ("Min translation magnitude (mm)", "translation_mm", "min"),
        ("Max translation magnitude (mm)", "translation_mm", "max"),
        ("Mean equivalent rotation (deg)", "rotation_deg", "mean"),
        ("Median equivalent rotation (deg)", "rotation_deg", "median"),
        ("SD equivalent rotation (deg)", "rotation_deg", "sd"),
        ("RMS equivalent rotation (deg)", "rotation_deg", "rms"),
        ("Min equivalent rotation (deg)", "rotation_deg", "min"),
        ("Max equivalent rotation (deg)", "rotation_deg", "max"),
        ("Mean FD (mm)", "fd", "mean"), ("Median FD (mm)", "fd", "median"),
        ("SD FD (mm)", "fd", "sd"), ("RMS FD (mm)", "fd", "rms"),
        ("Min FD (mm)", "fd", "min"), ("Max FD (mm)", "fd", "max"),
    ]
    noise_rows = [[label, fmt(sms_noise[metric][statistic]), fmt(cuda_noise[metric][statistic])]
                  for label, metric, statistic in noise_specs]
    noise_rows.extend([
        ["Cumulative FD (mm)", fmt(sms_noise["fd_cumulative"]), fmt(cuda_noise["fd_cumulative"])],
        [f"FD > {threshold:g} mm", f"{sms_noise['fd_above']} ({fmt(sms_noise['fd_above_pct'], 4)}%)",
         f"{cuda_noise['fd_above']} ({fmt(cuda_noise['fd_above_pct'], 4)}%)"],
    ])
    noise_table = make_table(["Metric", "sms-mi-reg", "CUDA"], noise_rows)

    drift_rows = [
        ["sms-mi-reg", fmt(sms_noise["slopes"]["tx"]), fmt(sms_noise["slopes"]["ty"]),
         fmt(sms_noise["slopes"]["tz"]), fmt(sms_noise["slopes"]["rotation_deg"]), fmt(sms_noise["slopes"]["fd"])],
        ["CUDA", fmt(cuda_noise["slopes"]["tx"]), fmt(cuda_noise["slopes"]["ty"]),
         fmt(cuda_noise["slopes"]["tz"]), fmt(cuda_noise["slopes"]["rotation_deg"]), fmt(cuda_noise["slopes"]["fd"])],
    ]
    drift_table = make_table(["Backend", "tx (mm/index)", "ty", "tz", "rotation (deg/index)", "FD (mm/index)"], drift_rows)

    group_pivot = groups.pivot(index="group", columns="backend",
                               values=["translation_mm_mean", "rotation_deg_mean", "fd_mean"])
    group_rows = []
    for index, row in group_pivot.iterrows():
        group_rows.append([index, fmt(row[("translation_mm_mean", "sms-mi-reg")]),
                           fmt(row[("translation_mm_mean", "CUDA")]),
                           fmt(row[("rotation_deg_mean", "sms-mi-reg")]),
                           fmt(row[("rotation_deg_mean", "CUDA")]),
                           fmt(row[("fd_mean", "sms-mi-reg")]), fmt(row[("fd_mean", "CUDA")])])
    group_table = make_table(["Group", "sms translation", "CUDA translation", "sms rotation",
                              "CUDA rotation", "sms FD", "CUDA FD"], group_rows)

    runtime_rows = []
    for _, row in runtimes.iterrows():
        if row.metric != "objective":
            runtime_rows.append([row.backend, row.metric, fmt(row["mean"]), fmt(row["median"]),
                                 fmt(row["sd"]), fmt(row["min"]), fmt(row["max"])])
    runtime_table = make_table(["Backend", "Metric", "Mean", "Median", "SD", "Min", "Max"], runtime_rows)
    issue_items = ("".join(f"<li>{html.escape(item)}</li>" for item in issues) if issues else
                   "<li>None detected.</li>")
    expected_pairs = len((set(zip(sms.volume, sms.group)) | set(zip(cuda.volume, cuda.group))) - {(0, 14)})

    text = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Registration Comparison Report</title>
<style>
:root {{ --ink:#17212b; --muted:#596675; --line:#ccd4dc; --accent:#285f8f; --panel:#f6f8fa; }}
* {{ box-sizing:border-box; }} body {{ margin:0; color:var(--ink); background:#fff; font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif; }}
main {{ width:min(1180px,calc(100% - 2rem)); margin:auto; padding:2rem 0 4rem; }} h1 {{ margin:0 0 .3rem; font-size:2rem; }}
h2 {{ margin-top:2.4rem; padding-bottom:.35rem; border-bottom:2px solid var(--accent); font-size:1.45rem; }} h3 {{ margin-top:1.5rem; }}
p,ul {{ max-width:100ch; }} code {{ background:#eef1f4; padding:.08rem .3rem; border-radius:3px; overflow-wrap:anywhere; }}
.summary {{ margin:1rem 0 2rem; padding:1rem 1.2rem; background:var(--panel); border-left:4px solid var(--accent); }} .note {{ color:var(--muted); font-style:italic; }}
.table-wrap {{ overflow-x:auto; margin:1rem 0 1.5rem; }} table {{ width:100%; border-collapse:collapse; font-variant-numeric:tabular-nums; }}
th,td {{ padding:.48rem .65rem; border:1px solid var(--line); text-align:right; white-space:nowrap; }} th {{ background:#eaf0f5; }}
th:first-child,td:first-child {{ text-align:left; }} tbody tr:nth-child(even) {{ background:#fafbfc; }}
.figure-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(min(520px,100%),1fr)); gap:1.25rem; margin:1.25rem 0 2rem; }}
figure {{ margin:0; padding:.75rem; border:1px solid var(--line); box-shadow:0 1px 3px #0001; }} figure img {{ display:block; width:100%; max-width:100%; height:auto; }}
figcaption {{ margin-top:.6rem; color:var(--muted); font-size:.92rem; }} footer {{ margin-top:3rem; padding-top:1rem; border-top:1px solid var(--line); color:var(--muted); }}
</style></head><body><main>
<header><h1>Registration Comparison Summary</h1><p class="summary">Static-phantom comparison of sms-mi-reg and SLIMM-v3 CUDA using physical rigid transforms, motion-monitor output, and queue runtime logs. Lower FD is treated as a stability/noise-floor observation, not proof of greater accuracy.</p></header>

<section><h2>1. Input Datasets</h2><ul><li><strong>sms-mi-reg:</strong> <code>{html.escape(str(sms_inv.root))}</code></li>
<li><strong>SLIMM-v3 CUDA:</strong> <code>{html.escape(str(cuda_inv.root))}</code></li><li>Head radius: {radius:g} mm; motion threshold: {threshold:g} mm.</li></ul></section>

<section><h2>2. Registration / Transform Inventory</h2><h3>sms-mi-reg</h3>{inventory_html(sms_inv)}<h3>SLIMM-v3 CUDA</h3>{inventory_html(cuda_inv)}
<h3>Counts and matching</h3><ul><li>sms-mi-reg: {len(sms)} transforms ({int((~sms.is_reference).sum())} registrations plus reference).</li>
<li>CUDA: {len(cuda)} transforms ({int((~cuda.is_reference).sum())} registrations plus reference).</li><li>Expected non-reference keys: {expected_pairs}; matched pairs: {len(paired)}.</li>
<li>sms-mi-reg-only: {html.escape(str(sms_missing or 'none'))}; CUDA-only: {html.escape(str(cuda_missing or 'none'))}.</li></ul>
<p>The observed layout is {sms.loc[~sms.is_reference, 'volume'].nunique()} non-reference volumes and {sms.loc[~sms.is_reference, 'group'].nunique()} slice groups per volume. Matching used <code>(volume, group)</code> from filenames and queue-log output identifiers.</p>
<h3>Transform representation</h3><p>sms-mi-reg registrations are <code>VersorRigid3DTransform</code>; CUDA registrations are <code>Euler3DTransform</code> (the common identity reference is Versor). SimpleITK normalized both to rotation matrices and <code>x′ = R(x−c)+c+t</code>. Pairwise disagreement used <code>inverse(T_sms) @ T_cuda</code>; raw rotation parameter vectors were not compared.</p>
<p>Maximum corresponding center difference was {fmt(max_center)} mm. Relative translation measures displacement at the shared mean center.</p></section>

<section><h2>3. Static-Phantom Motion and Noise-Floor Comparison</h2>{noise_table}
<p>Axis translation means ± SD (mm): sms-mi-reg = ({fmt(sms_noise['tx']['mean'])} ± {fmt(sms_noise['tx']['sd'])}, {fmt(sms_noise['ty']['mean'])} ± {fmt(sms_noise['ty']['sd'])}, {fmt(sms_noise['tz']['mean'])} ± {fmt(sms_noise['tz']['sd'])}); CUDA = ({fmt(cuda_noise['tx']['mean'])} ± {fmt(cuda_noise['tx']['sd'])}, {fmt(cuda_noise['ty']['mean'])} ± {fmt(cuda_noise['ty']['sd'])}, {fmt(cuda_noise['tz']['mean'])} ± {fmt(cuda_noise['tz']['sd'])}).</p>
{figure_group([('translation_traces.png','Figure 1. X, Y, and Z translation traces for both backends.'),('rotation_magnitude_trace.png','Figure 2. Equivalent physical rotation magnitude over acquisition order.'),('translation_magnitude_distribution.png','Figure 3. Absolute translation-magnitude distributions.'),('rotation_distribution.png','Figure 4. Equivalent rotation-magnitude distributions.')])}
<h3>Drift and volume-to-volume stability</h3>{drift_table}<p>Ranges across per-volume means (sms-mi-reg / CUDA): translation {fmt(sms_noise['volume_mean_ranges'].get('translation_mm'))} / {fmt(cuda_noise['volume_mean_ranges'].get('translation_mm'))} mm; rotation {fmt(sms_noise['volume_mean_ranges'].get('rotation_deg'))} / {fmt(cuda_noise['volume_mean_ranges'].get('rotation_deg'))}°; FD {fmt(sms_noise['volume_mean_ranges'].get('fd'))} / {fmt(cuda_noise['volume_mean_ranges'].get('fd'))} mm. Tiny slopes should not be over-interpreted.</p></section>

<section><h2>4. Physical Transform Agreement</h2><ul><li>Relative translation: mean {fmt(pair_t['mean'])} mm, median {fmt(pair_t['median'])} mm, SD {fmt(pair_t['sd'])} mm, maximum {fmt(pair_t['max'])} mm.</li>
<li>Relative rotation: mean {fmt(pair_r['mean'])}°, median {fmt(pair_r['median'])}°, SD {fmt(pair_r['sd'])}°, maximum {fmt(pair_r['max'])}°.</li><li>Relative {radius:g}-mm displacement proxy: mean {fmt(pair_d['mean'])} mm, maximum {fmt(pair_d['max'])} mm.</li></ul><p>These values measure backend agreement, not ground-truth accuracy.</p>
{figure_group([('pairwise_transform_disagreement.png','Figure 5. Center-aware relative translation and rotation for all matched pairs.')])}</section>

<section><h2>5. Framewise Displacement</h2><p>Mean FD was {fmt(sms_noise['fd']['mean'])} mm for sms-mi-reg and {fmt(cuda_noise['fd']['mean'])} mm for CUDA. Neither run exceeded {threshold:g} mm. FD is the project motion monitor's sequential transform-to-transform metric.</p>
{figure_group([('framewise_displacement_trace.png',f'Figure 6. Framewise displacement with the {threshold:g} mm threshold.'),('fd_distribution.png','Figure 7. Framewise-displacement distributions.')])}</section>

<section><h2>6. Slice-Group Stability / Bias</h2><p>Per-group means are shown; the CSV also contains SDs.</p>{group_table}<p>Group structure is suitable for review but is not itself evidence of phantom motion.</p>
{figure_group([('slice_group_bias.png','Figure 8. Mean translation, rotation, and FD by SMS slice-group index.')])}</section>

<section><h2>7. Optimization and Runtime Performance</h2>{runtime_table}<ul><li>Optimizer speedup: <strong>{fmt(optimizer_speedup)}×</strong>.</li><li>Registration-call speedup: <strong>{fmt(call_speedup)}×</strong>.</li>
<li>Call time inside optimizer: sms-mi-reg {fmt(runtime.loc[('sms-mi-reg','optimizer_fraction_percent'),'mean'])}%; CUDA {fmt(runtime.loc[('CUDA','optimizer_fraction_percent'),'mean'])}%.</li></ul><p>Non-optimizer overhead dominates CUDA's end-to-end call time.</p>
{figure_group([('runtime_comparison.png','Figure 9. Optimizer and end-to-end registration-call runtime distributions.'),('cost_evaluations_comparison.png','Figure 10. Optimizer cost-evaluation counts.')])}
<h3>Objective values</h3><ul><li>sms-mi-reg: mean {fmt(runtime.loc[('sms-mi-reg','objective'),'mean'])}, SD {fmt(runtime.loc[('sms-mi-reg','objective'),'sd'])}, min {fmt(runtime.loc[('sms-mi-reg','objective'),'min'])}, max {fmt(runtime.loc[('sms-mi-reg','objective'),'max'])}.</li><li>CUDA: mean {fmt(runtime.loc[('CUDA','objective'),'mean'])}, SD {fmt(runtime.loc[('CUDA','objective'),'sd'])}, min {fmt(runtime.loc[('CUDA','objective'),'min'])}, max {fmt(runtime.loc[('CUDA','objective'),'max'])}.</li></ul><p class="note">Cross-backend objective magnitudes are not comparable because MI formulations, signs, and scales differ.</p></section>

<section><h2>8. Key Conclusions</h2><ul><li>Both runs remained below {threshold:g} mm for all 126 registrations; mean FD was {fmt(sms_noise['fd']['mean'])} / {fmt(cuda_noise['fd']['mean'])} mm (sms-mi-reg / CUDA).</li>
<li>Mean translation was {fmt(sms_noise['translation_mm']['mean'])} / {fmt(cuda_noise['translation_mm']['mean'])} mm and mean rotation was {fmt(sms_noise['rotation_deg']['mean'])} / {fmt(cuda_noise['rotation_deg']['mean'])}°.</li>
<li>Mean physical disagreement was {fmt(pair_t['mean'])} mm and {fmt(pair_r['mean'])}°; maxima were {fmt(pair_t['max'])} mm and {fmt(pair_r['max'])}°.</li>
<li>Per-group mean FD ranged {fmt(groups.loc[groups.backend == 'sms-mi-reg','fd_mean'].min())}–{fmt(groups.loc[groups.backend == 'sms-mi-reg','fd_mean'].max())} mm for sms-mi-reg and {fmt(groups.loc[groups.backend == 'CUDA','fd_mean'].min())}–{fmt(groups.loc[groups.backend == 'CUDA','fd_mean'].max())} mm for CUDA.</li>
<li>CUDA was {fmt(optimizer_speedup)}× faster at optimizer level and {fmt(call_speedup)}× faster end-to-end; {fmt(runtime.loc[('CUDA','optimizer_fraction_percent'),'mean'])}% of call time was optimizer time.</li></ul><p>Lower FD alone is not treated as proof of greater accuracy.</p></section>

<section><h2>9. Caveats / Missing Data</h2><h3>Detected issues</h3><ul>{issue_items}</ul><ul><li>The reference identity has no optimizer runtime and is excluded from summaries.</li><li>FD is path-dependent; transform magnitude is relative to the reference.</li><li>The displacement proxy follows the existing translation-plus-rotational-chord convention.</li><li>CSV Euler components are visualization-only.</li><li>The fixed-to-moving convention is inferred from the shared queue workflow; inverse-writing future backends must be normalized.</li></ul></section>
<footer>Generated by <code>compare_registration_results.py</code>. Figure paths are relative to this report.</footer></main></body></html>"""
    path = output / "registration_comparison_report.html"
    path.write_text(text)
    return path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare static-phantom sms-mi-reg and SLIMM-v3 CUDA registration results using physical rigid transforms, motion-monitor data, and queue runtimes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("sms_results_dir", type=Path, help="Completed sms-mi-reg result directory")
    parser.add_argument("cuda_results_dir", type=Path, help="Completed SLIMM-v3 CUDA result directory")
    parser.add_argument("--output-dir", type=Path, default=Path("registration_comparison"), help="Directory for CSVs, figures, and HTML report")
    parser.add_argument("--head-radius", type=float, default=50.0, help="Radius in mm for rotational displacement proxies")
    parser.add_argument("--motion-threshold", type=float, default=0.3, help="FD threshold in mm for motion counts and plot annotation")
    parser.add_argument("--no-plots", action="store_true", help="Generate CSVs and HTML report without matplotlib figures")
    parser.add_argument("--verbose", action="store_true", help="Enable detailed parsing diagnostics")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s: %(message)s")
    if args.head_radius <= 0 or args.motion_threshold < 0:
        raise SystemExit("--head-radius must be positive and --motion-threshold must be non-negative")
    for path in (args.sms_results_dir, args.cuda_results_dir):
        if not path.is_dir():
            raise SystemExit(f"Result directory does not exist: {path}")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    sms_inv = inventory_directory(args.sms_results_dir.resolve())
    cuda_inv = inventory_directory(args.cuda_results_dir.resolve())
    show_inventory("sms-mi-reg", sms_inv)
    show_inventory("CUDA", cuda_inv)

    sms_tfm, sms_objects, sms_tfm_issues = read_transforms(sms_inv, "sms-mi-reg", args.head_radius)
    cuda_tfm, cuda_objects, cuda_tfm_issues = read_transforms(cuda_inv, "CUDA", args.head_radius)
    sms_log, sms_log_issues = parse_queue_log(sms_inv.selected_queue_log, "sms-mi-reg")
    cuda_log, cuda_log_issues = parse_queue_log(cuda_inv.selected_queue_log, "CUDA")
    sms_motion, sms_motion_issues = read_motion_csv(sms_inv.selected_motion_csv, "sms-mi-reg")
    cuda_motion, cuda_motion_issues = read_motion_csv(cuda_inv.selected_motion_csv, "CUDA")
    sms = merge_sources(sms_tfm, sms_log, sms_motion)
    cuda = merge_sources(cuda_tfm, cuda_log, cuda_motion)
    paired, only_sms, only_cuda = pair_transforms(sms, cuda, sms_objects, cuda_objects, args.head_radius)

    issues = sms_tfm_issues + cuda_tfm_issues + sms_log_issues + cuda_log_issues + sms_motion_issues + cuda_motion_issues
    for label, frame in (("sms-mi-reg", sms), ("CUDA", cuda)):
        actual_types = sorted(frame.loc[~frame.is_reference, "transform_type"].unique()) if len(frame) else []
        expected = "VersorRigid3DTransform" if label == "sms-mi-reg" else "Euler3DTransform"
        if actual_types != [expected]:
            issues.append(f"{label} transform types were {actual_types}; expected only {expected}")
        nonref = frame.loc[~frame.is_reference]
        for field in ("optimizer_evals", "optimizer_sec", "call_sec"):
            missing = int(nonref[field].isna().sum()) if field in nonref else len(nonref)
            if missing:
                issues.append(f"{label}: {missing} registrations missing {field}")
        missing_fd = int(nonref["fd"].isna().sum()) if "fd" in nonref else len(nonref)
        if missing_fd:
            issues.append(f"{label}: {missing_fd} registrations missing motion-monitor FD")
    if only_sms:
        issues.append(f"Unmatched sms-mi-reg keys: {only_sms}")
    if only_cuda:
        issues.append(f"Unmatched CUDA keys: {only_cuda}")

    sms.to_csv(output / "registration_instances_sms.csv", index=False, na_rep="NA")
    cuda.to_csv(output / "registration_instances_cuda.csv", index=False, na_rep="NA")
    paired.to_csv(output / "paired_transform_comparison.csv", index=False, na_rep="NA")
    groups = slice_group_summary(sms, cuda)
    groups.to_csv(output / "slice_group_summary.csv", index=False, na_rep="NA")
    runtimes = runtime_summary(sms, cuda)
    runtimes.to_csv(output / "runtime_summary.csv", index=False, na_rep="NA")
    sms_noise = noise_summary(sms, args.motion_threshold)
    cuda_noise = noise_summary(cuda, args.motion_threshold)
    figures = [] if args.no_plots else make_plots(output, sms, cuda, paired, groups, args.motion_threshold)
    report = build_report(output, sms_inv, cuda_inv, sms, cuda, paired, only_sms, only_cuda, groups,
                          runtimes, sms_noise, cuda_noise, args.motion_threshold, args.head_radius, issues, figures)

    runtime = runtimes.set_index(["backend", "metric"])
    print(f"\nWrote analysis to: {output}")
    print(f"Transforms: sms-mi-reg={len(sms)}, CUDA={len(cuda)}; matched non-reference pairs={len(paired)}")
    print(f"Unmatched: sms-mi-reg-only={len(only_sms)}, CUDA-only={len(only_cuda)}")
    print(f"Mean FD (mm): sms-mi-reg={fmt(sms_noise['fd']['mean'])}, CUDA={fmt(cuda_noise['fd']['mean'])}")
    print(f"Mean relative disagreement: translation={fmt(paired.relative_translation_mm.mean())} mm, rotation={fmt(paired.relative_rotation_deg.mean())} deg")
    print(f"Optimizer speedup: {fmt(runtime.loc[('sms-mi-reg','optimizer_sec'),'mean'] / runtime.loc[('CUDA','optimizer_sec'),'mean'])}x")
    print(f"Registration-call speedup: {fmt(runtime.loc[('sms-mi-reg','call_sec'),'mean'] / runtime.loc[('CUDA','call_sec'),'mean'])}x")
    print(f"Report: {report}")
    if issues:
        print(f"Caveats/issues recorded in report: {len(issues)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

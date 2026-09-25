#!/usr/bin/env python3
"""Isolated benchmark for the production motion-dashboard rendering path."""

import argparse
import csv
import glob
import logging
import os
from pathlib import Path
import statistics
import sys
import tempfile
import time

import cv2
import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import SimpleITK as sitk


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
MOTION_MONITOR_DIR = REPOSITORY_ROOT / "apps" / "motion_monitor"
sys.path.insert(0, str(MOTION_MONITOR_DIR))

from generate_motion_plots import plot_motion_dashboard  # noqa: E402


PHASE_FIELDS = (
    "plot_data_extract_ms",
    "figure_create_ms",
    "plot_population_ms",
    "annotation_ms",
    "tsnr_panel_ms",
    "layout_ms",
    "savefig_ms",
    "cleanup_ms",
    "atomic_replace_ms",
    "published_read_ms",
    "plot_total_ms",
    "total_ms",
)


def one_match(directory, pattern, description):
    matches = sorted(glob.glob(os.path.join(directory, pattern)))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one {description} matching {pattern!r} in {directory}; "
            f"found {len(matches)}"
        )
    return matches[0]


def protocol_from_csv(path):
    basename = os.path.basename(path)
    prefix = "motionMonitor_data_"
    if not basename.startswith(prefix) or len(basename) < len(prefix) + 20:
        raise RuntimeError(f"Cannot parse protocol name from {basename}")
    return basename[len(prefix):-20]


def reference_volume(directory):
    candidates = sorted(
        path for path in glob.glob(os.path.join(directory, "*volume_0000*.nhdr"))
        if "_upsampled.nhdr" not in path
    )
    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected one non-upsampled volume-0 NHDR in {directory}; "
            f"found {len(candidates)}"
        )
    return sitk.GetArrayFromImage(sitk.ReadImage(candidates[0]))


def render_once(motion_df, tsnr_volume, protocol, expected_volumes, output_dir):
    temporary_path = os.path.join(output_dir, f".{protocol}.tmp.jpg")
    published_path = os.path.join(output_dir, f"{protocol}.jpg")
    timings = {}
    started_ns = time.perf_counter_ns()
    plot_motion_dashboard(
        motion_df,
        output_filename=temporary_path,
        protocol_name=protocol,
        threshold=0.3,
        num_expected_volumes=expected_volumes,
        num_moved_volumes=int(
            motion_df.loc[motion_df["Motion_flag"] != 0, "Volume_index"].nunique()
        ),
        livestream_enabled=False,
        tsnr_volume=tsnr_volume,
        tsnr_count=max(20, int(motion_df["Volume_index"].max()) + 1),
        tsnr_min_samples=10,
        tsnr_display_max=100.0,
        profile_timings=timings,
    )
    if not os.path.isfile(temporary_path):
        raise RuntimeError("Production plotting path did not create its temporary JPEG")

    phase_ns = time.perf_counter_ns()
    os.replace(temporary_path, published_path)
    timings["atomic_replace_ms"] = (time.perf_counter_ns() - phase_ns) / 1e6

    phase_ns = time.perf_counter_ns()
    image = cv2.imread(published_path)
    timings["published_read_ms"] = (time.perf_counter_ns() - phase_ns) / 1e6
    if image is None:
        raise RuntimeError("OpenCV could not read the published dashboard")
    timings["total_ms"] = (time.perf_counter_ns() - started_ns) / 1e6
    return timings, image.shape


def percentile(values, percentile_value):
    return float(np.percentile(np.asarray(values, dtype=float), percentile_value))


def summarize(rows):
    print("stage     samples  mean_ms  median_ms  min_ms  max_ms  p95_ms")
    for stage in ("early", "middle", "late"):
        values = [row["total_ms"] for row in rows if row["stage"] == stage]
        print(
            f"{stage:<9} {len(values):>7} "
            f"{statistics.mean(values):>8.3f} "
            f"{statistics.median(values):>10.3f} "
            f"{min(values):>7.3f} {max(values):>7.3f} "
            f"{percentile(values, 95):>7.3f}"
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark the unchanged production motion-dashboard plot/save path "
            "with early, middle, and late subsets of an archived motion CSV."
        )
    )
    parser.add_argument("acquisition_dir", help="Saved acquisition directory")
    parser.add_argument("--warmup", type=int, default=3, help="Warm-up renders per stage (default: 3)")
    parser.add_argument("--iterations", type=int, default=10, help="Measured renders per stage (default: 10)")
    parser.add_argument("--output", help="Output CSV (default: /tmp/motion_dashboard_profile_benchmark.csv)")
    args = parser.parse_args()

    acquisition_dir = os.path.abspath(args.acquisition_dir)
    if not os.path.isdir(acquisition_dir):
        parser.error(f"acquisition directory not found: {acquisition_dir}")
    if args.warmup < 0 or args.iterations < 1:
        parser.error("--warmup must be nonnegative and --iterations must be positive")

    motion_csv = one_match(acquisition_dir, "motionMonitor_data_*.csv", "motion CSV")
    motion_df = pd.read_csv(motion_csv)
    if motion_df.empty:
        raise RuntimeError(f"Motion CSV is empty: {motion_csv}")
    protocol = protocol_from_csv(motion_csv)
    tsnr_volume = reference_volume(acquisition_dir)
    expected_volumes = int(motion_df["Volume_index"].max()) + 1
    output_path = os.path.abspath(
        args.output
        or os.path.join(tempfile.gettempdir(), "motion_dashboard_profile_benchmark.csv")
    )

    stages = (
        ("early", max(1, round(len(motion_df) * 0.25))),
        ("middle", max(1, round(len(motion_df) * 0.50))),
        ("late", len(motion_df)),
    )
    rows = []
    logging.getLogger().setLevel(logging.WARNING)
    with tempfile.TemporaryDirectory(prefix="motion_dashboard_profile_") as directory:
        for stage, sample_count in stages:
            subset = motion_df.iloc[:sample_count]
            for _ in range(args.warmup):
                render_once(
                    subset,
                    tsnr_volume,
                    protocol,
                    expected_volumes,
                    directory,
                )
            for iteration in range(1, args.iterations + 1):
                timings, image_shape = render_once(
                    subset,
                    tsnr_volume,
                    protocol,
                    expected_volumes,
                    directory,
                )
                row = {
                    "stage": stage,
                    "iteration": iteration,
                    "motion_sample_count": sample_count,
                    "image_height": image_shape[0],
                    "image_width": image_shape[1],
                    **timings,
                }
                rows.append(row)

    fieldnames = (
        "stage",
        "iteration",
        "motion_sample_count",
        "image_height",
        "image_width",
        *PHASE_FIELDS,
    )
    with open(output_path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Acquisition: {acquisition_dir}")
    print(f"Motion rows: {len(motion_df)}")
    print(f"TSNR panel surrogate: volume-0 image, shape={tsnr_volume.shape}")
    print(f"Results: {output_path}")
    summarize(rows)


if __name__ == "__main__":
    main()

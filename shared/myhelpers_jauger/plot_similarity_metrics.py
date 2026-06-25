#!/usr/bin/env python3
"""
Helper script to plot similarity metrics from an Excel file.

Expected table headers:
"Instance" "MI_reference" "NMSE_reference" "SSIM_reference"
"MI_nomoco" "NMSE_nomoco" "SSIM_nomoco"
"MI_moco_initial" "MI_moco_final" "NMSE_moco" "SSIM_moco"

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Date of creation: October 14, 2025
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_ssim(df, save_dir=None):
    """
    Plot SSIM_reference, SSIM_nomoco, and SSIM_moco on the same plot.
    """
    instances = df["Instance"]

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the curves
    ax.plot(instances, df["SSIM_reference"], color="black", alpha=0.4, linewidth=2, label="No motion")
    ax.plot(instances, df["SSIM_nomoco"], color="blue", linewidth=2, marker='o', label="Head shake")
    ax.plot(instances, df["SSIM_moco"], color="orange", linewidth=2, marker='s', label="Head shake PMC")

    # Titles and labels with larger font
    ax.set_xlabel("Volume Index", fontsize=20)
    ax.set_ylabel("SSIM", fontsize=20)

    # Set axes
    ax.set_ylim(0.5, 1.05)
    ax.grid(True, linestyle="--", alpha=0.5)
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(axis='both', which='major', length=8, width=3, labelsize=16)
    ax.tick_params(axis='both', which='minor', length=4, width=1)
    ax.legend(fontsize=12, markerscale=1.1, loc="lower right")
    ax.set_xlim(left=0, right=len(instances[0:100]) + 0)

    # Save or show
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        outfile = os.path.join(save_dir, "SSIM_comparison_timeseries.png")
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        print(f"Saved: {outfile}")
    else:
        plt.show()

def plot_mutual_information(df, save_dir=None):
    """
    Plot Mutual Information (MI) comparisons:
      1. MI_nomoco vs MI_reference
      2. MI_moco_initial + MI_moco_final vs MI_reference
    """
    instances = df["Instance"]

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the curves
    ax.plot(instances, df["MI_reference"], color="black", alpha=0.4, linewidth=2, label="No motion")
    ax.plot(instances, df["MI_nomoco"], color="blue", linewidth=2, marker='o', label="Head shake")
    ax.plot(instances, df["MI_moco_final"], color="orange", linewidth=2, marker='s', label="Head shake PMC")

    # Titles and labels
    ax.set_xlabel("Instance", fontsize=20)
    ax.set_ylabel("MI", fontsize=20)

    # Axis limits
    ax.set_ylim(0.75, 3.00)
    ax.grid(True, linestyle="--", alpha=0.5)
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(axis='both', which='major', length=8, width=3, labelsize=16)
    ax.tick_params(axis='both', which='minor', length=4, width=1)
    ax.legend(fontsize=12, markerscale=1.1, loc="upper right")
    ax.set_xlim(left=0, right=len(instances) + 9)

    # Save or show
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        outfile = os.path.join(save_dir, "MI_comparison_timeseries.png")
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        print(f"Saved: {outfile}")
    else:
        plt.show()



def main():
    parser = argparse.ArgumentParser(description="Plot similarity metrics from Excel file.")
    parser.add_argument("--infile", required=True, help="Path to input .xlsx file")
    parser.add_argument("--outdir", default=None, help="Optional directory to save plots")
    args = parser.parse_args()

    if not os.path.exists(args.infile):
        raise FileNotFoundError(f"Input file not found: {args.infile}")

    df = pd.read_excel(args.infile)
    print(f"Loaded {len(df)} rows from {args.infile}")

    # Verify expected columns
    expected_cols = [
        "Instance", "MI_reference", "NMSE_reference", "SSIM_reference",
        "MI_nomoco", "NMSE_nomoco", "SSIM_nomoco",
        "MI_moco_initial", "MI_moco_final", "NMSE_moco", "SSIM_moco"
    ]
    missing = [c for c in expected_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing expected columns: {missing}")

    # Plot configurations
    plot_mutual_information(df, save_dir=args.outdir)
    plot_ssim(df, save_dir=args.outdir)

if __name__ == "__main__":
    main()

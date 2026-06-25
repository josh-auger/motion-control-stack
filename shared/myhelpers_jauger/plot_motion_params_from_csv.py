#!/usr/bin/env python3
"""
Title: plot_motion_params_from_csv.py

Description:
Helper script to re-plot the motion parameters from the motion-monitor data table (*.csv)

Usage:
    python plot_motion_params_from_csv.py --input_dir /path/to/csv_directory

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Date of creation: October 28, 2025
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_parameters(motion_parameters, series_name="", output_filename=""):
    """
    Plot combined motion parameters with rotations (x, y, z) on the top subplot
    and translations (x, y, z) on the bottom subplot.
    Rotations are converted from radians to degrees.
    """

    smooth_window = 15
    def smooth_signal(signal, window):
        if window < 2:
            return signal
        kernel = np.ones(window) / window
        return np.convolve(signal, kernel, mode='same')

    fig, axes = plt.subplots(3, 1, figsize=(16, 16), sharex=True)
    subplot_colors = ['b', 'g', 'r']  # x, y, z colors
    fig.suptitle(f"Motion Parameters : {series_name}", y=0.95, fontsize=14)

    # Extract parameter arrays for each axis
    param_array = np.array(motion_parameters)
    rotations = np.degrees(param_array[:, 0:3])  # First three are rotations
    translations = param_array[:, 3:6]           # Next three are translations
    x_values = np.arange(rotations.shape[0])

    # Compute magnitudes and smooth with moving average
    rotation_magnitude = np.linalg.norm(rotations, axis=1)
    translation_magnitude = np.linalg.norm(translations, axis=1)
    rotation_magnitude_smooth = smooth_signal(rotation_magnitude, smooth_window)
    translation_magnitude_smooth = smooth_signal(translation_magnitude, smooth_window)

    marker_style = dict(marker='o', linestyle='-', markersize=5, alpha=0.5)

    # --- Rotation subplot ---
    ax_rot = axes[0]
    for i, label in enumerate(['X', 'Y', 'Z']):
        ax_rot.plot(x_values, rotations[:, i], color=subplot_colors[i], label=f'{label}', **marker_style)
    ax_rot.set_ylabel('Rotation (deg)', fontsize=24)
    ax_rot.grid(True, linestyle='-', linewidth=0.5, color='gray', alpha=0.5)
    ax_rot.legend(fontsize=20, markerscale=2, loc='upper right')
    # ax_rot.set_xlim(left=0, right=len(x_values)+20)
    ax_rot.tick_params(axis='both', labelsize=20)
    ax_rot.set_ylim([-8.5, 8])
    ax_rot.set_yticks(np.arange(-8, 9, 2))

    # --- Translation subplot ---
    ax_trans = axes[1]
    for i, label in enumerate(['X', 'Y', 'Z']):
        ax_trans.plot(x_values, translations[:, i], color=subplot_colors[i], label=f'{label}', **marker_style)
    # ax_trans.set_xlabel('Slice Index', fontsize=24)
    ax_trans.set_ylabel('Translation (mm)', fontsize=24)
    ax_trans.grid(True, linestyle='-', linewidth=0.5, color='gray', alpha=0.5)
    ax_trans.legend(fontsize=20, markerscale=2, loc='upper right')
    # ax_trans.set_xlim(left=0, right=len(x_values)+20)
    ax_trans.tick_params(axis='both', labelsize=20)
    ax_trans.set_ylim([-6.5, 6.5])
    ax_trans.set_yticks(np.arange(-6, 7, 2))

    # --- Magnitude subplot ---
    ax_mag = axes[2]
    ax_mag.plot(x_values, rotation_magnitude_smooth, color='blue', label='Rotation (deg)', linewidth=2)
    ax_mag.plot(x_values, translation_magnitude_smooth, color='red', label='Translation (mm)', linewidth=2)
    ax_mag.set_xlabel('Acquisition Group Index', fontsize=22)
    ax_mag.set_ylabel('Magnitude (deg / mm)', fontsize=22)
    ax_mag.grid(True, linestyle='-', linewidth=0.5, color='gray', alpha=0.5)
    ax_mag.legend(fontsize=16, loc='upper right')
    ax_mag.tick_params(axis='both', labelsize=20, width=2, length=8)
    ax_mag.set_ylim([-0.5, 8.5])
    ax_mag.set_yticks(np.arange(0, 9, 2))

    ax_mag.set_xlim(left=0, right=len(x_values) + 2)
    ax_trans.set_xticks(np.arange(0, len(x_values) + 1, 160))

    # --- Make spines and ticks thicker ---
    for ax in axes:
        for spine in ax.spines.values():
            spine.set_linewidth(2)  # Thicker outer bounding box
        ax.tick_params(width=2, length=10)  # Thicker ticks

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave room for the title
    plt.savefig(output_filename, dpi=300)
    print(f"Combined parameters plot saved as: {output_filename}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot motion parameters from a .csv file.")
    parser.add_argument("--csvfile", required=True, help="Path to input CSV file containing motion parameters.")
    parser.add_argument("--series_name", default="", help="Name of the motion series (optional).")
    parser.add_argument("--output", default="", help="Output filename for the plot (optional).")
    args = parser.parse_args()

    if not os.path.exists(args.csvfile):
        print(f"Error: CSV file not found: {args.csvfile}")
        exit(1)

    # Read the CSV
    df = pd.read_csv(args.csvfile)

    # Expected columns
    expected_cols = [
        "X_rotation(rad)", "Y_rotation(rad)", "Z_rotation(rad)",
        "X_translation(mm)", "Y_translation(mm)", "Z_translation(mm)",
        "Displacement(mm)", "Cumulative_displacement(mm)",
        "Volume_number", "Motion_flag"
    ]

    for col in expected_cols:
        if col not in df.columns:
            print(f"Error: Missing expected column in CSV: {col}")
            exit(1)

    # Extract motion parameter subset (first 6 columns)
    motion_parameters = df[[
        "X_rotation(rad)", "Y_rotation(rad)", "Z_rotation(rad)",
        "X_translation(mm)", "Y_translation(mm)", "Z_translation(mm)"
    ]].to_numpy()

    # Default output filename
    if args.output == "":
        base = os.path.splitext(os.path.basename(args.csvfile))[0]
        args.output = f"{base}_motion_params.png"

    # Plot
    plot_parameters(
        motion_parameters,
        series_name=args.series_name,
        output_filename=args.output
    )


if __name__ == "__main__":
    main()
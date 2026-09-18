#!/usr/bin/env python3

# Title: generate_motion_plots.py

# Description:
# Helper script to generate plots of motion parameters and framewise displacement.

# Created on: December 2025
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import textwrap
import logging
import csv
import os

try:
    from .running_tsnr import calculate_mean_tsnr, get_tsnr_display_slices
except ImportError:  # Support execution from the motion_monitor application directory.
    from running_tsnr import calculate_mean_tsnr, get_tsnr_display_slices


DASHBOARD_PIXEL_WIDTH = 1600
DASHBOARD_PIXEL_HEIGHT = 1050
DASHBOARD_DPI = 100


def add_tsnr_dashboard_panel(
    fig,
    subplot_spec,
    *,
    tsnr_volume=None,
    tsnr_count=0,
    min_samples=20,
    display_max=100.0,
):
    """Add the shared-orientation raw TSNR panel to a dashboard figure."""
    if tsnr_count < min_samples:
        axis = fig.add_subplot(subplot_spec)
        axis.set_title("raw TSNR", pad=4)
        axis.axis("off")
        message = axis.text(
            0.5,
            0.5,
            f"Waiting for TSNR\nN = {tsnr_count} / {min_samples}",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color="0.3",
        )
        return {"state": "waiting", "axis": axis, "text": message}

    if tsnr_volume is None:
        axis = fig.add_subplot(subplot_spec)
        axis.set_title("raw TSNR", pad=4)
        axis.axis("off")
        message = axis.text(
            0.5,
            0.5,
            f"TSNR unavailable\nN = {tsnr_count}",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color="0.3",
        )
        return {"state": "unavailable", "axis": axis, "text": message}

    display_slices = get_tsnr_display_slices(
        tsnr_volume,
        display_max=display_max,
    )
    mean_tsnr = calculate_mean_tsnr(tsnr_volume)

    title_axis = fig.add_subplot(subplot_spec)
    title_axis.set_title("raw TSNR", pad=4)
    title_axis.axis("off")
    title_axis.set_zorder(-1)

    panel_grid = subplot_spec.subgridspec(
        2,
        3,
        width_ratios=[1, 1, 0.08],
        hspace=0.015,
        wspace=0.025,
    )
    image_axes = [
        fig.add_subplot(panel_grid[row, column])
        for row in range(2)
        for column in range(2)
    ]
    image_artist = None
    for axis, display_slice in zip(image_axes, display_slices):
        image_artist = axis.imshow(
            display_slice,
            cmap="gray",
            vmin=0.0,
            vmax=display_max,
            interpolation="nearest",
            aspect="equal",
        )
        axis.axis("off")

    overlay = image_axes[2].text(
        0.04,
        0.04,
        f"N = {tsnr_count}\nMean TSNR = {mean_tsnr:.1f}",
        transform=image_axes[2].transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        color="white",
        bbox=dict(facecolor="black", alpha=0.55, edgecolor="none", pad=2),
    )

    colorbar_axis = fig.add_subplot(panel_grid[:, 2])
    colorbar = fig.colorbar(image_artist, cax=colorbar_axis, orientation="vertical")
    colorbar.ax.tick_params(labelsize=8)

    return {
        "state": "ready",
        "title_axis": title_axis,
        "image_axes": image_axes,
        "image_artist": image_artist,
        "colorbar": colorbar,
        "overlay": overlay,
        "mean_tsnr": mean_tsnr,
    }


def load_motion_csv(csv_filepath):
    """
    Load motion tracking CSV into a pandas DataFrame.
    """
    df = pd.read_csv(csv_filepath)

    required_columns = [
        "X_rotation(rad)", "Y_rotation(rad)", "Z_rotation(rad)",
        "X_translation(mm)", "Y_translation(mm)", "Z_translation(mm)",
        "Displacement(mm)", "Cumulative_displacement(mm)", "Motion_flag"
    ]

    missing = set(required_columns) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")

    return df


def plot_parameters_combined(motion_df, output_filename="", protocol_name=""):
    """
    motion_df = pandas DataFrame loaded from the motion CSV file
    Plot combined motion parameters with rotations (x, y, z) on one subplot and translations (x, y, z) on another.
    Rotations are converted from radians to degrees.
    """

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))
    subplot_colors = ['b', 'g', 'r']  # x, y, z colors
    fig.suptitle(f"Motion Parameters: {protocol_name}", y=0.98)

    # Extract parameter arrays from motion dataframe
    rotations_rad = motion_df[["X_rotation(rad)", "Y_rotation(rad)", "Z_rotation(rad)"]].to_numpy()
    translations = motion_df[["X_translation(mm)", "Y_translation(mm)", "Z_translation(mm)"]].to_numpy()
    # Convert rotations to degrees
    rotations = np.degrees(rotations_rad)

    # Rotation plot
    ax_rot = axes[0]
    for i, label in enumerate(['X', 'Y', 'Z']):
        ax_rot.plot(rotations[:, i], marker='o', linestyle='-', color=subplot_colors[i],
                    alpha=0.6, label=f'{label}-axis')
    ax_rot.set_title('Rotations')
    ax_rot.set_xlabel('Acquisition Index')
    ax_rot.set_ylabel('Rotation (deg)')
    ax_rot.grid(True, linestyle='-', linewidth=0.5, color='gray', alpha=0.5)
    ax_rot.set_xlim(left=0)
    ax_rot.legend(loc='upper left')

    # Translation plot
    ax_trans = axes[1]
    for i, label in enumerate(['X', 'Y', 'Z']):
        ax_trans.plot(translations[:, i], marker='o', linestyle='-', color=subplot_colors[i],
                      alpha=0.6, label=f'{label}-axis')
    ax_trans.set_title('Translations')
    ax_trans.set_xlabel('Acquisition Index')
    ax_trans.set_ylabel('Translation (mm)')
    ax_trans.grid(True, linestyle='-', linewidth=0.5, color='gray', alpha=0.5)
    ax_trans.set_xlim(left=0)
    ax_trans.legend(loc='upper left')

    plt.tight_layout()
    plt.savefig(output_filename)
    logging.info(f"Parameters plot saved as : {output_filename}")
    plt.ion()
    plt.close(fig)


def plot_displacements(motion_df, output_filename="", protocol_name="", threshold=None, num_expected_volumes=None, num_moved_volumes=None):
    displacements = motion_df["Displacement(mm)"].to_numpy()
    cumulative = motion_df["Cumulative_displacement(mm)"].to_numpy()
    volume_index = motion_df["Volume_index"]
    motion_flags = motion_df["Motion_flag"]

    fig, axs = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [3.5, 1]})

    # Plot framewise displacements
    axs[0].plot(displacements, marker='o', alpha=0.7, label="Framewise displacement")
    if threshold is not None:
        axs[0].axhline(threshold, color='r', linestyle='--', label=f"Threshold = {threshold} mm")

    axs[0].set_title(f"Framewise Displacements: {protocol_name}")
    axs[0].set_xlabel("Index")
    axs[0].set_ylabel("Displacement (mm)")
    axs[0].legend(loc='upper left')
    axs[0].grid(True)

    # Plot motion summary and displacements boxplot
    num_volumes = max(volume_index)
    if num_moved_volumes is None:
        logging.info(f"Number of motion-corrupt volumes not provided. Treating each motion flag separately.")
        num_moved_volumes = motion_flags.sum()
    num_motion_free_volumes = num_volumes - num_moved_volumes
    counts = [num_motion_free_volumes, num_moved_volumes]
    if num_expected_volumes is None:
        logging.info(f"Numer of expected volumes (sequence repetitions) not provided. Using current volume count.")
        num_expected_volumes = num_volumes

    text = (
        f"Number of acquisitions: {len(displacements)}\n"
        f"Motion flags: {motion_flags.sum()}\n"
        f"Cumulative displacement (mm): {cumulative[-1]:.3f}\n"
        f"Number of volumes: {num_volumes}\n"
        f"Motion-free volumes: {num_motion_free_volumes}\n"
        f"Motion-corrupt volumes: {num_moved_volumes}"
    )
    axs[1].text(0.5, 0.9, text, transform=axs[1].transAxes, ha='center', va='center', multialignment='left',
                bbox=dict(facecolor='white', alpha=0.8))
    axs[1].bar([f"Motion-free", f"Motion-corrupt"], counts, color=["tab:blue", "tab:orange"])
    ymax = num_expected_volumes
    axs[1].set_ylim(0, ymax * 1.4)
    axs[1].set_ylabel("N volumes")
    axs[1].set_title("Volume motion summary")
    for i, val in enumerate([num_motion_free_volumes, num_moved_volumes]):
        axs[1].text(i, val + 0.5, str(val), ha='center', va='bottom')
    axs[1].grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_filename)
    logging.info(f"Displacements plot saved as : {output_filename}")
    plt.close(fig)


def plot_cumulative_displacement(motion_df, output_filename="", threshold=None):
    cumulative = motion_df["Cumulative_displacement(mm)"].to_numpy()

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(cumulative, marker='o', alpha=0.7, label="Cumulative displacement")

    if threshold is not None:
        ax.plot(np.arange(len(cumulative)), threshold * np.arange(len(cumulative)),
                linestyle='--', color='r', label="Acceptable accumulation")

    ax.set_title(f"Cumulative Displacement")
    ax.set_xlabel("Acquisition Instance")
    ax.set_ylabel("Cumulative Displacement (mm)")
    ax.legend(loc='upper left')
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(output_filename)
    logging.info(f"Cumulative displacement plot saved as : {output_filename}")
    plt.close(fig)


def plot_motion_dashboard(motion_df, output_filename="", protocol_name="", threshold=None, num_expected_volumes=None,
                          num_moved_volumes=None, host_ip=None, host_port=None, livestream_enabled=False,
                          tsnr_volume=None, tsnr_count=0, tsnr_min_samples=20,
                          tsnr_display_max=100.0):
    """
    motion_df = pandas DataFrame loaded from the motion CSV file
    Generate a combined motion report figure (framewise displacement, motion summary, motion parameters) to push to a
    live stream of ongoing motion results.
    """

    # --- Extract data from  motion DataFrame ---
    rotations_rad = motion_df[["X_rotation(rad)", "Y_rotation(rad)", "Z_rotation(rad)"]].to_numpy()
    rotations = np.degrees(rotations_rad)
    translations = motion_df[["X_translation(mm)", "Y_translation(mm)", "Z_translation(mm)"]].to_numpy()
    displacements = motion_df["Displacement(mm)"].to_numpy()
    cumulative = motion_df["Cumulative_displacement(mm)"].to_numpy()
    volume_index = motion_df["Volume_index"]
    motion_flags = motion_df["Motion_flag"]

    # --- Figure layout and axes ---
    target_pixel_width = DASHBOARD_PIXEL_WIDTH
    target_pixel_height = DASHBOARD_PIXEL_HEIGHT
    target_dpi = DASHBOARD_DPI
    fig = plt.figure(figsize=(target_pixel_width / target_dpi, target_pixel_height / target_dpi),dpi=target_dpi)
    gs = fig.add_gridspec(
        nrows=3, ncols=2,
        height_ratios=[1.0, 1.0, 1.15],
        width_ratios=[3.7, 1.25],
        hspace=0.30, wspace=0.12,
    )
    ax_trans = fig.add_subplot(gs[0, 0])         # translations (top left)
    ax_status = fig.add_subplot(gs[0, 1])        # motion summary (top right)
    ax_status.axis("off")
    ax_rot = fig.add_subplot(gs[1, 0])           # rotations (middle left)
    ax_sum = fig.add_subplot(gs[1, 1])           # volume tracker (middle right)
    ax_disp = fig.add_subplot(gs[2, 0])          # framewise displacement (bottom left)

    # --- Translation parameters (top left) ---
    subplot_colors = ['b', 'g', 'r']
    for i, label in enumerate(['X', 'Y', 'Z']):
        ax_trans.plot(
            translations[:, i],
            linestyle='-',
            color=subplot_colors[i],
            alpha=0.75,
            label=f'{label}'
        )

    ax_trans.set_title("Motion Parameters", pad=4)
    ax_trans.set_ylabel("Translation (mm)")
    ax_trans.grid(True, linestyle='-', linewidth=0.5, color='gray', alpha=0.5)
    ax_trans.set_xlim(left=0)
    ax_trans.legend(loc='upper left')

    # --- Rotation parameters (middle left) ---
    for i, label in enumerate(['X', 'Y', 'Z']):
        ax_rot.plot(
            rotations[:, i],
            linestyle='-',
            color=subplot_colors[i],
            alpha=0.75,
            label=f'{label}'
        )

    # ax_rot.set_xlabel("Acquisition Index")
    ax_rot.set_ylabel("Rotation (deg)")
    ax_rot.grid(True, linestyle='-', linewidth=0.5, color='gray', alpha=0.5)
    ax_rot.set_xlim(left=0)
    ax_rot.legend(loc='upper left')

    # --- Framewise displacement (bottom left) ---
    ax_disp.plot(displacements, marker='o', alpha=0.7, label="Displacement")
    if threshold is not None:
        ax_disp.axhline(threshold, color='r', linestyle='--', label=f"Threshold = {threshold} mm")
    ax_disp.set_title("Framewise Displacement", pad=4)
    ax_disp.set_xlabel("Acquisition Index")
    ax_disp.set_ylabel("Displacement (mm)")
    ax_disp.legend(loc='upper left')
    ax_disp.grid(True, alpha=0.4)
    ax_disp.set_xlim(left=0)

    # --- Volume tracker (middle right) ---
    num_volumes = int(max(volume_index))
    if num_moved_volumes is None:
        num_moved_volumes = int(motion_flags.sum())

    num_motion_free_volumes = num_volumes - num_moved_volumes
    counts = [num_motion_free_volumes, num_moved_volumes]

    if num_expected_volumes is None:
        num_expected_volumes = num_volumes

    text = (
        f"Volumes: {num_volumes}/{num_expected_volumes}\n"
        f"Motion-free: {num_motion_free_volumes}\n"
        f"Motion-corrupt: {num_moved_volumes}"
    )

    ax_sum.text(
        0.5, 0.95, text,
        transform=ax_sum.transAxes,
        ha='center', va='top',
        multialignment='left',
        fontsize=10,
        bbox=dict(facecolor='white', alpha=0.85),
        clip_on=True
    )
    ax_sum.bar(["Motion-free", "Motion-corrupt"], counts, color=["tab:blue", "tab:orange"])
    ax_sum.set_ylim(0, num_expected_volumes * 1.4)
    ax_sum.set_ylabel("N volumes")
    ax_sum.set_title("Volume Tracker", pad=4)
    for i, val in enumerate(counts):
        ax_sum.text(i, val + 0.5, str(val), ha='center', va='bottom')
    ax_sum.grid(axis='y', alpha=0.3)

    # --- Status panel (top right) ---
    wrapped_protocol = "\n".join(textwrap.wrap(protocol_name, width=36)) \
        if protocol_name else "N/A"

    if livestream_enabled:
        stream_info = (
            f"Host: {host_ip}\n"
            f"Port: {host_port}\n"
            f"http://{host_ip}:{host_port}/stream.mjpg\n\n"
        )
    else:
        stream_info = (
            "Livestream: Disabled\n\n"
        )

    status_text = (
        stream_info +
        f"Protocol name:\n{wrapped_protocol}\n\n"
        f"Acquisitions: {len(displacements)}\n"
        f"Motion flags: {motion_flags.sum()}\n"
        f"Cumulative disp (mm): {cumulative[-1]:.3f}\n"
        f"Current volumes: {num_volumes} / {num_expected_volumes}\n"
        f"Motion-free volumes: {num_motion_free_volumes}\n"
        f"Motion-corrupt volumes: {num_moved_volumes}"
    )

    ax_status.text(
        0.01, 0.99,
        status_text,
        transform=ax_status.transAxes,
        ha="left", va="top",
        multialignment="left",
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray")
    )
    ax_status.set_title("MOTION MONITOR", pad=4)

    # --- Raw TSNR (bottom right) ---
    try:
        add_tsnr_dashboard_panel(
            fig,
            gs[2, 1],
            tsnr_volume=tsnr_volume,
            tsnr_count=tsnr_count,
            min_samples=tsnr_min_samples,
            display_max=tsnr_display_max,
        )
    except Exception:
        logging.exception("Failed to render TSNR dashboard panel.")
        add_tsnr_dashboard_panel(
            fig,
            gs[2, 1],
            tsnr_volume=None,
            tsnr_count=max(tsnr_count, tsnr_min_samples),
            min_samples=tsnr_min_samples,
            display_max=tsnr_display_max,
        )

    # --- Finalize figure ---
    fig.subplots_adjust(
        left=0.05,
        right=0.98,
        bottom=0.06,
        top=0.96,
        wspace=0.12,
        hspace=0.30
    )
    if output_filename:
        plt.savefig(output_filename, dpi=target_dpi, bbox_inches=None,pad_inches=0)
        logging.info(f"Motion dashboard figure saved as : {output_filename}")
    plt.show(block=False)
    plt.close(fig)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate motion plots from motion CSV")
    parser.add_argument("motion_csv", help="Path to motion tracking CSV")
    parser.add_argument("--series-name", default="", help="Series name for plot titles")
    parser.add_argument("--outdir", default=".", help="Output directory for plots")
    parser.add_argument("--threshold", type=float, default=None, help="Motion threshold (mm)")
    args = parser.parse_args()

    motion_df = load_motion_csv(args.motion_csv)

    plot_parameters_combined(
        motion_df,
        output_filename=os.path.join(args.outdir, "motion_parameters_combined.png")
    )

    plot_displacements(
        motion_df,
        output_filename=os.path.join(args.outdir, "framewise_displacement.png"),
        threshold=args.threshold
    )

    plot_cumulative_displacement(
        motion_df,
        output_filename=os.path.join(args.outdir, "cumulative_displacement.png"),
        threshold=args.threshold
    )

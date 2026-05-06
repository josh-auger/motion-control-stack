#!/usr/bin/env python3
"""
Title: display_tsnr_results.py

Description:
Generate a mid-slice comparison figure for input volume, mean signal, standard deviation, and tSNR map.
The resulting figure is saved in the same directory as the tSNR map file.

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Date of creation: October 17, 2025
"""

import os
import argparse
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np


def add_cross_hair(ax, x, y, color='red', linewidth=0.5):
    """Draw cross-hairs at (x, y) on the given axis."""
    ax.axhline(y, color=color, linewidth=linewidth)
    ax.axvline(x, color=color, linewidth=linewidth)

def preprocess_for_display(volume):
    """Flip nifti image volume about y-axis for display."""
    volume = np.flip(volume, axis=1)
    return volume


def display_tsnr_components(reference_img, mean_signal, std_signal, tsnr_map, save_dir):
    """Display mid-slice comparison for input volume, mean signal, std-dev, and tSNR with colorbars."""
    # Load reference image volume
    volume = sitk.GetArrayFromImage(reference_img)
    # Preprocess arrays for display orientation
    volume = preprocess_for_display(volume)
    mean_signal = preprocess_for_display(mean_signal)
    std_signal = preprocess_for_display(std_signal)
    tsnr_map = preprocess_for_display(tsnr_map)

    fig, axes = plt.subplots(1, 4, figsize=(22, 8), constrained_layout=True)

    slice_idx = 31 #volume.shape[0] // 2    # slice 30 (or 31) for 10/7 motion scans, slice 27 for 10/14 no-motion reference
    center_x = volume.shape[2] // 2
    center_y = volume.shape[1] // 2

    # Compute mean tSNR of axial slice
    tsnr_slice = tsnr_map[slice_idx, :, :]
    mean_tsnr_slice = np.mean(tsnr_slice[tsnr_slice > 0])
    print(f"Mean tSNR for slice {slice_idx}: {mean_tsnr_slice:.2f}")

    titles = [f"Input Volume (slice {slice_idx})", "Mean Signal", "Standard Deviation", "tSNR"]
    images = [volume, mean_signal, std_signal, tsnr_map]
    cmaps = ['gray', 'viridis', 'viridis', 'plasma']

    for ax, title, vol, cmap in zip(axes, titles, images, cmaps):
        # Cap the tSNR map at 100 for consistent visualization
        if title == "tSNR":
            im = ax.imshow(vol[slice_idx, :, :], cmap=cmap, vmin=0, vmax=100)
            # im = ax.imshow(vol[slice_idx, :, :], cmap=cmap)
        else:
            im = ax.imshow(vol[slice_idx, :, :], cmap=cmap)
        ax.set_title(title, fontsize=14)
        add_cross_hair(ax, center_x, center_y)
        ax.set_xticks([])
        ax.set_yticks([])

        # Orientation labels
        ax.text(0.5, 0.03, 'P', transform=ax.transAxes, ha='center', va='bottom', fontsize=14, color='white')
        ax.text(0.5, 0.97, 'A', transform=ax.transAxes, ha='center', va='top', fontsize=14, color='white')
        ax.text(0.03, 0.5, 'R', transform=ax.transAxes, ha='left', va='center', fontsize=14, color='white')
        ax.text(0.97, 0.5, 'L', transform=ax.transAxes, ha='right', va='center', fontsize=14, color='white')

        # Add colorbar
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    plt.show()

    # Save figure to same directory as tsnr map
    plot_filename = os.path.join(save_dir, "figure_tsnr_results.png")
    fig.savefig(plot_filename)
    print(f"Figure saved: {plot_filename}")


def display_tsnr_slices(tsnr_map, save_dir, slice_indices=None, vmax=None):
    """
    Display 4 axial slices from the tSNR volume with a consistent color scale.

    Parameters:
        tsnr_map (np.ndarray): 3D tSNR volume (z, y, x)
        save_dir (str): Directory to save the figure
        slice_indices (list or None): List of 4 slice indices. If None, auto-select.
        vmax (float): Upper bound for color scale (default=100)
    """
    # Preprocess orientation
    tsnr_map = preprocess_for_display(tsnr_map)
    nz = tsnr_map.shape[0]
    center_x = tsnr_map.shape[2] // 2
    center_y = tsnr_map.shape[1] // 2

    # Auto-select 4 evenly spaced slices if not provided
    if slice_indices is None:
        slice_indices = [
            int(nz * 0.2),
            int(nz * 0.4),
            int(nz * 0.6),
            int(nz * 0.8),
        ]

    fig, axes = plt.subplots(1, 4, figsize=(20, 6), constrained_layout=True)
    # Use consistent color scaling
    vmin = 0
    if vmax is None:
        # vmax = np.percentile(tsnr_map[tsnr_map > 0], 99)
        vmax = 100

    for ax, slice_idx in zip(axes, slice_indices):
        slice_data = tsnr_map[slice_idx, :, :]

        im = ax.imshow(slice_data, cmap='plasma', vmin=vmin, vmax=vmax)
        ax.set_title(f"tSNR Slice {slice_idx}", fontsize=14)
        add_cross_hair(ax, center_x, center_y)
        ax.set_xticks([])
        ax.set_yticks([])

        # Orientation labels
        ax.text(0.5, 0.03, 'P', transform=ax.transAxes, ha='center', va='bottom', fontsize=14, color='white')
        ax.text(0.5, 0.97, 'A', transform=ax.transAxes, ha='center', va='top', fontsize=14, color='white')
        ax.text(0.03, 0.5, 'R', transform=ax.transAxes, ha='left', va='center', fontsize=14, color='white')
        ax.text(0.97, 0.5, 'L', transform=ax.transAxes, ha='right', va='center', fontsize=14, color='white')

    # Single shared colorbar
    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("tSNR", fontsize=12)

    plt.show()

    # Save figure
    plot_filename = os.path.join(save_dir, "figure_tsnr_axialslices.png")
    fig.savefig(plot_filename)
    print(f"Figure saved: {plot_filename}")



def main():
    parser = argparse.ArgumentParser(description="Display tSNR results for a given dataset.")
    parser.add_argument("--ref", required=True, help="Path to reference volume file (e.g., NIfTI)")
    parser.add_argument("--mean", required=True, help="Path to mean signal file")
    parser.add_argument("--std", required=True, help="Path to standard deviation signal file")
    parser.add_argument("--tsnr", required=True, help="Path to tSNR map file")
    args = parser.parse_args()

    # Load images
    reference_img = sitk.ReadImage(args.ref)
    mean_signal = sitk.GetArrayFromImage(sitk.ReadImage(args.mean))
    std_signal = sitk.GetArrayFromImage(sitk.ReadImage(args.std))
    tsnr_map = sitk.GetArrayFromImage(sitk.ReadImage(args.tsnr))

    save_dir = os.path.dirname(args.tsnr)

    display_tsnr_components(reference_img, mean_signal, std_signal, tsnr_map, save_dir)
    display_tsnr_slices(tsnr_map, save_dir)

    mean_tsnr = np.mean(tsnr_map[tsnr_map > 0])
    print(f"Mean tSNR of volume: {mean_tsnr:.4f}")


if __name__ == "__main__":
    main()
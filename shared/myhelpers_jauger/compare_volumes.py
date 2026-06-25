
# Title: compare_volumes.py

# Description:
# Read in two volumes (NRRD file format)
# Compare them by computing the DICE similarity score
# Also compare by simple subtraction of the two volumes

# Created on: April 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import argparse
import os
import sys
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from skimage.metrics import structural_similarity as ssim
import csv
import datetime


def print_image_stats(itk_image, title):
    stats = sitk.StatisticsImageFilter()
    stats.Execute(itk_image)

    min_pixel_value = stats.GetMinimum()
    max_pixel_value = stats.GetMaximum()
    mean_pixel_value = stats.GetMean()

    origin = itk_image.GetOrigin()
    spacing = itk_image.GetSpacing()
    direction = itk_image.GetDirection()

    print(f"Statistics for {title}:")
    print(f"  Min Pixel Value: {min_pixel_value}")
    print(f"  Max Pixel Value: {max_pixel_value}")
    print(f"  Mean Pixel Value: {mean_pixel_value}")
    print(f"  Origin: {origin}")
    print(f"  Spacing: {spacing}")
    print(f"  Direction Cosine Matrix:")
    print(direction)

def binarize_image_otsu(itk_image):
    otsu_filter = sitk.OtsuThresholdImageFilter()
    otsu_filter.SetInsideValue(0)
    otsu_filter.SetOutsideValue(1)

    binary_itk_image = otsu_filter.Execute(itk_image)
    threshold_value = otsu_filter.GetThreshold()
    return binary_itk_image, threshold_value

def calculate_dice(image1, image2):
    binary_image1, threshold_value = binarize_image_otsu(image1)  # get Otsu threshold for reference volume
    print(f"Otsu threshold value: {threshold_value:.6f}")
    max_pixel_value_image2 = np.max(sitk.GetArrayFromImage(image2))
    # apply same Otsu threshold to target volume for comparison consistency
    binary_image2 = sitk.BinaryThreshold(
        image2,
        lowerThreshold=int(threshold_value),
        upperThreshold=int(max_pixel_value_image2),
        insideValue=1,
        outsideValue=0,
    )
    binary_volume_array1 = np.flipud(sitk.GetArrayFromImage(binary_image1))
    binary_volume_array2 = np.flipud(sitk.GetArrayFromImage(binary_image2))

    intersection = np.sum(binary_volume_array1[binary_volume_array2 == 1])
    total_voxels = np.sum(binary_volume_array1) + np.sum(binary_volume_array2)
    dice = 2 * intersection / total_voxels
    return dice

def calculate_nmse(volume1, volume2):
    # Ensure same shape
    if volume1.shape != volume2.shape:
        raise ValueError("Volumes must have the same shape for NMSE calculation.")

    mse = np.mean((volume1 - volume2) ** 2)
    norm_factor = np.mean(volume1 ** 2)
    nmse = mse / norm_factor
    return nmse

def calculate_ssim(volume1, volume2):
    # Ensure same shape
    if volume1.shape != volume2.shape:
        raise ValueError("Volumes must have the same shape for SSIM calculation.")

    # Normalize volumes to [0, 1] range for numerical stability
    v1 = (volume1 - np.min(volume1)) / (np.max(volume1) - np.min(volume1))
    v2 = (volume2 - np.min(volume2)) / (np.max(volume2) - np.min(volume2))

    # Compute SSIM slice-by-slice (since skimage.ssim is 2D)
    ssim_values = []
    for i in range(v1.shape[0]):
        ssim_val, _ = ssim(v1[i, :, :], v2[i, :, :], data_range=1.0, full=True)
        ssim_values.append(ssim_val)
    return np.mean(ssim_values)

def subtract_volumes(volume1, volume2):
    return (volume1 - volume2)

def add_cross_hair(x, y, color='red', linewidth=0.5):
    plt.axhline(y, color=color, linewidth=linewidth)
    plt.axvline(x, color=color, linewidth=linewidth)

def preprocess_for_display(volume):
    """Flip nifti image volume about y-axis for display."""
    volume = np.flip(volume, axis=1)
    return volume

def display_subtraction_results(volume1, volume2, subtracted_volume, dice_score, nmse_value, ssim_value):
    fig, axes = plt.subplots(1, 3, figsize=(16, 8))
    plt.subplots_adjust(wspace=0.05)

    slice_idx = volume1.shape[0] // 2
    center_x = volume1.shape[2] // 2
    center_y = volume1.shape[1] // 2

    titles = [f"Ref Volume 1, slice {slice_idx}", f"Volume 2", f"Normalized subtraction (ref-target)\nDice = {dice_score:.4f}, NMSE = {nmse_value:.4f}, SSIM = {ssim_value:.4f}"]
    images = [volume1, volume2, subtracted_volume]
    cmaps = ['gray', 'gray', 'coolwarm']

    for ax, title, vol, cmap in zip(axes, titles, images, cmaps):
        if title.startswith("Normalized subtraction"):
            vmax = np.nanmax(np.abs(vol))
            im = ax.imshow(vol[slice_idx, :, :], cmap=cmap, vmin=-1, vmax=1)
        else:
            im = ax.imshow(vol[slice_idx, :, :], cmap=cmap)
        ax.set_title(title)
        ax.axhline(center_y, color='red', linewidth=1.25)
        ax.axvline(center_x, color='red', linewidth=1.25)
        ax.set_xticks([])
        ax.set_yticks([])

        # Anatomic labels inside the image bounds
        ax.text(0.5, 0.03, 'P', transform=ax.transAxes,
                ha='center', va='bottom', fontsize=20, color='white')
        ax.text(0.5, 0.97, 'A', transform=ax.transAxes,
                ha='center', va='top', fontsize=20, color='white')
        ax.text(0.03, 0.5, 'R', transform=ax.transAxes,
                ha='left', va='center', fontsize=20, color='white')
        ax.text(0.97, 0.5, 'L', transform=ax.transAxes,
                ha='right', va='center', fontsize=20, color='white')

    bbox = axes[2].get_position()
    cbar_ax = fig.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, bbox.height])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(width=3, length=8, labelsize=20)
    cbar.set_ticks(np.arange(-1.0, 1.01, 0.5))

    # output_filename = f"scan000_YYYYMMDD_volume_comparison_0000_to_0001.png"
    # plt.savefig(output_filename, dpi=300)
    # print(f"Combined parameters plot saved as: {output_filename}")
    plt.show()


def compute_volume_comparison(file1, file2, show_subtraction=False):
    # Load images
    image1 = sitk.ReadImage(file1)
    image2 = sitk.ReadImage(file2)
    print_image_stats(image1, "Ref Volume 1")
    print_image_stats(image2, "Volume 2")

    # Convert images to numpy arrays, and flip along y-axis
    volume_array1 = np.flipud(sitk.GetArrayFromImage(image1))
    volume_array1 = preprocess_for_display(volume_array1)
    volume_array2 = np.flipud(sitk.GetArrayFromImage(image2))
    volume_array2 = preprocess_for_display(volume_array2)

    # Compute Dice score, NMSE, and SSIM
    dice_score = calculate_dice(image1, image2)
    nmse_value = calculate_nmse(volume_array1, volume_array2)
    ssim_value = calculate_ssim(volume_array1, volume_array2)

    print(f"Dice Similarity Coefficient (DSC):     {dice_score:.6f}")
    print(f"Normalized Mean Square Error (NMSE):   {nmse_value:.6f}")
    print(f"Structural Similarity Index Measure:   {ssim_value:.6f}")

    # Subtract grayscale volumes and display results
    if show_subtraction:
        subtracted_volume = subtract_volumes(volume_array1, volume_array2)
        normalized_subtracted_volume = subtracted_volume / np.max(volume_array1)
        display_subtraction_results(volume_array1, volume_array2, normalized_subtracted_volume, dice_score, nmse_value, ssim_value)

    return dice_score, nmse_value, ssim_value

def plot_similarity_metrics(results, input_dir):
    """Plot Dice, NMSE, and SSIM across time and save results to CSV."""
    instances = [r["instance"] for r in results]
    dice_vals = [r["dice"] for r in results]
    nmse_vals = [r["nmse"] for r in results]
    ssim_vals = [r["ssim"] for r in results]

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_filename = os.path.join(input_dir, f"timeseries_similarity_metrics_{timestamp}.csv")
    plot_filename_prefix = f"timeseries"

    # Save CSV
    with open(csv_filename, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Instance", "Dice", "NMSE", "SSIM"])
        for r in results:
            writer.writerow([r["instance"], r["dice"], r["nmse"], r["ssim"]])
    print(f"\nSaved similarity metrics to {csv_filename}")

    # Plot Dice
    plt.figure(figsize=(8, 5))
    plt.plot(instances, dice_vals, marker='o', color='black')
    plt.title("Dice Similarity over Time")
    plt.xlabel("Instance Number")
    plt.ylabel("Dice Coefficient")
    plt.ylim(0.75, 1.00)
    plt.grid(True)
    plt.tight_layout()
    # Save figure to input directory
    plt.savefig(os.path.join(input_dir, f"{plot_filename_prefix}_dice.png"), dpi=300)
    plt.show()

    # Plot NMSE
    plt.figure(figsize=(8, 5))
    plt.plot(instances, nmse_vals, marker='o', color='red')
    plt.title("Normalized Mean Square Error over Time")
    plt.xlabel("Instance Number")
    plt.ylabel("NMSE")
    plt.ylim(0.00, 0.10)
    plt.grid(True)
    plt.tight_layout()
    # Save figure to input directory
    plt.savefig(os.path.join(input_dir, f"{plot_filename_prefix}_nmse.png"), dpi=300)
    plt.show()

    # Plot SSIM
    plt.figure(figsize=(8, 5))
    plt.plot(instances, ssim_vals, marker='o', color='green')
    plt.title("Structural Similarity Index Measure over Time")
    plt.xlabel("Instance Number")
    plt.ylabel("SSIM")
    plt.ylim(0.5, 1.00)
    plt.grid(True)
    plt.tight_layout()
    # Save figure to input directory
    plt.savefig(os.path.join(input_dir, f"{plot_filename_prefix}_ssim.png"), dpi=300)
    plt.show()


def main(input_dir):
    """Loop through directory of volumes, using the first as reference."""
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)

    volume_files = sorted([os.path.join(input_dir, f) for f in os.listdir(input_dir)
                           if f.lower().endswith((".nii", ".nii.gz", ".nrrd", ".nhdr"))])

    if len(volume_files) < 2:
        print("Error: Need at least two image volumes in the directory.")
        sys.exit(1)

    ref_file = volume_files[0]  # set first volume to be the reference volume for comparisons
    print(f"Using first volume as reference: {os.path.basename(ref_file)}")

    results = []
    for idx, file2 in enumerate(volume_files[1:], start=1):
        dice, nmse, ssim_val = compute_volume_comparison(ref_file, file2)
        results.append({"instance": idx, "dice": dice, "nmse": nmse, "ssim": ssim_val})
        print(f"Compared {os.path.basename(file2)} --> Dice={dice:.4f}, NMSE={nmse:.4f}, SSIM={ssim_val:.4f}")

    # Plot and save
    plot_similarity_metrics(results, input_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare volumes: single or directory mode")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--inputdir", help="Directory of image volumes (first volume as reference)")
    group.add_argument("--file1", help="Reference image file for single comparison")
    parser.add_argument("--file2", help="Target image file for single comparison (required if --file1 is used)")
    args = parser.parse_args()

    if args.inputdir:
        main(args.inputdir)
    elif args.file1 and args.file2:
        compute_volume_comparison(args.file1, args.file2, show_subtraction=True)


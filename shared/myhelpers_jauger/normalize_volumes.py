#!/usr/bin/env python3
"""
Title: normalize_volumes.py

Description:
Perform intensity normalization on a directory of NIfTI image volumes using SimpleITK's intensity normalization
techniques (min-max, z-score, histogram matching).

Usage:
    python normalize_volumes.py \
        --indir /path/to/input_volumes \
        --outdir /path/to/output_volumes \
        --method minmax

Example:
    python normalize_volumes.py --indir ./raw_volumes --outdir ./norm_volumes --method histmatch

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Created on: 21 October 2025
"""

import os
import sys
import argparse
import SimpleITK as sitk
import numpy as np
import cv2


def normalize_minmax(image):
    """Min-max normalization to [0,1]."""
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    min_val, max_val = stats.GetMinimum(), stats.GetMaximum()
    if max_val == min_val:
        return image * 0  # Avoid division by zero
    return (image - min_val) / (max_val - min_val)


def normalize_zscore(image):
    """Z-score normalization (mean=0, std=1)."""
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    mean_val, std_val = stats.GetMean(), stats.GetSigma()
    if std_val == 0:
        return image * 0
    return (image - mean_val) / std_val


def normalize_histmatch(image, reference):
    """Histogram matching normalization using the reference image."""
    matcher = sitk.HistogramMatchingImageFilter()
    matcher.SetNumberOfHistogramLevels(64)
    # matcher.SetNumberOfMatchPoints(50)
    matcher.ThresholdAtMeanIntensityOn()
    return matcher.Execute(image, reference)


def normalize_histeq(image, num_bins=64):
    """Perform global histogram equalization with specified number of bins."""
    arr = sitk.GetArrayFromImage(image)
    arr = arr.astype(np.float32)

    # Create a uniform reference histogram
    ref_hist = np.linspace(np.min(arr), np.max(arr), num_bins, dtype=np.float32)
    ref_img = sitk.GetImageFromArray(np.tile(ref_hist, arr.size // num_bins + 1)[:arr.size].reshape(arr.shape))
    ref_img.CopyInformation(image)

    equalized = sitk.HistogramMatching(image, ref_img, numberOfHistogramLevels=num_bins)
    return equalized


def normalize_histeq_opencv(image):
    """
    Perform histogram equalization slice-by-slice using OpenCV.
      - cv2.equalizeHist() only supports 8-bit single-channel images.
      - We first scale each slice to [0,255], convert to uint8, equalize, then rescale to [0,1].
      - Operates on each 2D slice independently (typically axial slices).
    """
    arr = sitk.GetArrayFromImage(image).astype(np.float32)
    arr_min, arr_max = np.min(arr), np.max(arr)
    if arr_max == arr_min:
        print("Warning: zero intensity range — skipping OpenCV equalization.")
        return image

    # Scale to [0,255] for OpenCV
    arr_scaled = ((arr - arr_min) / (arr_max - arr_min) * 255.0).astype(np.uint8)
    eq_arr = np.zeros_like(arr_scaled)

    # Apply OpenCV equalization per slice
    for z in range(arr_scaled.shape[0]):
        eq_arr[z] = cv2.equalizeHist(arr_scaled[z])

    # Rescale back to [0,1]
    eq_arr = eq_arr.astype(np.float32) / 255.0

    eq_img = sitk.GetImageFromArray(eq_arr)
    eq_img.CopyInformation(image)
    return eq_img

def normalize_mean_scaling(image, target_mean):
    """
    Scale the nonzero voxel mean of an image to a specified target value.
    Voxels with value 0 (e.g., background) are ignored in the scaling factor computation.
    """
    arr = sitk.GetArrayFromImage(image).astype(np.float32)
    nonzero_mask = arr != 0

    # Compute mean of nonzero voxels
    current_mean = np.mean(arr[nonzero_mask]) if np.any(nonzero_mask) else 0
    if current_mean == 0:
        print("Warning: Image has no nonzero voxels, skipping normalization.")
        return image

    scale_factor = target_mean / current_mean
    arr_scaled = arr * scale_factor
    # Keep zeros as zeros
    arr_scaled[~nonzero_mask] = 0

    normalized_image = sitk.GetImageFromArray(arr_scaled)
    normalized_image.CopyInformation(image)
    return normalized_image


def get_nifti_files(indir):
    """Return sorted list of .nii or .nii.gz files in directory."""
    files = sorted(
        [
            os.path.join(indir, f)
            for f in os.listdir(indir)
            if f.endswith(".nii") or f.endswith(".nii.gz")
        ]
    )
    if not files:
        sys.exit(f"Error: No NIfTI files found in directory: {indir}")
    return files


def strip_extension(filename):
    """Strip .nii or .nii.gz extension for consistent output naming."""
    if filename.endswith(".nii.gz"):
        return filename[:-7]
    elif filename.endswith(".nii"):
        return filename[:-4]
    else:
        return filename


def main():
    parser = argparse.ArgumentParser(description="Normalize NIfTI image volumes with SimpleITK.")
    parser.add_argument("--indir", required=True, help="Input directory containing NIfTI files")
    parser.add_argument("--outdir", help="Output directory for normalized NIfTI files (optional)")
    parser.add_argument(
        "--method",
        required=True,
        choices=["minmax", "zscore", "histmatch", "histeq", "histeq_opencv", "mean_scaling"],
        help="Normalization method: minmax | zscore | histmatch | histeq | histeq_opencv | mean_scaling",
    )
    parser.add_argument("--bins", type=int, default=64, help="Number of bins for histogram equalization")
    parser.add_argument("--targetmean", type=float, default=None, help="Target mean value for mean scaling (required if method=mean_scaling)."
    )
    args = parser.parse_args()

    indir = os.path.abspath(args.indir)
    method = args.method

    # Determine output directory
    if args.outdir:
        outdir = os.path.abspath(args.outdir)
    else:
        outdir = f"{indir.rstrip(os.sep)}_{method}_normalized"

    os.makedirs(outdir, exist_ok=True)
    print(f"Output directory set to: {outdir}")

    nifti_files = get_nifti_files(indir)

    # Load reference image if histogram matching
    reference_image = None
    if method == "histmatch":
        print(f"Using first volume as reference for histogram matching: {nifti_files[0]}")
        reference_image = sitk.ReadImage(nifti_files[0], sitk.sitkFloat32)

    # Process each volume
    for i, infile in enumerate(nifti_files):
        base_name = os.path.basename(infile)
        name_root = strip_extension(base_name)
        outfile = os.path.join(outdir, f"{name_root}_{method}_normalized.nii")

        print(f"Processing ({i + 1}/{len(nifti_files)}): {base_name}")
        image = sitk.ReadImage(infile, sitk.sitkFloat32)

        if method == "minmax":
            norm_img = normalize_minmax(image)
        elif method == "zscore":
            norm_img = normalize_zscore(image)
        elif method == "histmatch":
            norm_img = normalize_histmatch(image, reference_image)
        elif method == "histeq":
            norm_img = normalize_histeq(image, num_bins=args.bins)
        elif method == "histeq_opencv":
            norm_img = normalize_histeq_opencv(image)
        elif args.method == "mean_scaling":
            if args.targetmean is None:
                parser.error("--targetmean must be specified when method=mean_scaling")
            norm_img = normalize_mean_scaling(image, args.targetmean)
        else:
            sys.exit(f"Unknown normalization method: {method}")

        sitk.WriteImage(norm_img, outfile)
        print(f" -> Saved: {outfile}")

    print(f"\nAll volumes processed successfully. Output saved to: {outdir}")


if __name__ == "__main__":
    main()
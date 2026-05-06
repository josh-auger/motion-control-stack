"""
Title: calculate_tsnr.py

Description:
Compute the temporal signal-to-noise ratio (tSNR) of a directory of image volume files and generate the equivalent 3D SNR volume map.
Library dependencies: nibabel, matplotlib and scipy

References:
    Adapted from https://gist.github.com/arokem/937534/0dec33cc0c302292642b10bce0a5855737550b1f

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Date of creation: March 19, 2025
"""

import os
import sys
import shutil
import numpy as np
from glob import glob
from tqdm import tqdm
import SimpleITK as sitk
import matplotlib.pyplot as plt
from scipy.io import savemat
import argparse
from run_synthstrip_docker import run_synthstrip
from display_tsnr_results import display_tsnr_components, display_tsnr_slices

def usage():
    print("Usage: python3 calculate_tsnr.py [option flags] <input_directory>")


def help_me():
    help_text = """ 
    This script calculates the temporal Signal-to-Noise Ratio (tSNR) for a directory of image volume files.

    Usage:
        python3 calculate_tsnr.py [option flags] <input_directory>

    Input(s):
        input_directory   Path to the directory containing 3D image files.

    Options:
        -h, --help        Display help message and exit script.
        -v, --verbose     Verbose mode for more detailed logging during processing.
        -p, --plot        Generate and save a bar plot of tSNR values for each processed dataset.

    Output(s):
        - Generates sub-directory '/TSNR' inside input directory for all output files.
        - For each processed dataset, a tSNR image volume file is saved.
        - A MATLAB (.mat) file of the tSNR values is also saved for each dataset.
        - If --plot is specified, a bar plot of all tSNR values is saved as 'plot_tsnr.png'.
        
    Additional note(s) :
        - Image data is read in using SimpleITK GetArrayFromImage() which does NOT account for spatial metadata within
            the image header. All image data is assumed to reside on the same cartesian lattice. Confirm that all images
            maintain the same spatial metadata (image frame) or resample image data to the same image frame prior to 
            tSNR calculation.
    """
    print(help_text)
    sys.exit(0)


def load_image(file):
    """
    Read image file using SimpleITK and load image data as a matrix array.
    NOTE: this does NOT account for any spatial metadata of the image data! Each cartesian pixel in image space is read
    in as a matrix index.
    """
    sitk_img = sitk.ReadImage(file)  # SimpleITK automatically detects format
    data = sitk.GetArrayFromImage(sitk_img)  # numpy array [z,y,x]
    return data, sitk_img


def process_image_files(input_path, plot, verbose, mask):
    """Processes a directory of image volume files and calculates tSNR."""
    if not os.path.exists(input_path):
        sys.exit(f"Error: Directory '{input_path}' not found.")

    # Ingest all image files in input directory, sort alphabetically by filename
    image_files = sorted(
        glob(os.path.join(input_path, "*.nii*")) +
        glob(os.path.join(input_path, "*.nrrd")) +
        glob(os.path.join(input_path, "*.nhdr"))
    )
    if not image_files:
        sys.exit(f"Error: No image files found in directory : {input_path}")
    if verbose:
        print(f"Processing {len(image_files)} image files from {input_path}...")

    # --- Create reference ---
    first_file = image_files.pop(0)
    first_data, reference_img = load_image(first_file)

    data_list = [first_data]
    base_shape = first_data.shape

    for file in tqdm(image_files, desc="Concatenating image files", unit="file"):
        new_data, _ = load_image(file)
        if new_data.shape != base_shape:
            print(f"Skipping {file}: Shape mismatch {new_data.shape} != {base_shape}")
            continue
        data_list.append(new_data)

    combined_data = np.stack(data_list, axis=-1)
    if verbose:
        print(f"Concatenated data shape : {combined_data.shape}")

    return combined_data, reference_img, first_file


def apply_mask_to_image(image_file, mask_file, output_file):
    """ Apply binary mask (mask_file) to an image (image_file) using SimpleITK."""
    img = sitk.ReadImage(image_file)
    mask = sitk.ReadImage(mask_file)
    # Ensure mask is binary
    mask = sitk.Cast(mask > 0, sitk.sitkUInt8)
    # Convert mask to float for multiplication
    mask = sitk.Cast(mask, sitk.sitkFloat32)

    # Apply mask to image
    masked_img = sitk.Cast(img, sitk.sitkFloat32) * mask
    masked_img.CopyInformation(img)

    sitk.WriteImage(masked_img, output_file)
    return masked_img


def calculate_tsnr(data, verbose=False):
    """Compute temporal SNR (tSNR) from an image volume dataset and save outputs with SimpleITK."""
    # last axis dimension should be time
    mean_data = np.mean(data, axis=-1)
    std_data = np.std(data, axis=-1)

    # Use np.divide to safely handle division by zero
    tsnr_map = np.divide(
        mean_data,
        std_data,
        out=np.zeros_like(mean_data, dtype=np.float32),  # default where std=0
        where=std_data != 0
    )

    return tsnr_map, mean_data, std_data


def summarize_tsnr(tsnr_map, label="TSNR Summary"):
    tsnr_nonzero = tsnr_map[tsnr_map > 0]
    mean_tsnr = np.mean(tsnr_nonzero)
    std_tsnr = np.std(tsnr_nonzero)
    median_tsnr = np.median(tsnr_nonzero)
    num_high_tsnr = np.sum(tsnr_nonzero > 50)
    perc_high_tsnr = (num_high_tsnr / tsnr_nonzero.size) * 100

    print(f"\n======== {label} ========")
    print(f"Mean volume tSNR: {mean_tsnr:.4f}")
    print(f"Std-dev volume tSNR: {std_tsnr:.4f}")
    print(f"Median volume tSNR: {median_tsnr:.4f}")
    print(f"Voxels with tSNR > 50: {num_high_tsnr:,} ({perc_high_tsnr:.2f}%)")


def save_tsnr_results(reference_img, tsnr_map, mean_data, std_data, output_path, verbose=False):
    # Convert numpy arrays to SimpleITK images
    tsnr_img = sitk.GetImageFromArray(tsnr_map.astype(np.float32))
    mean_img = sitk.GetImageFromArray(mean_data.astype(np.float32))
    std_img = sitk.GetImageFromArray(std_data.astype(np.float32))

    # Copy spacing, origin, direction from reference image
    tsnr_img.CopyInformation(reference_img)
    mean_img.CopyInformation(reference_img)
    std_img.CopyInformation(reference_img)

    # Save as NIfTI images
    sitk.WriteImage(tsnr_img, os.path.join(output_path, "tsnr_map.nii"))
    sitk.WriteImage(mean_img, os.path.join(output_path, "mean_signal.nii"))
    sitk.WriteImage(std_img, os.path.join(output_path, "stddev_signal.nii"))

    # Save as .mat files
    savemat(os.path.join(output_path, "tsnr_values.mat"), {'tsnr': tsnr_map})
    savemat(os.path.join(output_path, "mean_signal_values.mat"), {'mean': mean_data})
    savemat(os.path.join(output_path, "stddev_signal_values.mat"), {'std': std_data})

    if verbose:
        print(f"\ntSNR image saved: {output_path}_tsnr.nii.gz")
        print(f"Mean signal image saved: {output_path}_mean_signal.nii.gz")
        print(f"Standard deviation image saved: {output_path}_stddev_signal.nii.gz")
        print(f"tSNR values saved: {output_path}_tsnr.mat")
        print(f"Mean signal values saved: {output_path}_mean_signal.mat")
        print(f"Standard deviation values saved: {output_path}_std_signal.mat")


def main():
    parser = argparse.ArgumentParser(add_help=False)  # Disable default help
    parser.add_argument("input_directory", nargs="?", help="Path to the directory containing image files.")
    parser.add_argument("-p", "--plot", action="store_true", help="Generate and save a plot of mean tSNR values of each dataset.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging mode.")
    parser.add_argument("-m", "--mask", action="store_true", help="Call FSL SynthStrip to mask the brain.")
    parser.add_argument("-h", "--help", action="store_true", help="Show help message and exit.")
    args = parser.parse_args()

    if args.help:
        help_me()

    if not args.input_directory:
        print("Error: No input directory provided.\n")
        help_me()

    # Create /TSNR sub-directory for output files
    tsnr_path = os.path.join(args.input_directory, "TSNR")
    os.makedirs(tsnr_path, exist_ok=True)

    # Read in all image volumes and concatenate image data
    combined_data, reference_img, first_file = process_image_files(args.input_directory, args.plot, args.verbose, args.mask)

    # Calculate TSNR from image data
    tsnr_map, mean_data, std_data = calculate_tsnr(combined_data, args.verbose)
    save_tsnr_results(reference_img, tsnr_map, mean_data, std_data, tsnr_path, args.verbose)

    # Mask the brain for all tsnr results using skull-stripped reference volume
    if args.mask:
        if args.verbose:
            print(f"\nStripping skull from tSNR results using SynthStrip...")
        mean_file = os.path.join(tsnr_path, "mean_signal.nii")
        std_file = os.path.join(tsnr_path, "stddev_signal.nii")
        tsnr_file = os.path.join(tsnr_path, "tsnr_map.nii")

        # Call SynthStrip on first volume reference to get brain mask
        run_synthstrip(first_file, output_mask=True, erode_mask=True)

        # Extract the base filename (handle .nii.gz safely)
        first_base = os.path.basename(first_file)
        if first_base.endswith('.nii.gz'):
            first_base = first_base[:-7]  # remove '.nii.gz'
        else:
            first_base = os.path.splitext(first_base)[0]


        strip_file = os.path.join(args.input_directory, f"{first_base}_strip.nii")
        mask_file = os.path.join(args.input_directory, f"{first_base}_strip_mask.nii")
        # Define new destinations inside /TSNR
        strip_dest = os.path.join(tsnr_path, os.path.basename(strip_file))
        mask_dest = os.path.join(tsnr_path, os.path.basename(mask_file))
        # Move the SynthStrip outputs to /TSNR
        if os.path.exists(strip_file):
            shutil.move(strip_file, strip_dest)
            if args.verbose:
                print(f"Moved stripped image to: {strip_dest}")
        if os.path.exists(mask_file):
            shutil.move(mask_file, mask_dest)
            if args.verbose:
                print(f"Moved mask file to: {mask_dest}")

        # Use the relocated mask file from /TSNR for applying to tSNR outputs
        mask_file = mask_dest

        if args.verbose:
            print(f"Applying reference volume mask to mean signal, stddev signal, and tsnr map...")
        # Apply the same mask to stddev and tsnr maps
        mean_masked = apply_mask_to_image(mean_file, mask_file, os.path.join(tsnr_path, "mean_signal_strip.nii"))
        std_masked = apply_mask_to_image(std_file, mask_file, os.path.join(tsnr_path, "stddev_signal_strip.nii"))
        tsnr_masked = apply_mask_to_image(tsnr_file, mask_file, os.path.join(tsnr_path, "tsnr_map_strip.nii"))

        # Load masked data arrays for summary/plot
        tsnr_map_stripped = sitk.GetArrayFromImage(tsnr_masked)
        mean_data_stripped = sitk.GetArrayFromImage(mean_masked)
        std_data_stripped = sitk.GetArrayFromImage(std_masked)

        if args.verbose:
            summarize_tsnr(tsnr_map_stripped, "(Masked) TSNR Summary")

        if args.plot:
            display_tsnr_components(reference_img, mean_data_stripped, std_data_stripped, tsnr_map_stripped, tsnr_path)
    else:
        if args.verbose:
            summarize_tsnr(tsnr_map, "TSNR Summary")
        if args.plot:
            display_tsnr_components(reference_img, mean_data, std_data, tsnr_map, tsnr_path)
            display_tsnr_slices(tsnr_map, tsnr_path)


if __name__ == "__main__":
    main()
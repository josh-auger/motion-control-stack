#!/usr/bin/env python3
"""
Title: convert_3Dnrrd_to_4Dnifti.py

Description:
Convert a directory of 3D NRRD (.nhdr) image volumes into a single 4D NIfTI (.nii.gz) file for further processing.

Usage:
    python convert_3Dnrrd_to_4Dnifti.py --inputdir /path/to/nrrd_dir --outfile /path/to/output/merged_func.nii.gz

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Date of creation: October 17, 2025
"""

import os
import sys
import argparse
import glob
import SimpleITK as sitk
import numpy as np
import datetime

def check_input_directory(indir):
    """Check that directory exists and contains .nhdr files."""
    if not os.path.isdir(indir):
        sys.exit(f"Error : Input directory '{indir}' not found.")
    nrrd_files = sorted(glob.glob(os.path.join(indir, "*.nhdr")))
    if len(nrrd_files) == 0:
        sys.exit(f"Error : No *.NHDR files found in {indir}")
    return nrrd_files

def load_nrrd_volumes(file_list, verbose=False):
    """Read each NRRD 3D volume and return as numpy array list."""
    volumes = []
    for idx, fpath in enumerate(file_list):
        img = sitk.ReadImage(fpath)
        arr = sitk.GetArrayFromImage(img)  # shape: (z, y, x)
        volumes.append(arr)
        if verbose:
            print(f"Loaded {fpath} with shape {arr.shape}")
    return volumes, img  # Return last image for metadata (affine, spacing, etc.)

def stack_to_4d(volumes):
    """Stack list of 3D volumes into a 4D numpy array."""
    vol4d = np.stack(volumes, axis=0)  # (time, z, y, x)
    vol4d = np.transpose(vol4d, (1, 2, 3, 0))  # -> (z, y, x, t)
    return vol4d

def save_nifti(vol4d, ref_img, outfile):
    """Save 4D NIfTI image."""
    nii_img = sitk.GetImageFromArray(vol4d)  # SimpleITK expects (t,z,y,x) unless transposed
    nii_img.CopyInformation(ref_img)
    sitk.WriteImage(nii_img, outfile)
    print(f"Saved merged 4D NIfTI: {outfile}")

def main():
    parser = argparse.ArgumentParser(description="Convert a directory of 3D NRRD volumes into a single 4D NIfTI file.")
    parser.add_argument("--inputdir", required=True, help="Input directory containing .nhdr files")
    parser.add_argument("--outfile", help="Output .nii.gz file path")
    parser.add_argument("--verbose", action="store_true", help="Print progress messages")
    args = parser.parse_args()

    nrrd_files = check_input_directory(args.inputdir)

    # Default output path → parent of input directory
    if not args.outfile:
        parent_dir = os.path.dirname(os.path.abspath(args.inputdir.rstrip("/")))
        timestamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
        base_name = os.path.basename(os.path.normpath(args.inputdir))
        args.outfile = os.path.join(parent_dir, f"{base_name}_merged_func_{timestamp}.nii.gz")

    vols, ref_img = load_nrrd_volumes(nrrd_files, verbose=args.verbose)
    vol4d = stack_to_4d(vols)
    save_nifti(vol4d, ref_img, args.outfile)

    print(f"4D NIfTI shape: {vol4d.shape}")
    print("Conversion complete.")

if __name__ == "__main__":
    main()

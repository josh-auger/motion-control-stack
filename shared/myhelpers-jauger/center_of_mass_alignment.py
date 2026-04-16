#!/usr/bin/env python3

# Title: center_of_mass_alignment.py

# Description:
# Calculate intensity-weighted center of mass of two image volumes and return the translation parameters to align the
# reference volume center of mass with the target volume center of mass.
#
# Created on: March 2026
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import SimpleITK as sitk
import numpy as np
import argparse
import sys
import os


def read_image(image_path):
    """Read image using SimpleITK."""
    if not os.path.exists(image_path):
        raise FileNotFoundError(f"File not found: {image_path}")
    return sitk.ReadImage(image_path)


def compute_center_of_mass(image, normalize=False):
    """
    Fast, vectorized center of mass in physical space.
    Works for 3D images (including thin-slice volumes).
    """
    array = sitk.GetArrayFromImage(image).astype(np.float64)

    if normalize:
        total = np.sum(array)
        if total > 0:
            array /= total

    spacing = np.array(image.GetSpacing())      # (x, y, z)
    origin = np.array(image.GetOrigin())        # (x, y, z)
    direction = np.array(image.GetDirection()).reshape(3, 3)

    # Create index grid (z, y, x)
    z, y, x = np.indices(array.shape)
    # Convert to (x, y, z)
    indices = np.stack([x, y, z], axis=-1)
    # Scale to physical mm units
    scaled = indices * spacing
    # Apply direction matrix
    physical = np.tensordot(scaled, direction.T, axes=1)
    # Add origin
    physical += origin

    weights = array
    total_mass = np.sum(weights)
    if total_mass == 0:
        raise ValueError("Image has zero total intensity.")

    com = np.sum(physical * weights[..., None], axis=(0, 1, 2)) / total_mass
    return com, total_mass


def compute_combined_com(images, normalize=False):
    """
    Compute global COM from multiple 3D images (e.g., slices).
    """
    weighted_sum = np.zeros(3)
    total_mass = 0.0

    for i, img in enumerate(images):
        com, mass = compute_center_of_mass(img, normalize=normalize)
        print(f"Slice {i}: COM = {com}, Mass = {mass}")
        weighted_sum += com * mass
        total_mass += mass

    if total_mass == 0:
        raise ValueError("Total mass across all slices is zero.")

    combined_com = weighted_sum / total_mass
    return combined_com


def compute_translation(ref_image, target_images, normalize=False):
    """
    Compute translation vector to move ref -> target.
    target_images: list of one or more SimpleITK images
    """
    com_ref, _ = compute_center_of_mass(ref_image, normalize=normalize)

    if len(target_images) == 1:
        com_target, _ = compute_center_of_mass(target_images[0], normalize=normalize)
    else:
        com_target = compute_combined_com(target_images, normalize=normalize)

    translation = com_target - com_ref
    return translation, com_ref, com_target


def main(args):
    ref_img = read_image(args.refvolume)
    target_imgs = [read_image(p) for p in args.slicetarget]
    print("\nLoaded target images:")
    for i, p in enumerate(args.slicetarget):
        print(f"  [{i}] {p}")

    translation, com_ref, com_target = compute_translation(ref_img, target_imgs, normalize=args.normalize)
    print("\nCenter of Mass (Reference):", com_ref)
    print("Center of Mass (Target):   ", com_target)
    print("\nTranslation (X, Y, Z):     ", translation)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute center-of-mass translation between reference and target(s).")
    parser.add_argument("--refvolume", required=True, help="Path to reference image volume")
    parser.add_argument("--slicetarget",required=True,nargs='+',help="Path to 1+ target image files (space-separated)")
    parser.add_argument("--normalize",action="store_true",help="Normalize images by total intensity before COM calculation")
    args = parser.parse_args()
    main(args)
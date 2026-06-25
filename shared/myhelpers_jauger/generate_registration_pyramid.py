# Title: generate_registration_pyramid.py

# Description:
# Build a multi-resolution pyramid of the input image.
# Each level first applies a 1-voxel low-pass filter to smooth the image data and prevent aliasing. Then, the image is
# downsampled by a specified factor for each level. Process is repeated for N levels.
#

# Created on: March 2026
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import os
import argparse
import SimpleITK as sitk


def compute_new_size(original_size, original_spacing, new_spacing):
    """Compute image size required to maintain the same physical extent."""
    return [
        int(round(osz * ospc / nspc))
        for osz, ospc, nspc in zip(original_size, original_spacing, new_spacing)
    ]

def gaussian_smooth(image, sigma):
    """Apply Gaussian smoothing with sigma in physical units."""
    return sitk.SmoothingRecursiveGaussian(image, sigma)

def resample_to_spacing(image, new_spacing):
    """Resample image to new spacing."""
    original_spacing = image.GetSpacing()
    original_size = image.GetSize()
    new_size = compute_new_size(original_size, original_spacing, new_spacing)

    resampled = sitk.Resample(
        image,
        new_size,
        sitk.Transform(),
        sitk.sitkLinear,
        image.GetOrigin(),
        new_spacing,
        image.GetDirection(),
        0,
        image.GetPixelID(),
    )
    return resampled

def build_registration_pyramid(image, levels=3, shrink_factor=2.0, include_upsample=True):
    """
    Build a multi-resolution image processing pyramid.
    Each level is resampled from the original image to avoid cumulative interpolation artifacts.
    """
    pyramid = []
    original_spacing = image.GetSpacing()
    print(f"Original spacing = {original_spacing}")

    # Optionally add one upsample level
    level_indices = list(range(levels))
    if include_upsample:
        level_indices = [-1] + level_indices

    for level in level_indices:
        factor = shrink_factor ** level
        new_spacing = [s * factor for s in original_spacing]

        if factor < 1.0:
            # Upsampling: no smoothing
            sigma = 0.0
            smoothed = image
        elif level == 0:
            # Original image (no smoothing, no resampling)
            pyramid.append(image)
            print(f"Level {level}")
            print(f"\tFactor = {factor}")
            print(f"\tSigma = 0.0 (original image)")
            print(f"\tSpacing = {original_spacing}")
            continue
        else:
            # Downsampling: apply anti-alias smoothing
            sigma = new_spacing[0] / 2.0
            smoothed = gaussian_smooth(image, sigma)

        level_image = resample_to_spacing(smoothed, new_spacing)
        pyramid.append(level_image)

        print(f"Level {level}")
        print(f"\tFactor = {factor}")
        print(f"\tSigma = {sigma}")
        print(f"\tNew spacing = {new_spacing}")

    return pyramid


def main():
    parser = argparse.ArgumentParser(description="Build a registration-grade multi-resolution pyramid.")
    parser.add_argument("--inputimage", required=True, help="Input image.")
    parser.add_argument("--levels", type=int, default=3, help="Number of pyramid levels.")
    parser.add_argument("--factor", type=float, default=2.0, help="Shrink factor per level.")
    parser.add_argument("--outputprefix", help="Prefix for output files.")
    args = parser.parse_args()

    image = sitk.ReadImage(args.inputimage)

    pyramid = build_registration_pyramid(image, levels=args.levels, shrink_factor=args.factor, include_upsample=True)

    if args.outputprefix:
        prefix = args.outputprefix
    else:
        base = os.path.basename(args.inputimage)
        prefix = os.path.splitext(base)[0]

    for i, level_img in enumerate(pyramid):
        filename = f"{prefix}_pyramidLevel{i}_2.nii"
        sitk.WriteImage(level_img, filename)
        print(f"Saved level {i} as : {filename}")


if __name__ == "__main__":
    main()
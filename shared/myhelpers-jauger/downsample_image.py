# Title: downsample_image.py

# Description:
# Downsample the input image along all axes by a factor X.
#

# Created on: March 2026
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import os
import argparse
import SimpleITK as sitk


def downsample_image(itk_image, factor=2):
    """
    Downsample an image by a specified factor in all axes.
    """
    original_spacing = itk_image.GetSpacing()
    original_size = itk_image.GetSize()

    new_spacing = [s * factor for s in original_spacing]
    new_size = [int(sz / factor) for sz in original_size]

    resampled_image = sitk.Resample(
        itk_image,
        new_size,
        sitk.Transform(),
        sitk.sitkLinear,
        itk_image.GetOrigin(),
        new_spacing,
        itk_image.GetDirection(),
        0,
        itk_image.GetPixelID(),
    )

    return resampled_image


def main():
    parser = argparse.ArgumentParser(description="Downsample an ITK image.")
    parser.add_argument("--inputimage", required=True, help="Path to the input image file.")
    parser.add_argument("--outputimage", required=False, help="Filename for the output image.")
    parser.add_argument("--factor", type=float, default=2.0, help="Downsampling factor applied to all axes (default: 2.0).")
    args = parser.parse_args()

    input_image = sitk.ReadImage(args.inputimage)

    # Apply downsampling
    downsampled_image = downsample_image(input_image, factor=args.factor)

    # Determine output filename
    if args.outputimage:
        output_filename = args.outputimage
    else:
        input_dir = os.path.dirname(args.inputimage)
        base_name = os.path.basename(args.inputimage)
        name_without_ext = os.path.splitext(base_name)[0]
        output_filename = os.path.join(input_dir, f"{name_without_ext}_downsampled.nii")

    sitk.WriteImage(downsampled_image, output_filename)

    print(f"Downsampled image saved as: {output_filename}")


if __name__ == "__main__":
    main()
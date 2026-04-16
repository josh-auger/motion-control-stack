# Title: low_pass_filter_image.py

# Description:
# Applies a low-pass filter to an ITK image using SimpleITK
#

# Created on: January 2025
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import os
import SimpleITK as sitk
import argparse


def apply_low_pass_filter(itk_image, sigma=1.0):
    # Create the Gaussian filter
    gaussian_filter = sitk.RecursiveGaussianImageFilter()
    gaussian_filter.SetSigma(sigma)
    # Apply the filter along each dimension
    smoothed_image = gaussian_filter.Execute(itk_image)
    return smoothed_image


def main():
    parser = argparse.ArgumentParser(description="Apply a low-pass filter to an ITK image.")
    parser.add_argument("--inputimage", required=True, help="Path to the input image file.")
    parser.add_argument("--outputimage", required=False, help="Filename for the output filtered image.")
    parser.add_argument("--sigma", type=float, default=1.0, help="Standard deviation for the Gaussian kernel (default: 1.0).")
    args = parser.parse_args()

    input_image = sitk.ReadImage(args.inputimage)

    # Apply the low-pass filter
    smoothed_image = apply_low_pass_filter(input_image, sigma=args.sigma)

    # Determine output filename
    if args.outputimage:
        output_filename = args.outputimage
    else:
        # Use input filename with "_filtered"
        base_name = os.path.basename(args.inputimage)
        name_without_ext = os.path.splitext(base_name)[0]
        output_filename = f"{name_without_ext}_filtered.nii"

    # Write the filtered image
    sitk.WriteImage(smoothed_image, output_filename)
    print(f"Filtered image saved as: {output_filename}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def read_image(image_path):
    """Read an image volume using SimpleITK."""
    if not os.path.exists(image_path):
        print(f"Error: File '{image_path}' not found.")
        sys.exit(1)

    try:
        image = sitk.ReadImage(image_path)
        print(f"Loaded image: {image_path}")
        print(f"Size: {image.GetSize()}")
        print(f"Pixel type: {image.GetPixelIDTypeAsString()}")
        return image
    except Exception as e:
        print(f"Error reading image file: {e}")
        sys.exit(1)

def plot_histogram(image, bins=100, output=None):
    """Generate and display/save histogram of pixel values."""
    array = sitk.GetArrayFromImage(image).flatten()
    print(f"Image array shape: {array.shape}")
    print(f"Value range: {np.min(array)} to {np.max(array)}")

    plt.figure(figsize=(8, 5))
    plt.hist(array, bins=bins, color='steelblue', edgecolor='black', alpha=0.7)
    plt.title("Image Histogram")
    plt.xlabel("Pixel Value")
    plt.ylabel("Frequency")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(0, 100)
    plt.ylim(0, 2000)

    if output is not None:
        plt.savefig(output)

    if output:
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"Histogram saved to {output}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot histogram of image pixel values using SimpleITK.")
    parser.add_argument("image_path", help="Path to image volume file (e.g., .nii, .nrrd, .mha, .dcm, etc.)")
    parser.add_argument("--bins", type=int, default=100, help="Number of histogram bins (default: 100)")
    parser.add_argument("--output", help="Optional path to save histogram plot (e.g., histogram.png)")
    args = parser.parse_args()

    image = read_image(args.image_path)
    plot_histogram(image, bins=args.bins, output=args.output)

if __name__ == "__main__":
    main()

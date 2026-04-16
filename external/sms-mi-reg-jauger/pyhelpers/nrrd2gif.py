import argparse
import SimpleITK as sitk
import imageio
import numpy as np


def read_image(file_path, slice_number=None):
    image = sitk.ReadImage(file_path)
    array = sitk.GetArrayFromImage(image)
    # If the image is 3D, take the specified slice or the middle slice by default
    if array.ndim == 3:
        if slice_number is None:
            slice_number = array.shape[0] // 2  # middle slice
        if slice_number < 0 or slice_number >= array.shape[0]:
            raise ValueError(f"Slice number {slice_number} is out of range for a 3D image with {array.shape[0]} slices")
        array = array[slice_number]
    return array


def normalize_image(image):
    # Normalize the image to the range [0, 255]
    image = image - np.min(image)
    image = image / np.max(image) * 255
    return image.astype(np.uint8)


def create_gif(image1, image2, output_file, duration=0.5):
    images = [image1, image2]
    imageio.mimsave(output_file, images, duration=duration)


def main(image1_path, image2_path, output_file, duration, slice_number):
    image1 = read_image(image1_path, slice_number)
    image2 = read_image(image2_path, slice_number)

    # Normalize images
    image1 = normalize_image(image1)
    image2 = normalize_image(image2)

    # Ensure the images have the same dimensions
    if image1.shape != image2.shape:
        raise ValueError("The input images must have the same dimensions")

    # Create the GIF
    create_gif(image1, image2, output_file, duration)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a GIF toggling between two medical images.")
    parser.add_argument('--imagepath1', type=str, required=True, help='Path to the first image file (.nrrd or .nii)')
    parser.add_argument('--imagepath2', type=str, required=True, help='Path to the second image file (.nrrd or .nii)')
    parser.add_argument('--outputfilename', type=str, required=True, help='Path to save the output GIF')
    parser.add_argument('--duration', type=float, default=0.5,
                        help='Duration for each frame in the GIF (default: 0.5 seconds)')
    parser.add_argument('--slice', type=int, default=None,
                        help='Slice number to extract from 3D images (default: middle slice)')
    args = parser.parse_args()

    main(args.imagepath1, args.imagepath2, args.outputfilename, args.duration, args.slice)


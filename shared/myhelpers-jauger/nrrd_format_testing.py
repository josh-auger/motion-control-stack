"""
Title: nrrd_format_testing.py

Description:
    Test script to understand NRRD format for saving byte-strings of image data to local disk.

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Created on: 2 April 2025

"""

import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt


def generate_image_data(width, height, depth):
    """
    Numpy array format follows a row-column (i,j) convention for matrices, where rows = height and columns = width
    """
    print(f"Width : {width}")   # x
    print(f"Height : {height}") # y
    print(f"Depth : {depth}")   # z
    volume_array = np.zeros((height, width, depth), dtype=np.uint8) # Pre-allocate volume array to hold slices, row-column (y x z) convention!
    slices = []
    for i in range(depth):
        img_array = np.random.randint(0, 256, size=(height, width), dtype=np.uint8) # 2D image array is row-column (y x) convention
        slices.append(img_array)
        volume_array[:, :, i] = img_array   # 3D matrix of all 2D image slices
    print(f"First image slice array shape : {slices[0].shape}")
    print(f"Volume array shape : {volume_array.shape}")  # shape = row-column convention! (y x z)
    return volume_array, slices


def save_as_raw_bytes(slices, filename_prefix):
    for i, img in enumerate(slices):
        with open(f"{filename_prefix}_{i}.raw", "wb") as f:
            f.write(img.tobytes())  # written as contiguous array of bytes in row-major order (x y)


def write_nhdr(width, height, depth, header_filename, filename_prefix):
    """
    NRRD format: https://teem.sourceforge.net/nrrd/format.html
    For NRRD format, the "sizes" field expects the axis order to be from fastest to slowest, where "fastest" is the axis
    associated with the coordinate which increments fastest as the samples are traversed in linear order in memory. NRRD
    does not assume any names for axes. The "space directions" field needs to match this same axis order.
        So, for an image with dimensions 256 x 256 x 5 where voxels are ordered in memory as x, then y, then z, the "sizes"
        field should be listed as "256 256 5".
    HOWEVER, when using the LIST feature for data, NRRD expects each file corresponds to one slice of the fastest changing
    dimension, that is, the first listed sizes axis!
    """
    dir_x = [1, 0, 0]
    dir_y = [0, 1, 0]
    dir_z = [0, 0, 1]

    print(f"dir_x : {dir_x}")
    print(f"dir_y : {dir_y}")
    print(f"dir_z : {dir_z}")

    with open(header_filename, "w") as hdr:
        hdr.write("NRRD0004\n")
        hdr.write("# Created by manual header writer\n")
        hdr.write("type: uint8\n")
        hdr.write("dimension: 3\n")
        hdr.write(f"sizes: {width} {height} {depth}\n") # LPS shape (x y z)
        hdr.write("space: left-posterior-superior\n")   # specify LPS space coordinate system
        hdr.write(f"space origin: (0.0, 0.0, 0.0)\n")
        hdr.write(f"space directions: ({dir_x[0]},{dir_x[1]},{dir_x[2]}) ({dir_y[0]},{dir_y[1]},{dir_y[2]}) ({dir_z[0]},{dir_z[1]},{dir_z[2]})\n")   # vector order reflects size order (x y z)
        hdr.write("encoding: raw\n")
        hdr.write(f"data file: LIST\n")
        for i in range(depth):
            hdr.write(f"{filename_prefix}_{i}.raw\n")


def print_nhdr(header_filename):
    with open(header_filename, "r") as hdr:
        header_contents = hdr.read()
        print("\n--- NRRD Header Contents ---")
        print(header_contents)
        print("--- End of Header ---\n")


def load_nrrd_data(header_filename):
    """
    SimpleITK expects LIST files to run along the last sizes axis.
    ReadImage() will interpret an ITK image in (x y z) physical space.
    HOWEVER, when reading an image into a Numpy array with GetArrayFromImage(), it will load the array in "C-style", with dimensions ordered (z y x)!
    This is because the Numpy memory layout convention is from slowest to fastest.
    In this way, SimpleITK expects NHDR sizes field as (z y x) order.
    """
    image_raw = sitk.ReadImage(header_filename)    # loaded as shape (x y z)
    print(f"SimpleITK image object shape : {image_raw.GetSize()}")
    # Convert back to NumPy array for comparison
    np_image_raw = sitk.GetArrayFromImage(image_raw) # array formed as shape (z y x)!
    print("Numpy converted array shape : ", np_image_raw.shape)
    # np_image_raw = np.transpose(np_image_raw, (2, 1, 0))
    np_image_raw = np.transpose(np_image_raw, (1, 2, 0))    # (z y x) back to (y x z) for numpy comparison
    print("Transposed numpy array shape:", np_image_raw.shape)
    print("Raw dtype : ", np_image_raw.dtype)
    print("Origin : ", image_raw.GetOrigin())
    print("Direction : ", image_raw.GetDirection())
    print("Spacing : ", image_raw.GetSpacing())
    return np_image_raw


def compare_arrays(array1, array2):
    print("Arrays equal?", np.array_equal(array1, array2))
    diff = np.abs(array1.astype(np.int16) - array2.astype(np.int16))
    print("Max pixel difference:", np.max(diff))
    print("Nonzero differences:", np.count_nonzero(diff))
    return diff


def plot_mid_slice(nrrd_volume, original_volume, diff, depth):
    middle_slice = depth // 2
    fig, axes = plt.subplots(3, 1, figsize=(5, 15))
    axes[0].imshow(nrrd_volume[:,:,middle_slice], cmap='gray')
    axes[0].set_title("Byte string NRRD data")
    axes[0].axis("off")

    axes[1].imshow(original_volume[:,:,middle_slice], cmap='gray')
    axes[1].set_title("Initial data array")
    axes[1].axis("off")

    axes[2].imshow(diff[:,:,middle_slice], cmap='hot')
    axes[2].set_title("Pixel-wise difference")
    axes[2].axis("off")

    plt.tight_layout()
    plt.show()


def main():
    width, height, depth = 128, 256, 5
    header_filename = "test_data/nrrd_format_testing/volume_header.nhdr"
    filename_prefix = "volume_bytes_slice"
    volume_array, slices = generate_image_data(width, height, depth)
    save_as_raw_bytes(slices, filename_prefix)
    write_nhdr(width, height, depth, header_filename, filename_prefix)
    print_nhdr(header_filename)
    np_image_raw = load_nrrd_data(header_filename)
    diff = compare_arrays(np_image_raw, volume_array)
    plot_mid_slice(np_image_raw, volume_array, diff, depth)

if __name__ == "__main__":
    main()

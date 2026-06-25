# Title: nifti2mrd.py

# Description:
# Convert a directory of 3D NIFTI volumes into a single MRD file containing all image data in order to replay the scan
# data stream offline.
#
# Example use:
# python3 nifti_to_mrd.py --inputdir /path/to/nifti_files --outputfilename output.mrd

# Created on: August 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital


import os
import sys
import h5py
import numpy as np
import SimpleITK as sitk
import ismrmrd
import argparse
import struct
import ctypes
from datetime import datetime


def check_dimensions(volume_data):
    """Check if the input NIFTI or NRRD file has 3 dimensions."""
    if len(volume_data.shape) != 3:
        raise ValueError(f"The input file {input_file} does not have 3 dimensions. It has dimensions: {volume_data.shape}")


def process_image_file(input_file):
    """Process a single NIFTI or NRRD file and return its slices."""
    sitk_image = sitk.ReadImage(input_file)
    volume_data = sitk.GetArrayFromImage(sitk_image)  # Shape will be [height, width, slices]
    check_dimensions(volume_data)  # Ensure the file is 3D

    return volume_data


def gather_input_files(input_dir):
    """Check if the directory exists and contains NIFTI or NRRD files."""
    if not os.path.isdir(input_dir):
        print(f"Error: The directory {input_dir} does not exist.")
        sys.exit(1)

    valid_extensions = ('.nii', '.nii.gz', '.nrrd')
    input_files = sorted([f for f in os.listdir(input_dir) if f.lower().endswith(valid_extensions)])

    if not input_files:
        print("Error: No NIFTI or NRRD files found in the specified directory.")
        sys.exit(1)

    print(f"Found {len(input_files)} image volume file(s) in the directory.")
    return input_files


def create_config_file():
    config_content = "dummy config content"
    config_bytes = config_content.encode('utf-8')

    return config_bytes[:1024].ljust(1024, b'\x00')


def create_mrd_header(first_volume):
    """Create MRD XML header based on the first volume."""
    xml_header = """<?xml version="1.0" encoding="UTF-8"?>
    <ismrmrdHeader>
        <version>1.0</version>
        <subjectInformation>
            <patientName>Doe^John</patientName>
        </subjectInformation>
        <acquisitionSystemInformation>
            <systemVendor>VendorName</systemVendor>
            <systemModel>ModelName</systemModel>
            <systemFieldStrength_T>3.0</systemFieldStrength_T>
        </acquisitionSystemInformation>
    </ismrmrdHeader>"""

    # print("MRD Header:\n", xml_header)
    xml_header_bytes = bytes(xml_header.encode('utf-8'))

    return xml_header_bytes


def create_image_header(slice_counter, matrix_size):
    """Create and serialize an ImageHeader."""
    header = ismrmrd.ImageHeader()  # initializes a 198-byte image header
    header.version = 1
    header.channels = 1
    header.matrix_size = matrix_size
    header.data_type = ismrmrd.DATATYPE_FLOAT
    header.image_index = slice_counter
    header.image_series_index = 0  # or another relevant index

    # print(f"Version: {header.version}")
    # print(f"Channels: {header.channels}")
    print(f"Matrix Size Tuple: {header.matrix_size[0]} , {header.matrix_size[1]} , {header.matrix_size[2]}")
    print(f"Data Type: {header.data_type}")
    print(f"Converted dtype: {ismrmrd.get_dtype_from_data_type(header.data_type)}")
    print(f"Image Index: {header.image_index}")
    # print(f"Image Series Index: {header.image_series_index}")

    header_bytes = bytes(header)

    # Does each header property need to be hard-coded into bytes and concatenated in the correct item order?
    # https://ismrmrd.readthedocs.io/en/latest/mrd_image_data.html
    # Something like this: header_bytes = bytes(header.version, unit16) + bytes(header.data_type) ... + bytes(header.measurement_uid, uint32) ...
    print(f"Header Bytes Length: {len(header_bytes)}")

    return header_bytes


def create_attributes(slice_counter):
    attributes = f"<imageAttributes><sliceCounter>{slice_counter}</sliceCounter></imageAttributes>"
    # attribute_bytes = attributes.encode('utf-8')
    attribute_bytes = bytes(attributes.encode('utf-8'))

    # ************ WORK-IN-PROGRESS REPORT ************
    # Tunc's offline FIRE save code is not properly reading the attributes string
    # The MRD file is correctly being written, but may still be missing some specific
    # Abandoned MRD conversion work in favor of directly feeding NIFTI image files into the registration code

    return attribute_bytes


def create_mrd_file(image_files, input_dir, output_file):
    """Create an MRD file with all slices from all NIFTI/NRRD volumes in the input directory."""
    slice_counter = 0  # Counter for image slices in the MRD file

    with h5py.File(output_file, 'w') as mrd_file:
        # Step 1: Add a config file (dummy content for example)
        config_file_bytes = create_config_file()
        mrd_file.create_dataset('Config File', data=bytearray(config_file_bytes))

        # Step 2: Create and save the MRD header
        first_volume = process_image_file(os.path.join(input_dir, image_files[0]))
        matrix_size = first_volume.shape[:2] + (1,)

        xml_header_bytes = create_mrd_header(first_volume)
        mrd_file.create_dataset('Metadata XML', data=bytearray(xml_header_bytes))

        # Process each image file in the list
        for filename in image_files:
            file_path = os.path.join(input_dir, filename)
            volume_data = process_image_file(file_path)
            num_slices = volume_data.shape[2]
            print(f"{filename} contains {num_slices} image slices.")

            # Add each slice to the MRD file
            for slice_idx in range(num_slices):
                # Extract the individual slice (2D image)
                slice_data = volume_data[:, :, slice_idx]
                header_bytes = create_image_header(slice_counter, matrix_size)
                attribute_bytes = create_attributes(slice_counter)
                slice_data_bytes = slice_data.astype(np.float32).tobytes()

                # Directly create datasets for each image slice
                mrd_file.create_dataset(f"image_{slice_counter}/header", data=bytearray(header_bytes))
                mrd_file.create_dataset(f"image_{slice_counter}/attribute", data=bytearray(attribute_bytes))
                mrd_file.create_dataset(f"image_{slice_counter}/data", data=bytearray(slice_data_bytes))

                slice_counter += 1

        print(f"Total image slices: {slice_counter}")


def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(
        description='Convert multiple NIFTI/NRRD files in a directory to a single MRD file.')
    parser.add_argument('--inputdir', type=str, required=True, help='Path to the directory containing NIFTI/NRRD files')
    parser.add_argument('--outputfilename', type=str, help='Name of the output MRD file (without directory)')

    # Parse the arguments
    args = parser.parse_args()

    # Check input directory and gather NIFTI/NRRD files
    input_files = gather_input_files(args.inputdir)

    # Determine the output path
    if args.outputfilename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_basename = f"output_{timestamp}.mrd"
    else:
        output_basename = args.outputfilename

    # Save the MRD file one directory up from the input directory
    parent_dir = os.path.dirname(args.inputdir)  # Get the parent directory of the input directory
    grandparent_dir = os.path.dirname(parent_dir)  # Get the parent directory of the parent directory
    output_file = os.path.join(grandparent_dir, output_basename)

    # Call the function to create the MRD file
    create_mrd_file(input_files, args.inputdir, output_file)

if __name__ == '__main__':
    main()




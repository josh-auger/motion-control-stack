# Title: check_mrd_file.py

# Description:
# Brief script to confirm contents of MRD file after converting a set of 3D NIFTI volumes.
#
# Example use:
# python3 check_mrd.py --mrdfile /path/to/input/file.mrd

# Created on: August 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import ismrmrd
import argparse

def check_mrd_file(mrd_file_path):
    try:
        # Open the MRD file
        dataset = ismrmrd.Dataset(mrd_file_path, create_if_needed=False)

        # Get the list of image groups
        image_groups = dataset.list()

        number_of_images = 0

        # Iterate over all groups and count images
        for group in image_groups:
            if 'image' in group:
                number_of_images += 1
                image = dataset.read_image(group)
                image_size = image.data.size
                print(f"Image {number_of_images}: Size = {image_size} pixels")

        print(f"Total number of images in the MRD file: {number_of_images}")

    except Exception as e:
        print(f"An error occurred: {e}")


def print_mrd_metadata(mrd_file_path):
    try:
        # Open the MRD file
        dataset = ismrmrd.Dataset(mrd_file_path, create_if_needed=False)

        # Print the XML header (metadata)
        xml_header = dataset.read_xml_header()
        print("MRD XML Metadata:")
        print(xml_header.decode('utf-8'))  # Decode from bytes to string

    except Exception as e:
        print(f"An error occurred while reading metadata: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check the number of images and their sizes in an MRD file, and print MRD XML metadata.")
    parser.add_argument("--mrdfile", required=True, help="Path to the input MRD file")
    parser.add_argument("--print_metadata", action="store_true", help="Print MRD XML metadata")

    args = parser.parse_args()

    check_mrd_file(args.mrdfile)

    if args.print_metadata:
        print_mrd_metadata(args.mrdfile)

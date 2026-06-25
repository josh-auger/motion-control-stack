"""
Title: update_image_frame_with_moco_transform.py

Description:
    Test script for updating NRRD image header metadata (origin and directions) with the applied moco feedback in order
    to correctly place the image data in the "true" position in global space.

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Created on: 24 September 2025

"""

import SimpleITK as sitk
import numpy as np
import argparse


def update_image_frame_with_moco_transform(origin, scaled_direction_matrix, transform):
    """
    Update the image header metadata to account for the moco feedback applied by the scanner, placing the image
    data in its true position in global space.
    Moco rotation gets applied about the device isocenter, altering the image origin location and image directions.
    Moco translation gets applied to the image origin.
    """
    print(f"Updating header metadata to account for applied moco feedback :")
    R = np.array(transform.GetMatrix()).reshape((3, 3))
    t = np.array(transform.GetTranslation())
    # Try to get center
    if hasattr(transform, "GetCenter"):
        c = np.array(transform.GetCenter())
    else:
        fixed_params = transform.GetFixedParameters()
        if len(fixed_params) >= 3:
            c = np.array(fixed_params[:3])
        else:
            c = np.zeros(3)

    # Apply moco rotation (R) and translation (t) to the image header direction cosine matrix and origin, respectively
    # R_inv = np.array((moco_transform.GetInverse()).GetMatrix()).reshape((3,3))   # Get inverse rotation matrix directly from transform object
    R_inv = np.linalg.inv(R)
    print(f"R_inv :")
    print(f"{R_inv}")
    print(f"About center of rotation : {c}\n")

    updated_direction_matrix = R_inv @ scaled_direction_matrix
    updated_origin = (origin - c) @ R + c - t
    print(f"updated_origin : {updated_origin}")
    print(f"updated_direction_matrix : {updated_direction_matrix}\n")

    return updated_origin, updated_direction_matrix

def main():
    parser = argparse.ArgumentParser(description="Apply a transform to the spatial positioning metadata in the image header.")
    parser.add_argument("--refimage", required=True, help="Input reference image used to define starting image frame origin and direction vectors.")
    parser.add_argument("--transform", required=True, help="Input transform to be applied to the image header.")
    args = parser.parse_args()

    # Read input reference image
    refimage = sitk.ReadImage(args.refimage)
    origin = np.array(refimage.GetOrigin())
    direction_matrix = np.array(refimage.GetDirection()).reshape(3, 3)
    spacing = np.array(refimage.GetSpacing())
    scaled_direction_matrix = direction_matrix * spacing
    dir_x = scaled_direction_matrix[:, 0]
    dir_y = scaled_direction_matrix[:, 1]
    dir_z = scaled_direction_matrix[:, 2]
    print(f"Reference origin : {origin}")
    print("Reference direction matrix :")
    print(direction_matrix)
    print(f"Reference spacing : {spacing}")
    print("Reference image direction vectors :")
    print(f"X direction : {dir_x}")
    print(f"Y direction : {dir_y}")
    print(f"Z direction : {dir_z}\n")

    # Read input transform
    read_transform = sitk.ReadTransform(args.transform)
    print(f"Input transform to apply :")
    print(read_transform)

    new_origin, new_directions = update_image_frame_with_moco_transform(origin, scaled_direction_matrix, read_transform)
    print(f"New image header metadata :")
    print(f"New space origin: ({new_origin[0]}, {new_origin[1]}, {new_origin[2]})")
    print(f"New space directions: ({new_directions[:, 0][0]}, {new_directions[:, 0][1]}, {new_directions[:, 0][2]})\
    ({new_directions[:, 1][0]}, {new_directions[:, 1][1]}, {new_directions[:, 1][2]})\
    ({new_directions[:, 2][0]}, {new_directions[:, 2][1]}, {new_directions[:, 2][2]})")


if __name__ == "__main__":
    main()
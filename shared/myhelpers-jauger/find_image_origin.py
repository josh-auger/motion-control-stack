"""
Title: calculate_image_origin.py

Description:
    Calculate the spatial location of the first (top, left) voxel of an ismrmrd image using image header metadata.

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Created on: 25 June 2025

"""

import SimpleITK as sitk
import numpy as np
import argparse
import ctypes
import time


def get_first_voxel_center(image_center, image_size, direction_matrix):
    # Calculate the spatial location of the center of the first voxel (top left) of an image
    center_indices = (image_size - 1) / 2.0   # indices of the center voxel
    print(f"\tCenter indices : {center_indices}")
    first_voxel_center = image_center - direction_matrix @ center_indices
    print(f"\tFirst voxel center : {first_voxel_center}")

    return first_voxel_center


def main():
    print(f"Calculating center of first voxel as origin...")
    center = np.array([2.4736549854278564, 1.1560499668121338, -101.07147979736328])
    size = np.array([86, 86, 1])
    spacing = np.array([2.9767441860465116, 2.9767441860465116, 3.0])

    print(f"\tCenter of image : {center}")
    print(f"\tSize of image : {size}")
    print(f"\tSpacing of image : {spacing}")

    dir_x = np.array([-2.969799130461937, 0.0, -0.2032223191372184])
    dir_y = np.array([0.0, 2.9767441860465116, 0.0])
    dir_z = np.array([-0.2048099935054779, 0.0, 2.9930006861686707])

    print(f"\tX direction of image : {dir_x}")
    print(f"\tY direction of image : {dir_y}")
    print(f"\tZ direction of image : {dir_z}")

    # direction x spacing matrix, shape (3,3)
    D = np.column_stack((
        dir_x,
        dir_y,
        dir_z))

    print(f"\tDirection matrix : {D}")

    space_origin = get_first_voxel_center(center, size, D)

if __name__ == "__main__":
    main()
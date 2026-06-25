"""
Title: convert_transform.py

Description:
    Helper script to convert transforms from one type to another.

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Created on: 8 August 2025

"""

import SimpleITK as sitk
import numpy as np
import argparse
import ctypes
import time
import os
import math
from scipy.spatial.transform import Rotation as R


def construct_transform_object(transform, transformtype):
    rotation_matrix = np.array(transform.GetMatrix()).reshape(3, 3)
    translation = transform.GetTranslation()
    center = transform.GetCenter()

    if transformtype == "Affine":
        new_transform = sitk.AffineTransform(3)
        new_transform.SetCenter(center)
        new_transform.SetMatrix(rotation_matrix.ravel())
        new_transform.SetTranslation(translation)

    elif transformtype == "VersorRigid3D":
        new_transform = sitk.VersorRigid3DTransform()
        new_transform.SetCenter(center)
        new_transform.SetMatrix(rotation_matrix.ravel())
        new_transform.SetTranslation(translation)

    elif transformtype == "Euler3D":
        new_transform = sitk.Euler3DTransform()
        new_transform.SetCenter(center)
        new_transform.SetMatrix(rotation_matrix.ravel())
        new_transform.SetTranslation(translation)

    else:
        raise ValueError(f"Unsupported transform type: {transformtype}")

    return new_transform


def main():
    parser = argparse.ArgumentParser(description="Convert a 3D rigid transform to an affine transform.")
    parser.add_argument("--transform", required=True, help="Input transform for conversion.")
    parser.add_argument("--convertto", required=True, choices=["Affine", "VersorRigid3D", "Euler3D"], help="Target transform type.")
    parser.add_argument("--outputtransform", required=False, help="Output filename for the converted transform.")
    args = parser.parse_args()

    # Read transform
    read_transform = sitk.ReadTransform(args.transform)

    # Format helper for printing floats to 6 decimals
    fmt = lambda arr: np.array(arr, dtype=float).round(6).tolist()

    print(f"Input transform type: {read_transform.GetName()}")
    print(f"Input transform center: {fmt(read_transform.GetCenter())}")
    print(f"Input transform parameters: {fmt(read_transform.GetParameters())}")

    # Convert transform to new type
    output_transform = construct_transform_object(read_transform, args.convertto)

    print(f"Converted transform type: {args.convertto}")
    print(f"Converted transform center: {fmt(output_transform.GetCenter())}")
    print(f"Converted transform parameters: {fmt(output_transform.GetParameters())}")

    # Write converted transform
    if args.outputtransform:
        output_filename = args.outputtransform
    else:
        base, ext = os.path.splitext(args.transform)
        output_filename = f"{base}_{args.convertto}.tfm"

        # Write converted transform
    sitk.WriteTransform(output_transform, output_filename)
    print(f"Transform written to: {output_filename}")

if __name__ == "__main__":
    main()
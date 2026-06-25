#!/usr/bin/env python3

import SimpleITK as sitk
import numpy as np
import argparse
import os

#
# This program reads in a transform and computes its inverse.
#

# Parse the arguments.
parser = argparse.ArgumentParser(description='invert an affine transform.')
parser.add_argument("--inputtransformfile", required=True)
parser.add_argument("--outputtransformfile", required=False)
args = parser.parse_args()


# Review the definition of transformations here:
# https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/22_Transforms.html
# https://simpleitk.org/SPIE2019_COURSE/01_spatial_transformations.html
#
# T(x)=A(x−c)+t+c
# A is a 3x3 matrix.
# t is a 3x1 vector.
# c is the center of rotation, a point represented by a 3x1 vector.

# input transform
read_transform1 = sitk.ReadTransform( args.inputtransformfile )
print("input transform:")
print(read_transform1)

# inverse transform, as defined by SimpleITK
read_transform1_inverse = read_transform1.GetInverse()
print("inverse transform (from SimpleITK):")
print(read_transform1_inverse)
print(read_transform1_inverse)

# inverse transform as "equal and opposite"
R_inverse = read_transform1_inverse.GetMatrix() # inverse rotation of original transform
t_original = read_transform1.GetTranslation()    # same translation as original transform
t_opposite = tuple(-x for x in t_original)        # flip the sign
read_transform1_opposite = read_transform1_inverse
read_transform1_opposite.SetMatrix(R_inverse)
read_transform1_opposite.SetTranslation(t_opposite)
print("opposite transform (inverse rotation, flipped translation):")
print(read_transform1_opposite)

# Write converted transform
if args.outputtransform:
    output_filename = args.outputtransformfile
else:
    base, ext = os.path.splitext(args.inputtransformfile)
    output_filename = f"{base}_inverse.tfm"

sitk.WriteTransform( read_transform1_opposite, output_filename )

exit(1)
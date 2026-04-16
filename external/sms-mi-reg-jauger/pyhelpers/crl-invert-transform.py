#!/usr/bin/env python3

import SimpleITK as sitk
import numpy as np
import argparse

#
# This program reads in a transform and computes its inverse.
#

# Parse the arguments.
parser = argparse.ArgumentParser(description='invert an affine transform.')
parser.add_argument("--inputtransformfile", required=True)
parser.add_argument("--outputtransformfile", required=True)
args = parser.parse_args()


# Review the definition of transformations here:
# https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/22_Transforms.html
# https://simpleitk.org/SPIE2019_COURSE/01_spatial_transformations.html
#
# T(x)=A(x−c)+t+c
# A is a 3x3 matrix.
# t is a 3x1 vector.
# c is the center of rotation, a point represented by a 3x1 vector.

read_transform1 = sitk.ReadTransform( args.inputtransformfile )
print(read_transform1)

read_transform1_inverse = read_transform1.GetInverse()
print(read_transform1_inverse)

sitk.WriteTransform( read_transform1_inverse, args.outputtransformfile )

exit(1)


#!/usr/bin/env python3

import SimpleITK as sitk
import numpy as np
import argparse
import json

#
# This program reads in an affine transform and a desired center of rotation.
# It computes the equivalent affine transformation with the specified center
# of rotation.
# The new transform is written to a file.
#

# Example usage:
# crl-recenter-affine-transform.py --transformfile \
#      affine-ozero-to-identity-directions.tfm --center "[0,0,0]" \
#      --outputfile tst.tfm

# We believe that angles are always in radians in SimpleITK.

# Parse the arguments.
parser = argparse.ArgumentParser(description='Adjust transform to a specific center.')
parser.add_argument("--transformfile", required=True)
parser.add_argument("--center", required=True, type=json.loads)
parser.add_argument("--outputfile", required=True)
args = parser.parse_args()


# Review the definition of transformations here:
# https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/22_Transforms.html
# https://simpleitk.org/SPIE2019_COURSE/01_spatial_transformations.html
#
# T(x)=A(x−c)+t+c
# A is a 3x3 matrix.
# t is a 3x1 vector.
# c is the center of rotation, a point represented by a 3x1 vector.

read_transform = sitk.ReadTransform( args.transformfile )
print(read_transform)

A0 = np.asarray(read_transform.GetMatrix()).reshape(3,3)
c0 = np.asarray(read_transform.GetCenter())
t0 = np.asarray(read_transform.GetTranslation())

# Create a transform manually. 
#
# T_0(A0,c0,t0)(x) = A0(x-c0) + c0 + t0
# T_0(A0,c0,t0)(x) = A0(x) -A0c0 + c0 + t0
#
# T_1(A0,c1,t1)(x) = A0(x) -A0c1 + c1 + t1
#     t1 + c1 - A0c1 = t0 + c0 -A0c0
#     t1 = t0 + c0 - A0c0 + A0c1 - c1
#     t1 = t0 + c0 - c1 + A0c1 - A0c0
#     t1 = t0 + A0(c1 - c0) - (c1 - c0)

print(args.center)

combined_mat = A0
# combined_center =  np.zeros(3)
c1 = args.center
print(c1)
combined_center = c1
center_difference = c1 - c0
combined_translation = t0 + np.dot(A0, center_difference) - center_difference
combined_affine = sitk.AffineTransform(combined_mat.flatten(), combined_translation, combined_center)

print('Matrix A0')
print(A0)
print('Combined center')
print(combined_center)
print('Translation')
print(combined_translation)

print('combined affine transform is : ')
print(combined_affine)
print(combined_affine.GetParameters())
print(combined_affine.GetFixedParameters())

sitk.WriteTransform( combined_affine, args.outputfile )

versorrigid3d = sitk.VersorRigid3DTransform()
versorrigid3d.SetCenter(combined_center)
versorrigid3d.SetTranslation(combined_translation)
versorrigid3d.SetMatrix(combined_mat.flatten())
print(versorrigid3d)
print(versorrigid3d.GetParameters())

# First three parameters are rotation angles in radians.
# Second three parameters are translations.

euler3d = sitk.Euler3DTransform()
euler3d.SetCenter(combined_center)
euler3d.SetTranslation(combined_translation)
euler3d.SetMatrix(combined_mat.flatten())
print(euler3d)
print(euler3d.GetParameters())

exit(1)


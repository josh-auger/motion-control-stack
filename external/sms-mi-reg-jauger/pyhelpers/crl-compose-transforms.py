#!/usr/bin/env python3

import SimpleITK as sitk
import numpy as np
import argparse

#
# This program reads in two affine transforms.
# Each transform maps from a reference volume to a target volume.
# The second transform maps from a reference volume to a target volume.
# The first transform is inverted and composed with the second transform to 
#   create a transform that maps from the 
#       first target volume to the second target volume.
# The composed transform is written to a file.
#

# We believe that angles are always in radians in SimpleITK.

# Parse the arguments.
parser = argparse.ArgumentParser(description='compose two affine transforms.')
parser.add_argument("--transformfile1", required=True)
parser.add_argument("--transformfile2", required=True)
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

read_transform1 = sitk.ReadTransform( args.transformfile1 )
print(read_transform1)

read_transform1_inverse = read_transform1.GetInverse()
print(read_transform1_inverse)

read_transform2 = sitk.ReadTransform( args.transformfile2 )
print(read_transform2)

A0 = np.asarray(read_transform2.GetMatrix()).reshape(3,3)
c0 = np.asarray(read_transform2.GetCenter())
t0 = np.asarray(read_transform2.GetTranslation())

A1 = np.asarray(read_transform1_inverse.GetMatrix()).reshape(3,3)
c1 = np.asarray(read_transform1_inverse.GetCenter())
t1 = np.asarray(read_transform1_inverse.GetTranslation())

# Create a single transform manually. 
# this is a recipe for compositing any two global transformations
# into an affine transformation, T_0(T_1(x)):
# A = A=A0*A1
# c = c1
# t = A0*[t1+c1-c0] + t0+c0-c1
#
# T_0(A0,c0,t0)(x) = A0(x-c0) + c0 + t0
# T_1(A1,c1,t1)(x) = A1(x-c1) + c1 + t1
# T_0(T_1(x)) = A0( A1(x-c1) + c1 + t1 - c0 ) + c0 + t0
#             = A0A1(x-c1) + A0(c1 + t1 - c0) + c0 + t0
#             = A0A1(x-c1) + c1 - c1 + A0(c1 + t1 - c0) + c0 + t0
#             = A0A1(x-c1) + c1 + t
#               where t = A0(t1 + c1 - c0) + t0 + c0 - c1
#

combined_mat = np.dot(A0,A1)
combined_center = c1
combined_translation = np.dot(A0, t1+c1-c0) + t0+c0-c1
combined_affine = sitk.AffineTransform(combined_mat.flatten(), combined_translation, combined_center)

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

# Compute the displacement:
radius = 50
parms = np.asarray( euler3d.GetParameters() )
print(parms)
displacement = abs(parms[0]*radius) + abs(parms[1]*radius) + \
    abs(parms[2]*radius) + abs(parms[3]) + abs(parms[4]) + abs(parms[5])
print(displacement)


exit(1)


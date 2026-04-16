#!/usr/bin/env python3

#
# pip3 install SimpleITK
# pip3 install NumPy
# pip3 install matplotlib
# sudo pip3 install scipy
#

import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import sys
import argparse

parser = argparse.ArgumentParser(description="Rotate the space directions of an image.")
parser.add_argument("--input_image", required=True)
parser.add_argument("--transform", required=True)
parser.add_argument("--output_image", required=True)
args = parser.parse_args()

inVolume1 = sitk.ReadImage( args.input_image )
print("Origin is : " + str(inVolume1.GetOrigin()) )
D0 = inVolume1.GetDirection()
print("Directions are : " + str(D0))

D00 = D0[0:3]
D01 = D0[3:6]
D02 = D0[6:9]

print("D00 is : " + str(D00))
print("D01 is : " + str(D01))
print("D02 is : " + str(D02))

read_transform1 = sitk.ReadTransform( args.transform )
# print(read_transform1)

# Review the definition of transformations here:
# https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/22_Transforms.html
# https://simpleitk.org/SPIE2019_COURSE/01_spatial_transformations.html
#
# T(x)=A(x−c)+t+c
# A is a 3x3 matrix.
# t is a 3x1 vector.
# c is the center of rotation, a point represented by a 3x1 vector.

A0 = np.asarray(read_transform1.GetMatrix()).reshape(3,3)


print('Rotation is ' + str(A0))


Frame00 = np.dot(A0, D00)
print('Rotated e1 is ' + str(Frame00))

Frame01 = np.dot(A0, D01)
print('Rotated e2 is ' + str(Frame01))

Frame02 = np.dot(A0, D02)
print('Rotated e3 is ' + str(Frame02))

E0 = np.concatenate([Frame00, Frame01, Frame02], axis=None)

print("Rotated directions are : " + str(E0))


inVolume1.SetDirection( E0 )
E0 = inVolume1.GetDirection()

writer = sitk.ImageFileWriter()
writer.SetFileName( args.output_image )
writer.Execute( inVolume1 )

quit(0)


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

parser = argparse.ArgumentParser(description="Convert a 3D rigid transform to an affine transform.")
parser.add_argument("--transform", required=True)
parser.add_argument("--affineTransform", required=True)
args = parser.parse_args()

read_transform1 = sitk.ReadTransform( args.transform )

# Review the definition of transformations here:
# https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/22_Transforms.html
# https://simpleitk.org/SPIE2019_COURSE/01_spatial_transformations.html
#
# T(x)=A(x−c)+t+c
# A is a 3x3 matrix.
# t is a 3x1 vector.
# c is the center of rotation, a point represented by a 3x1 vector.

A0 = np.asarray(read_transform1.GetMatrix()).reshape(3,3)
print(A0)

affine_transform = sitk.AffineTransform( 3 )
affine_transform.SetIdentity()
affine_transform.SetCenter( read_transform1.GetCenter() )
affine_transform.SetTranslation( read_transform1.GetTranslation() )
affine_transform.SetMatrix( read_transform1.GetMatrix() )

sitk.WriteTransform( affine_transform, args.affineTransform )


quit(0)


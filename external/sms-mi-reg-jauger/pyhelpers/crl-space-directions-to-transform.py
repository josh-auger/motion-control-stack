#!/usr/bin/env python3

#
# pip3 install SimpleITK
# pip3 install NumPy
# pip3 install matplotlib
# sudo pip3 install scipy
#

import SimpleITK as sitk
import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser(description="Rotate the space directions of an image.")
parser.add_argument("--input", required=True)
parser.add_argument("--transformfile", required=True)
args = parser.parse_args()

inVolume1 = sitk.ReadImage( args.input )
print("Origin is : " + str(inVolume1.GetOrigin()) )
D0 = inVolume1.GetDirection()
print("Directions are : " + str(D0))

D00 = D0[0:3]
D01 = D0[3:6]
D02 = D0[6:9]

print("D00 is : " + str(D00))
print("D01 is : " + str(D01))
print("D02 is : " + str(D02))

transform = sitk.AffineTransform(3)
transform.SetIdentity()
transform.GetMatrix()

# T(x)=A(x−c)+t+c
# A is a 3x3 matrix.
# t is a 3x1 vector.
# c is the center of rotation, a point represented by a 3x1 vector.

A0 = np.asarray(transform.GetMatrix()).reshape(3,3)
print(A0)

transform.SetMatrix(D0)
print(transform)

sitk.WriteTransform( transform, args.transformfile )

quit(0)


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

parser = argparse.ArgumentParser(description="Translate the origin of an image.")
parser.add_argument("--input_image", required=True)
parser.add_argument("--transform", required=True)
parser.add_argument("--output_image", required=True)
args = parser.parse_args()

inVolume1 = sitk.ReadImage( args.input_image )
print("Origin is : " + str(inVolume1.GetOrigin()) )

read_transform1 = sitk.ReadTransform( args.transform )
# print(read_transform1)

t0 = np.asarray(read_transform1.GetTranslation())
print('Transformation translation is ' + str(t0))

newOrigin = read_transform1.TransformPoint( inVolume1.GetOrigin() )

print("new origin is : " + str(newOrigin) )

inVolume1.SetOrigin( newOrigin )

writer = sitk.ImageFileWriter()
writer.SetFileName( args.output_image )
writer.Execute( inVolume1 )

quit(0)


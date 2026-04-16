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

parser = argparse.ArgumentParser(description="Translate the origin of an image to be located at the geometrical center.")
parser.add_argument("--input_image", required=True)
parser.add_argument("--output_image", required=True)
args = parser.parse_args()

inVolume1 = sitk.ReadImage( args.input_image )
print("Origin is : " + str(inVolume1.GetOrigin()) )

# Center of first voxel.
image_origin = inVolume1.TransformContinuousIndexToPhysicalPoint(
    [0 for index in inVolume1.GetSize()] )
print('image_origin is : ' + str(image_origin) )

# Center of middle of volume.
image_center = inVolume1.TransformContinuousIndexToPhysicalPoint(
    [(index-1)/2.0 for index in inVolume1.GetSize()] )
print('image_center is : ' + str(image_center) )

# Center of furthest voxel.
image_far = inVolume1.TransformContinuousIndexToPhysicalPoint(
    [(index-1) for index in inVolume1.GetSize()] )
print('image_far is : ' + str(image_far) )

print('Set the new origin to the image center: ' + str(image_center) )
newOrigin = image_center

inVolume1.SetOrigin( newOrigin )

writer = sitk.ImageFileWriter()
writer.SetFileName( args.output_image )
writer.Execute( inVolume1 )

quit(0)


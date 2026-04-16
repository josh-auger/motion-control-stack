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

# Parse the arguments.
parser = argparse.ArgumentParser(description='resample volume to slice.')
parser.add_argument("--refvolume", required=True)
parser.add_argument("--slicetarget", required=True)
parser.add_argument("--transformfile", required=True)
parser.add_argument("--outputfile", required=True)
args = parser.parse_args()

refVolume = args.refvolume
sliceTarget = args.slicetarget
transformFile = args.transformfile
outputImageFileName = args.outputfile

reader = sitk.ImageFileReader()
reader.SetFileName( refVolume )
refImage = reader.Execute();

# Center of first voxel.
image_origin = refImage.TransformContinuousIndexToPhysicalPoint(
    [0 for index in refImage.GetSize()] )
print('image_origin is : ' + str(image_origin) )

# Center of middle of volume.
image_center = refImage.TransformContinuousIndexToPhysicalPoint(
    [(index-1)/2.0 for index in refImage.GetSize()] )
print('image_center is : ' + str(image_center) )
 
# Center of furthest voxel.
image_far = refImage.TransformContinuousIndexToPhysicalPoint(
    [(index-1) for index in refImage.GetSize()] )
print('image_far is : ' + str(image_far) )

sliceReader = sitk.ImageFileReader()
sliceReader.SetFileName( sliceTarget )
sliceImage = sliceReader.Execute();

transform = sitk.ReadTransform( transformFile )

if transform is None:
  transform = sitk.Transform()
  transform.SetIdentity()

output_pixel_type = sliceReader.GetPixelID()
resampleFilter = sitk.ResampleImageFilter()
resampleFilter.SetInterpolator( sitk.sitkLinear )
resampleFilter.SetTransform( transform )
resampleFilter.SetOutputPixelType( output_pixel_type )
resampleFilter.SetDefaultPixelValue( 0.0 )
resampleFilter.SetReferenceImage( sliceImage )
newImage = resampleFilter.Execute( refImage )

writer = sitk.ImageFileWriter()
writer.SetFileName( outputImageFileName )
writer.Execute( newImage )

quit(0)


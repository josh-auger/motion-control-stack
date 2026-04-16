#!/usr/bin/env python3

#
# pip3 install SimpleITK
# pip3 install NumPy
# pip3 install matplotlib
# sudo pip3 install scipy
#

import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy.matlib
import numpy as np
import scipy.stats as st
import sys
import argparse

# Parse the arguments.
parser = argparse.ArgumentParser(description='resample volume to slice.')
parser.add_argument("--refvolume", required=True)
parser.add_argument("--outputfile", required=True)
parser.add_argument("--outputtransformfile", required=True)
args = parser.parse_args()

def duplicate_image( imageName ):
  newImage = sitk.Image( imageName.GetSize(), imageName.GetPixelID(), 
                         imageName.GetNumberOfComponentsPerPixel() )
  newImage.SetOrigin( imageName.GetOrigin() )
  newImage.SetSpacing( imageName.GetSpacing() )
  # There is some question about what the assumed axes are e.g. RAI,LPS
  newImage.SetDirection( imageName.GetDirection() )
  return newImage


# When manipulating the space directions, it is important to account for
# the ordering of the data on disk. 
# This relates to the ordering of axes in space. 
# An ITK image with an identity direction cosine matrix is in LPS+
# (Left, Posterior, Superior) orientation as defined by the DICOM standard.
# https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1DICOMOrientImageFilter.html#details
# The DICOMOrient filter can be used to set the axes without changing the 
# coordinates of any voxel.

#
#
def create_image_same_size_default_directions( imageName ):
  newImage = sitk.Image( imageName.GetSize(), imageName.GetPixelID(), 
                         imageName.GetNumberOfComponentsPerPixel() )
  newImage.SetOrigin( imageName.GetOrigin() )
  newImage.SetSpacing( imageName.GetSpacing() )
  newImage.SetDirection( imageName.GetDirection() )
  # There is some question about what the assumed axes are e.g. RAI,LPS
  # We use the DICOMOrient filter to change the space directions.
  # This is done in such a way that the geometrical locations of the voxels
  # are unchanged.
  newImage = sitk.DICOMOrient( newImage, 'LPS')
  return newImage


refVolume = args.refvolume
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

# Now we clone the input image.
clonedImage = duplicate_image( refImage )
clonedImage = create_image_same_size_default_directions( refImage )
D0 = clonedImage.GetDirection()
print('Direction of cloned Image:')
print(D0)

D1 = np.matlib.identity(3, dtype = float)
print(D1)
D2 = np.asarray(D1).flatten()
print('Final directions:')
print(tuple(D2))
clonedImage.SetDirection( D2 )

transform = sitk.AffineTransform(3)
transform.SetIdentity()
transform.GetMatrix()

# Review the definition of transformations here:
# https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/22_Transforms.html
# https://simpleitk.org/SPIE2019_COURSE/01_spatial_transformations.html
#
# T(x)=A(x−c)+t+c
# A is a 3x3 matrix.
# t is a 3x1 vector.
# c is the center of rotation, a point represented by a 3x1 vector.

#1. Get the rotation matrix from the transform
A0 = np.asarray(transform.GetMatrix()).reshape(3,3)

print('Identity transform is : ')
print(A0)

#2. Record the space directions rotation.
print('SDm')
print('Space directions transform is : ')
SDm = np.asarray(D0).reshape(3,3)
print(SDm)
transform.SetMatrix(SDm.reshape(-1))

print('SDmi')
SDmi = np.linalg.inv(SDm)
print(SDmi)

print(np.matmul(SDmi, SDm))

A1 = SDmi
print('Inverse space directions transform is : ')
print(A1)

print(transform)

sitk.WriteTransform( transform, args.outputtransformfile )

output_pixel_type = refImage.GetPixelID()
resampleFilter = sitk.ResampleImageFilter()
resampleFilter.SetInterpolator( sitk.sitkLinear )
resampleFilter.SetTransform( transform )
resampleFilter.SetOutputPixelType( output_pixel_type )
resampleFilter.SetDefaultPixelValue( 0.0 )
resampleFilter.SetReferenceImage( clonedImage )
newImage = resampleFilter.Execute( refImage )

writer = sitk.ImageFileWriter()
writer.SetFileName( outputImageFileName )
writer.Execute( newImage )

quit(0)


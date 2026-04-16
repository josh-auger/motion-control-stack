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
import sys
import argparse
import logging

def construct_affine_identity_transform(refVolume, transformFile):

    reader = sitk.ImageFileReader()
    reader.SetFileName( refVolume )
    refImage = reader.Execute();

    # Center of first voxel.
    image_origin = refImage.TransformContinuousIndexToPhysicalPoint(
        [0 for index in refImage.GetSize()] )
    # print('image_origin is : ' + str(image_origin) )
    logging.info('\timage_origin is : ' + str(image_origin) )

    # Center of middle of volume.
    image_center = refImage.TransformContinuousIndexToPhysicalPoint(
        [(index-1)/2.0 for index in refImage.GetSize()] )
    # print('image_center is : ' + str(image_center) )
    logging.info('\timage_center is : ' + str(image_center) )

    # Center of furthest voxel.
    image_far = refImage.TransformContinuousIndexToPhysicalPoint(
        [(index-1) for index in refImage.GetSize()] )
    # print('image_far is : ' + str(image_far) )
    logging.info('\timage_far is : ' + str(image_far) )

    transform = sitk.AffineTransform( 3 )
    transform.SetIdentity()
    transform.SetCenter( image_center )

    sitk.WriteTransform( transform, transformFile )

    return transform


if __name__ == "__main__":
    # Parse the arguments.
    parser = argparse.ArgumentParser(description='construct an affine identity transform with a center of rotation at the geometric center of a reference volume.')
    parser.add_argument("--refvolume", required=True)
    parser.add_argument("--transformfile", required=True)
    args = parser.parse_args()

    construct_affine_identity_transform(args.refvolume, args.transformfile)
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
parser = argparse.ArgumentParser(description='Change the default orientation.')
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--order", default="LPS", required=False)
args = parser.parse_args()

inputVolume = args.input
outputImageFileName = args.output

reader = sitk.ImageFileReader()
reader.SetFileName( inputVolume )
inputImage = reader.Execute();

# Now we clone the input image.
reorientedImage = sitk.DICOMOrient( inputImage, args.order )

writer = sitk.ImageFileWriter()
writer.SetFileName( outputImageFileName )
writer.Execute( reorientedImage )

quit(0)


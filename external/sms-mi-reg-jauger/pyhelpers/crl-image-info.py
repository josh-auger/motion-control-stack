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

parser = argparse.ArgumentParser(description="Display image meta data.")
parser.add_argument("--input", nargs='+' )

args = parser.parse_args()

for imagename in args.input:

  print('Info for : ' + str(imagename))
  image = sitk.ReadImage( imagename )

  print('Size: ' + str(image.GetSize()) )
  print('Origin: ' + str(image.GetOrigin()) )
  print('Spacing: ' + str(image.GetSpacing()) )
  print('Direction: ' + str(image.GetDirection()) )
  print('Number of components per pixel: ' + str(image.GetNumberOfComponentsPerPixel()) )


quit(0)


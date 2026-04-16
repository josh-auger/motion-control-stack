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

parser = argparse.ArgumentParser(description="Checkerboard combine images.")
parser.add_argument("--input_image1")
parser.add_argument("--input_image2")
parser.add_argument("--output_image")
args = parser.parse_args()

inVolume1 = sitk.ReadImage( args.input_image1 )
inVolume2 = sitk.ReadImage( args.input_image2 )

inVolume2 = sitk.Cast( inVolume2, sitk.sitkInt16 )

newImage = sitk.CheckerBoard(inVolume1, inVolume2, [8,8,8])

writer = sitk.ImageFileWriter()
writer.SetFileName( args.output_image )
writer.Execute( newImage )

quit(0)


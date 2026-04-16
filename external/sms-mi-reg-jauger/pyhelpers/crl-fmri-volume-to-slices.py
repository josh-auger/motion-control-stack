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
from pathlib import Path
from pathlib import PurePath


parser = argparse.ArgumentParser(description="Extract a subregion of an image.")
parser.add_argument("input")
parser.add_argument("output")
args = parser.parse_args()

inVolume = sitk.ReadImage( args.input )
dims = inVolume.GetDimension()
if dims != 3:
  print('Expecting the number of dimensions to be 3.')
  exit()

sizes = inVolume.GetSize()
for i in range(0,sizes[2]):
  sliceIndex = i
  startIndex = (0, 0, sliceIndex)
  sizeROI = (sizes[0], sizes[1], 1)
  roiVolume = sitk.RegionOfInterest(inVolume, sizeROI, startIndex )
  outName = PurePath(args.output).stem + '-' + str(sliceIndex).zfill(3) + PurePath(args.output).suffix
  print(outName)

  writer = sitk.ImageFileWriter()
  writer.SetFileName( outName )
  writer.Execute( roiVolume )

quit(0)


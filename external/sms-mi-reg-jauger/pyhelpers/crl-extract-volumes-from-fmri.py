#!/usr/bin/env python3

#
# pip3 install SimpleITK
# pip3 install NumPy
# pip3 install matplotlib
# sudo pip3 install scipy
#

# Example run command:
# python3 crl-extract-volumes-from-fmri /path/to/4D/nifti/file output_filename.nii

import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import sys
import argparse
from pathlib import Path
from pathlib import PurePath

parser = argparse.ArgumentParser(description="Convert a 4D image into a sequence of 3D images.")
parser.add_argument("input", help="File path to the input 4D image (e.g., /local/path/to/input_4d_image.nii.gz)")
parser.add_argument("output", help="Base filename for output 3D images with desired file extension (e.g., output_3d_image.nrrd)")
args = parser.parse_args()

inVolume = sitk.ReadImage( args.input )
print(inVolume)
dims = inVolume.GetDimension()
if dims != 4:
    print('Error: Input image must be 4D.')
    sys.exit(1)

sizes = inVolume.GetSize()
smallerSizes = (sizes[0], sizes[1], sizes[2], 0)
print(smallerSizes)

for i in range(0,sizes[3]):
# Select same subregion using ExtractImageFilter
  extract = sitk.ExtractImageFilter()
  extract.SetSize( smallerSizes )
  setIndex = (0,0,0,i)
  extract.SetIndex( setIndex )
  smallerVolume = extract.Execute(inVolume)
  outName = PurePath(args.output).stem + '-' + str(setIndex[3]).zfill(4) + PurePath(args.output).suffix
  print(outName)
  writer = sitk.ImageFileWriter()
  writer.SetFileName( outName )
  writer.Execute( smallerVolume )
  print(smallerVolume)


quit(0)


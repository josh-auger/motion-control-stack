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

parser = argparse.ArgumentParser(description="Rotate the space directions of an image to be identity.")
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

inVolume1 = sitk.ReadImage( args.input )
print("Origin is : " + str(inVolume1.GetOrigin()) )
D0 = inVolume1.GetDirection()
print("Directions are : " + str(D0))

D1 = np.matlib.identity(3, dtype = float)
print(D1)
D2 = np.asarray(D1).flatten()
print('Final directions:')
print(tuple(D2))
inVolume1.SetDirection( D2 )

writer = sitk.ImageFileWriter()
writer.SetFileName( args.output )
writer.Execute( inVolume1 )

quit(0)


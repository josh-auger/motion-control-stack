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
import os
import argparse
import subprocess
from pathlib import Path
from pathlib import PurePath

parser = argparse.ArgumentParser(description="Prepare a DICOM fmri scan for motion analysis.")
parser.add_argument("input")
parser.add_argument("output")
args = parser.parse_args()

# INPUT: DICOM zip file
# OUTPUT: 
# 1. run dcm2niix on the DICOM file

# Convert the DICOM directory to a 4D nii file.
subprocess.run(["dcm2niix", "-b", "y", "-ba", "n", "-f", "fmri",
 str(args.input) ])

# 2. create 3D volumes from the FMRI.
subprocess.run(["crl-extract-volumes-from-fmri.py", "fmri.nii", 
    str(args.output) ])

# 3. extract slices from an FMRI volume.
subprocess.run(["crl-fmri-volume-to-slices.py", "align-0001.nii",
    "slices.nii" ])

# host computer : container directory alias
dirmapping = str(Path.cwd()) + ":" + "/data"
dockerprefix = ["docker","run","--rm", "-it", "--init", "-v", dirmapping,
    "--user", str(os.getuid())+":"+str(os.getgid())]

subprocess.run(dockerprefix + ["crl/crkit",
    "crlImageInfo", "align-0001.nii"] )



# 4.  Run some alignments:
subprocess.run([ "run-registration-example.py", 
  "--refvolume", "align-0000.nii", 
  "--outputtransformfile",  "trsf-align1",
  "--slicetarget",
  "/data/slices-000.nii",
  "/data/slices-015.nii",
  "/data/slices-030.nii",
  "/data/slices-045.nii"
  ])

subprocess.run([ "run-registration-example.py", 
  "--refvolume", "align-0000.nii", 
  "--outputtransformfile",  "trsf-align2",
  "--slicetarget",
  "/data/slices-001.nii",
  "/data/slices-016.nii",
  "/data/slices-031.nii",
  "/data/slices-046.nii"
  ])

subprocess.run([ "resample-ref-to-slice.py",
"--refvolume", "align-0000.nii",
"--slicetarget", "slices-001.nii",
"--transformfile", "sliceTransformtrsf-align2.tfm",
"--outputfile", "ref-aligned-to-slice-001.nii"
])


quit(0)


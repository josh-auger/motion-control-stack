#!/usr/bin/env python3 

import os
import subprocess
import tempfile
import sys
import logging
import glob
import argparse

dirmapping = os.getcwd() + ":" + "/data"

dockerprefix = ["sudo", "docker","run","--rm", "-it", "--init", "-v", dirmapping,
    "--user", str(os.getuid())+":"+str(os.getgid())]

print(' ')

fixedNumber = "003"
fixedName = "navigator" + fixedNumber + ".nrrd"

print('Generating the initial transformation with rotation in center of volume.')
subprocess.run( dockerprefix +
  [ "crl/sms-mi-reg", "crl-identity-transform-at-volume-center.py",
    "--refvolume", fixedName,
    "--transformfile", "identity-at-volume-center.tfm"] )

print('Running registration.')

for i in range(1,101,1):
  targetNumber = str(i).zfill(3)
  targetName="navigator" + targetNumber + ".nrrd"
  print(targetNumber)

  subprocess.run( dockerprefix +  
    [ "crl/crkit", "crlRigidRegistration",
    "-l", "identity-at-volume-center.tfm",
    fixedName,
    targetName,
    "rraligned-" + targetNumber + "-on-" + fixedNumber + ".nrrd",
    "trans-" + targetNumber + "-to-" + fixedNumber + ".tfm"
    ])

  subprocess.run( dockerprefix +  
  [ "crl/crkit", "crlResampler",
    "--dontOrientGeometryImage",
    targetName,
    "trans-" + targetNumber + "-to-" + fixedNumber + ".tfm",
    fixedName,
    "linear",
    "vol-" + targetNumber + "-on-" + fixedNumber + ".nrrd"
    ])

  subprocess.run( dockerprefix + 
    [ "crl/crkit", "crlImageAddMpyAdd", fixedName,
     "0.0", "1.0", "0.0", "rescaled-" + fixedName ])

# We can also run the registration in the other direction.
#subprocess.run( dockerprefix +  
#  [ "crl/crkit", "crlRigidRegistration",
#    "-l", "identity-at-volume-center.tfm",
#    "navigator002.nrrd",
#    "navigator001.nrrd",
#    "rraligned-001-on-002.nrrd",
#    "trans-001-to-002.tfm"
#    ])
#
#subprocess.run( dockerprefix +  
#  [ "crl/crkit", "crlResampler",
#    "--dontOrientGeometryImage",
#    "navigator001.nrrd",
#    "trans-001-to-002.tfm",
#    "navigator002.nrrd",
#    "bspline",
#     "vol-001-on-002.nrrd"
#    ])
#
#

#subprocess.run( ["/home/warfield/links/gitrepos/github/sms-mi-reg/pyhelpers/crl-checkerboard-images.py", "--input_image1", "navigator-c-001.nrrd",
#   "--input_image2", "vol-002-on-001.nrrd",
#   "--output_image",  "checker-002-on-001.nrrd" ])


exit(0)


#!/usr/bin/env python3

# Run with python3

# Example usage:
# run-registration-example.py  --refvolume volume-name.nii \
#  --inputtransformfile identity.tfm --outputtransformfile align \
# --slicetarget /data/slices-000.nii /data/slices-015.nii \
#    /data/slices-030.nii  /data/slices-045.nii

#  python3 run-registration-example.py --refvolume func-bold_task-rest480_run-01slimmon_volume_0001_20240328T184812.nrrd --inputtransformfile identity.tfm --outputtransformfile align --slicetarget func-bold_task-rest480_run-01slimmon_volume_0415_20240328T185235.nrrd

import os
import subprocess
import tempfile
import sys
import logging
import glob

import argparse

# Parse the arguments.
parser = argparse.ArgumentParser(description='align a volume to slices.')
parser.add_argument("--refvolume", required=True)
parser.add_argument("--slicetarget", nargs='+', required=True )
parser.add_argument("--outputtransformfile", required=True)
parser.add_argument("--inputtransformfile", required=False)
args = parser.parse_args()

print('args.slicetarget is ' ,  str(args.slicetarget) )

refVolumeName = args.refvolume
print('refVolumeName is ', str(refVolumeName) )

print('Initializing logging.')
# Initialize logging and set its level
logging.basicConfig()
log = logging.getLogger()
log.setLevel( logging.INFO )
log.setLevel( logging.DEBUG )

# host computer : container directory alias
dirmapping = os.getcwd() + ":" + "/data"
dockerprefix = ["docker","run","--rm", "-it", "--init", "-v", dirmapping,
    "--user", str(os.getuid())+":"+str(os.getgid())]
print(dockerprefix)
log.info("docker prefix is " + str(dockerprefix) )

if args.inputtransformfile is None:
    inputTransformFileName = "identity_centered.tfm"
    subprocess.run( dockerprefix +  [ "crl/sms-mi-reg", "crl-identity-transform-at-volume-center.py",
    "--refvolume", args.refvolume,
    "--transformfile", inputTransformFileName ] )
else:
    inputTransformFileName = args.inputtransformfile

print("input transform file is :", inputTransformFileName)

#if args.inputtransformfile is None:
#  inputTransformFileName = "identity.tfm"
#  subprocess.run( dockerprefix + ["crl/crkit", "/opt/crkit/bin/crlIdentityAffineTransform", inputTransformFileName ] )
#  print('Created ' + inputTransformFileName )
#else:
#  inputTransformFileName = args.inputtransformfile
#


print('Docker container command is:')
print(dockerprefix )
print(
[ "crl/sms-mi-reg", "sms-mi-reg", args.refvolume, inputTransformFileName, 
    args.outputtransformfile ]
)

# Run an example registration:
subprocess.run( dockerprefix + ["crl/sms-mi-reg", "sms-mi-reg",
    args.refvolume,
    inputTransformFileName,
    args.outputtransformfile ] + args.slicetarget
)

#subprocess.run(dockerprefix + ["crl/crkit", 
#    "crlResampler", 
#    "-d", 
#    DATADIR+"reference-volume.nii",
#    "sliceTransform"+"pad"+".tfm", 
#    DATADIR+"z_slice-15_vol-0.nii", 
#    "bspline", 
#    "volume-on-slice-15.nii" 
#])
#

# # Re-align volume to reference
# subprocess.run(dockerprefix + ["crl/crkit",
#    "crlResampler",
#    "-d",
#    args.refvolume,
#    "sliceTransform"+args.outputtransformfile+".tfm",
#    args.slicetarget[0],
#    "bspline",
#    "vol2vol_aligned"+args.outputtransformfile+".nrrd"
# ])


exit(0)


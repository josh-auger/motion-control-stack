#!/usr/bin/env python3

# Run with python3

# Example usage:
# run-registration-example.py  --refvolume volume-name.nii \
#  --inputtransformfile identity.tfm --outputtransformfile align \
# --slicetarget /data/slices-000.nii /data/slices-015.nii \
#    /data/slices-030.nii  /data/slices-045.nii  

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
parser.add_argument("--optimizer", required=False)
parser.add_argument("--maxiterations", required=False)
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

use_valgrind = False
if use_valgrind:
    # valgrind_command = ["valgrind", "--leak-check=full"]
    # valgrind_command = ["valgrind", "--tool=callgrind", "--dump-instr=yes", "--collect-jumps=yes", "--callgrind-out-file=/data/callgrind.out"]  # CPU profiling
    valgrind_command = ["valgrind", "--tool=massif", "--massif-out-file=/data/massif.out"]  # Heap memory profiling
else:
    valgrind_command = []

if args.inputtransformfile is None:
    inputTransformFileName = "identity-centered.tfm"
    subprocess.run( dockerprefix +  [ "crl/sms-mi-reg", "crl-identity-transform-at-volume-center.py", 
    "--refvolume", args.refvolume, 
    "--transformfile", inputTransformFileName ] )
else:
    inputTransformFileName = args.inputtransformfile

print("input transform file is :", inputTransformFileName)


print([ "crl/sms-mi-reg", "sms-mi-reg", args.refvolume, inputTransformFileName, args.outputtransformfile ])

# Run an example registration:
# Build base run command
cmd = dockerprefix + ["crl/sms-mi-reg"] + valgrind_command + ["sms-mi-reg",
    args.refvolume,
    inputTransformFileName,
    args.outputtransformfile
] + args.slicetarget
# Add optimizer argument if provided
if args.optimizer is not None:
    cmd += ["--optimizer", args.optimizer]
# Add maxiterations argument if provided
if args.maxiterations is not None:
    cmd += ["--maxiter", str(args.maxiterations)]

print("Full run command :", cmd)
subprocess.run(cmd, check=True)

# subprocess.run( dockerprefix + ["crl/sms-mi-reg"] + valgrind_command + ["sms-mi-reg",
#     args.refvolume,
#     inputTransformFileName,
#     args.outputtransformfile ] + args.slicetarget
# )

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
#subprocess.run(dockerprefix + ["crl/crkit", 
#    "crlResampler", 
#    "-d", 
#    DATADIR+"reference-volume.nii",
#    "sliceTransform"+"pad"+".tfm", 
#    DATADIR+"z_slice-29_vol-0.nii", 
#    "bspline", 
#    "volume-on-slice-29.nii" 
#])


exit(0)


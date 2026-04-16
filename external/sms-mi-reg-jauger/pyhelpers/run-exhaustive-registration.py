#!/usr/bin/env python3

# Title: run-exhaustive-registration.py

# Description:
# Iteratively run registration on two volumes while incrementing through a range of values on two input transform
# parameters. The specified input transform is updated with new parameter values with each registration iteration.
#
# Created on: August 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

# Run with python3

# Example usage:
# python3 run-exhaustive-registration.py --refvolume func-bold_task-rest480_run-01slimmon_volume_0001_20240503T173138.nrrd --slicetarget func-bold_task-rest480_run-01slimmon_volume_0415_20240503T174915.nrrd --outputtransformfile _MImap --inputtransformfile ./rest480_searchtransform_centered_versorrigid3d_0001.tfm


import os
import subprocess
import tempfile
import sys
import logging
import glob
import SimpleITK as sitk
import numpy as np

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

# Specify upper and lower bounds of exhaustive metric search
upper_bound = 1
lower_bound = -1
# Specify search increment
search_increment = 0.05

inputTransformFileName = args.inputtransformfile

transform = sitk.ReadTransform(inputTransformFileName)
trans_x_array = np.arange(lower_bound, upper_bound + search_increment, search_increment)
trans_y_array = np.arange(lower_bound, upper_bound + search_increment, search_increment)
trans_z_array = np.arange(lower_bound, upper_bound + search_increment, search_increment)
rot_x_array = np.arange(lower_bound, upper_bound + search_increment, search_increment)
rot_y_array = np.arange(lower_bound, upper_bound + search_increment, search_increment)

for rot_x in rot_x_array:
    for rot_y in rot_y_array:
        new_params = [float(rot_x), float(rot_y), 0.0, 0.0, 0.0, 0.0]
        log.info(f"new parameters : {new_params}")
        log.info(f"center of rotation : {transform.GetFixedParameters()}")
        transform.SetParameters(new_params)
        sitk.WriteTransform(transform, inputTransformFileName)

        output_file_string = f"{args.outputtransformfile}_xrot_{rot_x:.2f}_yrot_{rot_y:.2f}_"

        print('Docker container command is:')
        print(dockerprefix)
        print(["crl/sms-mi-reg", "sms-mi-reg", args.refvolume, inputTransformFileName, output_file_string])
        subprocess.run(dockerprefix + ["crl/sms-mi-reg", "sms-mi-reg",
                                       args.refvolume,
                                       inputTransformFileName,
                                       output_file_string] + args.slicetarget)

print(f"X rotation array : {rot_x_array}")
print(f"Y rotation array : {rot_y_array}")

exit(0)

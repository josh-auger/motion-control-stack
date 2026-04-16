
# Title: compute_displacement.py

# Description:
# Calculate the displacement value given the input transform parameters, following Tisdall et al. 2012.
# Input transform is assumed to be a VersorRigid3D transform file of six parameters where the first three parameters
# are rotation angles in radians (about x, y, and z-axis, respectively, and the second three parameters are x, y, and
# z-axis translations in millimeters.

# Created on: May 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital


import SimpleITK as sitk
import numpy as np
import argparse

# Compute the displacement:
def compute_displacement(transform_file, radius=50):
    # radius = 50
    versortransform = sitk.ReadTransform(transform_file)
    params = np.asarray( versortransform.GetParameters() )
    print("Transform parameters (VersorRigid3D) : ", params)

    theta = np.abs(np.arccos(0.5 * (-1 + np.cos(params[0]) * np.cos(params[1]) + \
                                    np.cos(params[0]) * np.cos(params[2]) + \
                                    np.cos(params[1]) * np.cos(params[2]) + \
                                    np.sin(params[0]) * np.sin(params[1]) * np.sin(params[2]))))
    drot = radius * np.sqrt((1 - np.cos(theta)) ** 2 + np.sin(theta) ** 2)
    dtrans = np.linalg.norm(params[3:])
    displacement = drot + dtrans

    print("Displacement (mm) : ", displacement)

    return displacement

def main():
    parser = argparse.ArgumentParser(description='Compute displacement from a transform file.')
    parser.add_argument('transformfile', type=str, help='The path to the transform file.')
    args = parser.parse_args()

    transform_file_path = args.transformfile
    displacement = compute_displacement(transform_file_path)


if __name__ == "__main__":
    main()
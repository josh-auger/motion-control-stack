#!/usr/bin/env python3

# Title: run-exhaustive-registration.py

# Description:
# Run an exhaustive search of the registration similarity metric (e.g. mutual information) for a given range of
# input transform parameters. Iterate through all permutations of parameter sets and run registration for each set.
#
# Created on: August 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

# Run with python3

# Example usage:
# python3 run-exhaustive-registration.py --refvolume /path/to/refvolume.nii --slicetarget /path/to/targetvolume.nii --outputstring _VVR_search --inputtransformfile /path/to/inputtransform.tfm


import os
import argparse
import logging
import subprocess
import numpy as np
import pandas as pd
import SimpleITK as sitk
import shutil
from datetime import datetime


def setup_logging():
    """Set up logging to both terminal and a log .txt file with a timestamp."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"exhaustive_search_log_{timestamp}.txt"

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(message)s',
                        handlers=[
                            logging.StreamHandler(),
                            logging.FileHandler(log_filename, mode='w')
                        ])
    log = logging.getLogger()
    log.info(f"Log initialized: {log_filename}")
    return log, timestamp

def validate_file(filepath, description):
    """Ensure a file exists."""
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"{description} not found: {filepath}")
    return filepath

def get_gaussian_sampled_parameters(mean_params, sigma_rot_deg, sigma_trans_mm, rng):
    """Generate 6 transform parameters (3 rot, 3 trans) from a random sampling of a Gaussian distribution."""
    # sample 3 rotations, in degrees
    rot_deg = rng.normal(loc=0.0, scale=sigma_rot_deg, size=3)
    # convert rotations to radians
    rot_rad = np.deg2rad(rot_deg)

    # sample 3 translations, in mm
    trans = rng.normal(loc=0.0, scale=sigma_trans_mm, size=3)

    new_params = [
        mean_params[0] + rot_rad[0],
        mean_params[1] + rot_rad[1],
        mean_params[2] + rot_rad[2],
        mean_params[3] + trans[0],
        mean_params[4] + trans[1],
        mean_params[5] + trans[2],
    ]
    return new_params, rot_deg, trans

def create_identity_transform(ref_volume, input_transform_file, docker_prefix):
    """Create an identity transform at reference volume center using crl/sms-mi-reg tools."""
    log.info("No input transform file provided. Creating an identity transform.")
    command = docker_prefix + [
        "crl/sms-mi-reg", "crl-identity-transform-at-volume-center.py",
        "--refvolume", ref_volume,
        "--transformfile", input_transform_file
    ]
    subprocess.run(command, check=True)
    log.info(f"Identity transform created : {input_transform_file}")

def construct_transform_object(transform, transformtype):
    rotation_matrix = np.array(transform.GetMatrix()).reshape(3, 3)
    translation = transform.GetTranslation()
    center = transform.GetCenter()

    if transformtype == "Affine":
        new_transform = sitk.AffineTransform(3)
        new_transform.SetCenter(center)
        new_transform.SetMatrix(rotation_matrix.ravel())
        new_transform.SetTranslation(translation)
    elif transformtype == "VersorRigid3D":
        new_transform = sitk.VersorRigid3DTransform()
        new_transform.SetCenter(center)
        new_transform.SetMatrix(rotation_matrix.ravel())
        new_transform.SetTranslation(translation)
    elif transformtype == "Euler3D":
        new_transform = sitk.Euler3DTransform()
        new_transform.SetCenter(center)
        new_transform.SetMatrix(rotation_matrix.ravel())
        new_transform.SetTranslation(translation)
    else:
        raise ValueError(f"Unsupported transform type: {transformtype}")
    return new_transform

def run_smsmireg_docker(docker_prefix, ref_volume, transform_file, output_file, slicetargets):
    """Run Docker command for sms-mi-reg."""
    command = (docker_prefix +
               ["crl/sms-mi-reg:initstep90_ftol0.0001",
                "sms-mi-reg",
                ref_volume,
                transform_file,
                output_file]
               + slicetargets +
               ["--optimizer", "LN_SBPLX",
                "--maxiter", "1000"]) # set maximum iterations = 1 to evaluate only the initialization transform
    subprocess.run(command, check=True)
    # log.debug(f"Docker command executed successfully : {command}")

def append_csv_to_master(master_df, counter, euler_params, specific_csv_filepath):
    """Keep a cumulative CSV file of all evaluated transform parameters and calculated MI values."""
    df = pd.read_csv(specific_csv_filepath)
    rename_map = {
        "EvaluatedTransformParams": "VersorTransformParams",
    }
    df = df.rename(columns=rename_map)

    # Construct new entry dict for master search trace
    entry = {
        "counter": counter,
        "RotX_Euler": euler_params[0],
        "RotY_Euler": euler_params[1],
        "RotZ_Euler": euler_params[2],
        "TransX": euler_params[3],
        "TransY": euler_params[4],
        "TransZ": euler_params[5],
    }

    # Add initial values from this search iteration
    first_row = df.iloc[0]
    if first_row is not None:
        entry["VersorTransformParams"] = first_row.get("VersorTransformParams", None)
        entry["FixedParams"] = first_row.get("FixedParams", None)
        entry["MutualInfo_initial"] = first_row.get("MutualInfo", None)
    else:
        entry["VersorTransformParams"] = None
        entry["MutualInfo_initial"] = None

    # Find evaluation row that contains the maximum MI value from this search case
    if "MutualInfo" in df.columns:
        best_row = df.loc[df["MutualInfo"].idxmax()]
    else:
        best_row = None

    if best_row is not None:
        entry["MutualInfo_final"] = best_row.get("MutualInfo", None)
        entry["VersorTransformParams_final"] = best_row.get("VersorTransformParams", None)
    else:
        entry["MutualInfo_final"] = None
        entry["VersorTransformParams_final"] = None

    # Add number of function evaluations to find optimum from this search iteration
    entry['n_evals'] = len(df)
    # Add runtime breakdown
    entry['total_init_time'] = df['InitRuntime'].sum()
    entry['total_parallel_time'] = df['ParallelRuntime'].sum()
    entry['total_serial_time'] = df['SerialRuntime'].sum()
    entry['total_reg_time'] = entry['total_init_time'] + entry['total_parallel_time'] + entry['total_serial_time']

    # Append new entry to current master dataframe
    new_entry_df = pd.DataFrame([entry])
    if master_df is None:
        return new_entry_df

    return pd.concat([master_df, new_entry_df], ignore_index=True)


def main(args):
    # Initialize logging
    global log
    log, timestamp = setup_logging()
    log.info("Starting exhaustive alignment exploration.")

    # Validate reference volume
    ref_volume = validate_file(args.refvolume, "Reference volume")

    # Docker prefix
    dirmapping = os.getcwd() + ":" + "/data"
    docker_prefix = ["docker", "run", "--rm", "-it", "--init", "-v", dirmapping, "--user", f"{os.getuid()}:{os.getgid()}"]
    log.info(f"Docker prefix : {docker_prefix}")

    # Handle input transform file
    if args.inputtransformfile:
        input_transform_file = validate_file(args.inputtransformfile, "Input transform file")
    else:
        input_transform_file = "./identity-centered.tfm"
        create_identity_transform(ref_volume, input_transform_file, docker_prefix)

    # Read in the input transform and convert to Euler3D
    input_transform = sitk.ReadTransform(input_transform_file)
    temp_transform_euler = construct_transform_object(input_transform, "Euler3D")
    temp_transform_file = "./temp_versor_transform.tfm"

    # Read the transform
    input_params = list(temp_transform_euler.GetParameters())
    if len(input_params) != 6:
        raise RuntimeError("Expected 6 parameters in transform. Check transform type (Affine vs Versor3dRigid)!")
    log.info(f"Starting 'mean' transform parameters (rads, mm) : {input_params}")

    # Exhaustive search parameters
    search_quota = 10000
    range_rot_deg = 10
    range_trans_mm = 10
    random_seed = 42
    rng = np.random.default_rng(random_seed)
    sigma_rot_deg = range_rot_deg / 3
    sigma_trans_mm = range_trans_mm / 3
    log.info(f"Search quota : {search_quota} evaluations")
    log.info(f"\tRotation search range (3*sigma) : +/-{range_rot_deg} deg")
    log.info(f"\tTranslation search range (3*sigma) : +/-{range_trans_mm} mm")
    master_df = None
    counter = 1
    for i in range(search_quota):
        if counter==1:
            new_params = input_params
            rot_offset = trans_offset = np.zeros(3)
        else:
            new_params, rot_offset, trans_offset = get_gaussian_sampled_parameters(mean_params=input_params,
                                                                                   sigma_rot_deg=sigma_rot_deg,
                                                                                   sigma_trans_mm=sigma_trans_mm,
                                                                                   rng=rng)
        log.info(f"Iteration {counter:04d}")
        log.info(f"\tRotation offsets (deg) : {rot_offset}")
        log.info(f"\tTranslation offsets (mm) : {trans_offset}")
        log.info(f"\tNew Euler params (rads, mm) : {new_params}")

        temp_transform_euler.SetParameters(new_params)
        temp_transform_versor = construct_transform_object(temp_transform_euler,"VersorRigid3D")
        log.info(f"\tNew Versor params : {temp_transform_versor.GetParameters()}")
        sitk.WriteTransform(temp_transform_versor, temp_transform_file)

        output_file = f"{args.outputstring}_{counter:04d}"
        try:
            run_smsmireg_docker(docker_prefix,
                                ref_volume,
                                temp_transform_file,
                                output_file,
                                args.slicetarget)
        except subprocess.CalledProcessError as e:
            log.error(f"Docker subprocess command failed : {e}")
            counter += 1
            continue

        master_df = append_csv_to_master(master_df, counter, new_params, f"regTrace_{output_file}.csv")
        # os.remove(f"./regTrace_{output_file}.csv")
        # os.remove(f"./alignTransform_{output_file}.tfm")
        counter += 1

    # Clean up temporary transform file
    os.remove(temp_transform_file)
    log.info(f"Temporary transform file removed : {temp_transform_file}")

    # Save master search trace file
    master_csv_path = f"./regTrace_combined_{args.outputstring}_{timestamp}.csv"
    if master_df is not None:
        master_df.to_csv(master_csv_path, index=False)
        log.info(f"Master search trace written to {master_csv_path}")
    else:
        log.warning("No successful registrations – master dataframe is empty.")

if __name__ == "__main__":
    # Parse the arguments.
    parser = argparse.ArgumentParser(description='Align a volume to slices.')
    parser.add_argument("--refvolume", required=True, help="Path to the reference volume file.")
    parser.add_argument("--slicetarget", nargs='+', required=True, help="Paths to slice target files.")
    parser.add_argument("--outputstring", required=True, help="Output file string for all outputs.")
    parser.add_argument("--inputtransformfile", required=False, help="Initial input transform file (optional).")
    args = parser.parse_args()
    main(args)

#!/usr/bin/env python3

# Title: run-registration-restart-from-failures.py

# Description:
# Re-run registration on failed cases from the prior exhaustive random registration experiment. Given a master table of
# all searched cases, identify those that failed to converge on the known solution, and re-run image registration
# starting from the final alignment transform of the first run.
#
# Created on: March 2026
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

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

def parse_transform_str(s):
    """Convert comma-separated string of transform params to a float list."""
    s_clean = str(s).strip("[]")
    return [float(x) for x in s_clean.split(" ")]

def get_failed_cases(master_csv, threshold=3.0):
    """Return dataframe of failed registration cases."""
    df = pd.read_csv(master_csv)
    failed_df = df[df["MutualInfo_final"] < threshold].copy()
    log.info(f"Total cases in CSV : {len(df)}")
    log.info(f"Failed cases (MI < {threshold}) : {len(failed_df)}")
    return failed_df

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
               ["crl/sms-mi-reg:initstep50_ftol0.0001",
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
        "MutualInfo": "MutualInfo_initial",
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
    for col in [
        "VersorTransformParams",
        "FixedParams",
        "MutualInfo_initial",
    ]:
        if col in first_row:
            entry[col] = first_row[col]

    # Find evaluation row that contains the maximum MI value from this search case
    if "MutualInfo_initial" in df.columns:
        best_row = df.loc[df["MutualInfo_initial"].idxmax()]
    else:
        best_row = None

    if best_row is not None:
        entry["MutualInfo_final"] = best_row.get("MutualInfo_initial", None)
        entry["VersorTransformParams_final"] = best_row.get(
            "VersorTransformParams", None
        )
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
    log.info("Starting restart-from-failure experiment.")

    # Validate reference volume
    ref_volume = validate_file(args.refvolume, "Reference volume")

    # Get failed cases from first run of exhaustive registration experiment
    mi_threshold = 1.70
    failed_cases = get_failed_cases(args.mastercsv, threshold=mi_threshold)
    if len(failed_cases) == 0:
        log.warning("No failed cases found.")
        return

    # Docker prefix
    dirmapping = os.getcwd() + ":" + "/data"
    docker_prefix = ["docker", "run", "--rm", "-it", "--init", "-v", dirmapping, "--user", f"{os.getuid()}:{os.getgid()}"]
    log.info(f"Docker prefix : {docker_prefix}")

    # Loop through all failed cases and restart-from-failure
    temp_transform_file = "./temp_versor_transform.tfm"
    master_df = None
    counter = 1
    for idx, row in failed_cases.iterrows():
        original_case_idx = int(row["counter"])
        versor_params = parse_transform_str(row["VersorTransformParams_final"])
        fixed_params = parse_transform_str(row["FixedParams"])
        log.info(f"Iteration {counter:04d}")

        temp_transform_versor = sitk.VersorRigid3DTransform()
        temp_transform_versor.SetParameters(versor_params)
        temp_transform_versor.SetFixedParameters(fixed_params)
        log.info(f"Retest of failed case {original_case_idx}")
        log.info(f"\tVersor params: {versor_params}")
        log.info(f"\tFixed params: {fixed_params}")
        sitk.WriteTransform(temp_transform_versor, temp_transform_file)

        output_file = f"{args.outputstring}_{(original_case_idx - 1):04d}_repeat"
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

        temp_transform_euler = construct_transform_object(temp_transform_versor,"Euler3D")
        new_params = list(temp_transform_euler.GetParameters())
        master_df = append_csv_to_master(master_df, original_case_idx, new_params, f"regTrace_{output_file}.csv")
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
    parser.add_argument("--mastercsv", required=True, help="Master CSV file from the first exhaustive experiment.")
    args = parser.parse_args()
    main(args)

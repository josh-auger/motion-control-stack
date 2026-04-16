# Title: reconstruct_master_trace_from_exhaustive_search.py

# Description:
# Construct the master trace table (.csv) from a directory of traces of single registration instances (.csv). In the
# event that the exhaustive registration experiment terminates prematurely, the master trace will not be exported. This
# helper script will collate the individual registration traces in the same manner.
#

# Created on: April 2026
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import os
import glob
import argparse
import pandas as pd
import numpy as np
import SimpleITK as sitk


# Transform Conversion
def versor_to_euler(versor_params, fixed_params):
    """
    Convert VersorRigid3D parameters + fixed params into Euler3D parameters.
    """
    versor_params = list(map(float, versor_params))
    fixed_params = list(map(float, fixed_params))

    if len(versor_params) != 6:
        raise ValueError("VersorTransformParams must have 6 elements")

    if len(fixed_params) != 3:
        raise ValueError("FixedParams must have 3 elements")

    # Build Versor transform
    versor_tx = sitk.VersorRigid3DTransform()
    versor_tx.SetParameters(versor_params)
    versor_tx.SetFixedParameters(fixed_params)

    # Convert to Euler3D transform
    euler_tx = sitk.Euler3DTransform()
    euler_tx.SetCenter(versor_tx.GetCenter())
    euler_tx.SetMatrix(versor_tx.GetMatrix())
    euler_tx.SetTranslation(versor_tx.GetTranslation())
    return list(euler_tx.GetParameters())  # [rx, ry, rz, tx, ty, tz]


# Parse parameter values from string entry
def parse_parameter_str(s):
    """Convert comma-separated string of transform params to a float list."""
    s_clean = str(s).strip("[]")
    return [float(x) for x in s_clean.split(" ")]

# Single CSV Processing
def process_trace_file(counter, csv_path):
    """
    Extract summary metrics from a single registration trace CSV.
    """
    df = pd.read_csv(csv_path)
    rename_map = {
        "EvaluatedTransformParams": "VersorTransformParams",
        "MutualInfo": "MutualInfo_initial",
    }
    df = df.rename(columns=rename_map)

    if len(df) == 0:
        print(f"Warning: Empty CSV skipped -> {csv_path}")
        return None

    first_row = df.iloc[0]

    entry = {"counter": counter,}

    # Parse initial transform parameters
    versor_params = None
    fixed_params = None
    try:
        if "VersorTransformParams" in first_row:
            versor_params = parse_parameter_str(first_row["VersorTransformParams"])
        if "FixedParams" in first_row:
            fixed_params = parse_parameter_str(first_row["FixedParams"])
    except Exception as e:
        print(f"Warning: Failed parsing params in {csv_path}: {e}")

    # Convert to Euler parameters
    if versor_params is not None and fixed_params is not None:
        try:
            euler_params = versor_to_euler(versor_params, fixed_params)
            entry["RotX_Euler"] = euler_params[0]
            entry["RotY_Euler"] = euler_params[1]
            entry["RotZ_Euler"] = euler_params[2]
            entry["TransX"] = euler_params[3]
            entry["TransY"] = euler_params[4]
            entry["TransZ"] = euler_params[5]
        except Exception as e:
            print(f"Warning: Euler conversion failed for {csv_path}: {e}")

    # Add initial transform params and MI to entry
    entry["VersorTransformParams"] = first_row["VersorTransformParams"]
    entry["FixedParams"] = first_row["FixedParams"]
    entry["MutualInfo_initial"] = first_row.get("MutualInfo_initial", None)

    # Find best (max MI)
    if "MutualInfo_initial" in df.columns and not df["MutualInfo_initial"].isna().all():
        best_row = df.loc[df["MutualInfo_initial"].idxmax()]
        entry["MutualInfo_final"] = best_row.get("MutualInfo_initial", None)
        entry["VersorTransformParams_final"] = best_row.get("VersorTransformParams", None)
    else:
        entry["MutualInfo_final"] = None
        entry["VersorTransformParams_final"] = None

    # Evaluation count
    entry["n_evals"] = len(df)

    # Runtime aggregation
    entry["total_init_time"] = df["InitRuntime"].sum() if "InitRuntime" in df else 0
    entry["total_parallel_time"] = df["ParallelRuntime"].sum() if "ParallelRuntime" in df else 0
    entry["total_serial_time"] = df["SerialRuntime"].sum() if "SerialRuntime" in df else 0
    entry["total_reg_time"] = (entry["total_init_time"] + entry["total_parallel_time"] + entry["total_serial_time"])
    return entry

def extract_counter_from_filename(filepath):
    """
    Extract numeric counter from filename like: regTrace_output_0001.csv
    """
    base = os.path.basename(filepath)
    try:
        return int(base.split("_")[-1].split(".")[0])
    except Exception:
        return None


def rebuild_master_csv(input_dir, output_csv):
    """
    Scan directory for trace CSVs and rebuild master table.
    """
    csv_files = glob.glob(os.path.join(input_dir, "regTrace_*.csv"))
    if not csv_files:
        raise RuntimeError(f"No regTrace CSV files found in: {input_dir}")

    # Sort by extracted counter if possible
    csv_files.sort(key=lambda x: extract_counter_from_filename(x) or 0)

    master_entries = []

    for idx, csv_path in enumerate(csv_files, start=1):
        counter = extract_counter_from_filename(csv_path) or idx
        print(f"Processing {counter:04d}: {os.path.basename(csv_path)}")
        entry = process_trace_file(counter, csv_path)
        if entry is not None:
            master_entries.append(entry)

    if not master_entries:
        raise RuntimeError("No valid entries processed.")

    master_df = pd.DataFrame(master_entries)
    master_df.to_csv(output_csv, index=False)

    print("\n======================================")
    print(f"Master CSV written to: {output_csv}")
    print(f"Total entries: {len(master_df)}")
    print("======================================")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rebuild master CSV from registration trace files.")
    parser.add_argument("--inputdir",required=True,help="Directory containing regTrace_*.csv files",)
    parser.add_argument("--outputcsv",default="regTrace_rebuilt.csv",help="Output master CSV filename",)
    args = parser.parse_args()

    rebuild_master_csv(args.inputdir, args.outputcsv)
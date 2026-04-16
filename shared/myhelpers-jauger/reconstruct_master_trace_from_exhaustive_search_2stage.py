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
from collections import defaultdict


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
def process_trace_file(counter, csv_path1, csv_path2):
    """
    Extract summary metrics from a single registration trace CSV.
    """
    # Read in stage 1 registration traces
    df_stage1 = pd.read_csv(csv_path1)
    rename_map = {
        "EvaluatedTransformParams": "VersorTransformParams",
    }
    df_stage1 = df_stage1.rename(columns=rename_map)

    # Read in stage 2 registration trace
    df_stage2 = pd.read_csv(csv_path2)
    rename_map = {
        "EvaluatedTransformParams": "VersorTransformParams",
    }
    df_stage2 = df_stage2.rename(columns=rename_map)

    first_row = df_stage1.iloc[0]
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
        print(f"Warning: Failed parsing params in {csv_path1}: {e}")

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
            print(f"Warning: Euler conversion failed for {csv_path1}: {e}")

    # Add initial transform params and MI to entry
    entry["VersorTransformParams"] = first_row["VersorTransformParams"]
    entry["FixedParams"] = first_row["FixedParams"]
    entry["MutualInfo_initial"] = first_row.get("MutualInfo", None)

    # STAGE 1: Find evaluation row that contains the maximum MI value from this search case
    if "MutualInfo" in df_stage1.columns:
        best_row_stage1 = df_stage1.loc[df_stage1["MutualInfo"].idxmax()]
    else:
        best_row_stage1 = None

    if best_row_stage1 is not None:
        entry["VersorTransformParams_stage1"] = best_row_stage1.get("VersorTransformParams", None)
        entry["MutualInfo_stage1"] = best_row_stage1.get("MutualInfo", None)
    else:
        entry["VersorTransformParams_stage1"] = None
        entry["MutualInfo_stage1"] = None

    # STAGE 2: Find evaluation row that contains the maximum MI value from this search case
    if "MutualInfo" in df_stage2.columns:
        best_row_stage2 = df_stage2.loc[df_stage2["MutualInfo"].idxmax()]
    else:
        best_row_stage2 = None

    if best_row_stage2 is not None:
        entry["VersorTransformParams_final"] = best_row_stage2.get("VersorTransformParams", None)
        entry["MutualInfo_final"] = best_row_stage2.get("MutualInfo", None)
    else:
        entry["VersorTransformParams_final"] = None
        entry["MutualInfo_final"] = None

    # STAGE 1: Add number of function evaluations and runtimes
    entry['n_evals1'] = len(df_stage1)
    entry['total_init_time1'] = df_stage1['InitRuntime'].sum()
    entry['total_parallel_time1'] = df_stage1['ParallelRuntime'].sum()
    entry['total_serial_time1'] = df_stage1['SerialRuntime'].sum()
    entry['total_reg_time1'] = entry['total_init_time1'] + entry['total_parallel_time1'] + entry['total_serial_time1']

    # STAGE 2: Add number of function evaluations and runtimes
    entry['n_evals2'] = len(df_stage2)
    entry['total_init_time2'] = df_stage2['InitRuntime'].sum()
    entry['total_parallel_time2'] = df_stage2['ParallelRuntime'].sum()
    entry['total_serial_time2'] = df_stage2['SerialRuntime'].sum()
    entry['total_reg_time2'] = entry['total_init_time2'] + entry['total_parallel_time2'] + entry['total_serial_time2']

    # Add total overall registration runtime
    entry['total_reg_time_all'] = entry['total_reg_time1'] + entry['total_reg_time2']
    return entry


def extract_counter_from_filename(filepath):
    """
    Extract counter from filename like: regTrace_output_0001_1.csv
    """
    base = os.path.basename(filepath)
    try:
        parts = base.replace(".csv", "").split("_")
        return int(parts[-2])  # second-to-last element
    except Exception:
        return None


def extract_stage_from_filename(filepath):
    """
    Extract stage (1 or 2) from filename.
    """
    base = os.path.basename(filepath)
    try:
        return int(base.replace(".csv", "").split("_")[-1])
    except Exception:
        return None


def group_csvs_by_counter(csv_files):
    """
    Returns: dict[counter] = {1: filepath_stage1, 2: filepath_stage2}
    """
    grouped = defaultdict(dict)
    for path in csv_files:
        counter = extract_counter_from_filename(path)
        stage = extract_stage_from_filename(path)
        if counter is None or stage not in [1, 2]:
            continue

        grouped[counter][stage] = path
    return grouped


def rebuild_master_csv(input_dir, output_csv):
    """
    Scan directory for trace CSVs and rebuild master table.
    """
    csv_files = glob.glob(os.path.join(input_dir, "regTrace_*.csv"))
    if not csv_files:
        raise RuntimeError(f"No regTrace CSV files found in: {input_dir}")

    grouped = group_csvs_by_counter(csv_files)
    if not grouped:
        raise RuntimeError("No valid grouped CSV pairs found.")

    master_entries = []

    for counter in sorted(grouped.keys()):
        files = grouped[counter]
        csv_path1 = files.get(1)
        csv_path2 = files.get(2)
        if csv_path1 is None or csv_path2 is None:
            print(f"Warning: Missing stage file(s) for counter {counter:04d}, skipping.")
            continue

        print(f"Processing {counter:04d}")
        print(f"\tStage 1 : {os.path.basename(csv_path1)}")
        print(f"\tStage 2 : {os.path.basename(csv_path2)}")

        entry = process_trace_file(counter, csv_path1, csv_path2)
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
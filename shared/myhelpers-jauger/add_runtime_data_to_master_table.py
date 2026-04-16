import os
import pandas as pd
import numpy as np

# === User settings ===
csv_dir = "/home/jauger/GitHubRepos/python-fire-server-jauger/myhelpers-jauger/testRandomSearch_20260305_from20260123_scan001_retroVVR_0000_to_0110/regTraces_scan001_20260123_0000_to_0000_restartFromFailure/"  # Directory with per-instance CSV files
master_csv_path = "/home/jauger/GitHubRepos/python-fire-server-jauger/myhelpers-jauger/testRandomSearch_20260305_from20260123_scan001_retroVVR_0000_to_0110/regTrace_combined_scan001_20260123_0000_to_0000_retestRandomSearch_20260312_105636.csv"

# === Step 1: Load master CSV ===
master_df = pd.read_csv(master_csv_path)
num_master_rows = len(master_df)

# === Step 2: List and sort CSV files ===
csv_files = sorted([f for f in os.listdir(csv_dir) if f.endswith(".csv")])
if len(csv_files) != num_master_rows:
    raise ValueError(f"Number of CSV files ({len(csv_files)}) does not match master CSV rows ({num_master_rows})")

# === Step 3: Compute measures for each CSV file ===
n_evals_list = []
total_init_time_list = []
total_parallel_time_list = []
total_serial_time_list = []
total_reg_time_list = []

for csv_file in csv_files:
    csv_path = os.path.join(csv_dir, csv_file)
    df = pd.read_csv(csv_path)

    n_evals = len(df)
    total_init_time = df['InitRuntime'].sum()
    total_parallel_time = df['ParallelRuntime'].sum()
    total_serial_time = df['SerialRuntime'].sum()
    total_reg_time = total_init_time + total_parallel_time + total_serial_time

    n_evals_list.append(n_evals)
    total_init_time_list.append(total_init_time)
    total_parallel_time_list.append(total_parallel_time)
    total_serial_time_list.append(total_serial_time)
    total_reg_time_list.append(total_reg_time)

# === Step 4: Add new columns to master CSV ===
master_df['n_evals'] = n_evals_list
master_df['total_init_time'] = total_init_time_list
master_df['total_parallel_time'] = total_parallel_time_list
master_df['total_serial_time'] = total_serial_time_list
master_df['total_reg_time'] = total_reg_time_list

# === Step 5: Save updated master CSV ===
updated_master_path = master_csv_path.replace(".csv", "_with_timings.csv")
master_df.to_csv(updated_master_path, index=False)
print(f"Updated master CSV saved to: {updated_master_path}")
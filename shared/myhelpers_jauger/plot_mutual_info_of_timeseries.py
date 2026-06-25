"""
Title: plot_mutual_info_of_timeseries.py

Description:
Reads all *.csv files in a given input directory, sorting by filename.
Extracts the first and max value from the 'MutualInfo' column, and plots both values versus instance in the timeseries.

Usage:
    python plot_max_mutual_info.py --input_dir /path/to/csv_directory

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Date of creation: October 14, 2025
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_rel

def extract_first_and_max_mutual_info(csv_file):
    """Extract the first and maximum MutualInfo values from a CSV file."""
    try:
        df = pd.read_csv(csv_file)
        if "MutualInfo" not in df.columns:
            print(f"Warning: 'MutualInfo' column not found in {csv_file}")
            return None, None
        first_val = df["MutualInfo"].iloc[0]
        max_val = df["MutualInfo"].max()
        return first_val, max_val
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return None, None


def plot_mutual_info(input_dir):
    """Find all CSVs, extract first and max MI values, and plot them."""
    csv_files = sorted(
        [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".csv")]
    )

    if not csv_files:
        print(f"No CSV files found in {input_dir}")
        return

    first_values = []
    max_values = []
    filenames = []
    for i, csv_file in enumerate(csv_files):
        first_val, max_val = extract_first_and_max_mutual_info(csv_file)
        if first_val is not None and max_val is not None:
            first_values.append(first_val)
            max_values.append(max_val)
            filenames.append(os.path.basename(csv_file))
            print(f"[{i}] {os.path.basename(csv_file)} -> "
                  f"First MI = {first_val:.6f}, Max MI = {max_val:.6f}")

    if not max_values:
        print("No valid MutualInfo values extracted.")
        return None, None

    # convert to numpy arrays
    first_values = np.array(first_values)
    max_values = np.array(max_values)

    # save MI values to csv
    summary_path = os.path.join(input_dir, "timeseries_mutualinformation.csv")
    summary_df = pd.DataFrame({
        "Filename": filenames,
        "First_MI": first_values,
        "Max_MI": max_values
    })
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSaved timeseries MI data to: {summary_path}")

    # Plot both first and max MI values
    plt.figure(figsize=(8, 5))
    instances = range(1, len(max_values) + 1)
    plt.plot(instances, first_values, marker='o', linestyle='--', label='Initial MI')
    plt.plot(instances, max_values, marker='s', linestyle='-', label='Final MI')

    plt.title("Mutual Information vs Instance")
    plt.xlabel("Instance (Time Series Index)")
    plt.ylabel("Mutual Information")
    plt.ylim(0.75, 3.00)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # Save the plot as PNG
    plot_path = os.path.join(input_dir, "timeseries_mutualinformation.png")
    plt.savefig(plot_path, dpi=300)
    plt.show()
    print(f"Saved timeseries MI plot to: {plot_path}")

    # --- Compute Statistics ---
    print("\n=== Mutual Information Statistics ===")
    print(f"Mean (std-dev) initial MI value = {np.mean(first_values):.4f} ({np.std(first_values):.4f})")
    print(f"Mean (std-dev) max MI value     = {np.mean(max_values):.4f} ({np.std(max_values):.4f})")

    if len(first_values) > 9:
        first_trimmed = first_values[9:]
        max_trimmed = max_values[9:]
        print(f"\n(Excluding first 9 instances)")
        print(f"Mean (std-dev) initial MI value = {np.mean(first_trimmed):.4f} ({np.std(first_trimmed):.4f})")
        print(f"Mean (std-dev) max MI value     = {np.mean(max_trimmed):.4f} ({np.std(max_trimmed):.4f})")
    else:
        print("\nNot enough data points to exclude the first 9 instances.")

    return first_values, max_values


def run_paired_ttest(first_values, max_values):
    """Run paired t-tests between first and max MI values (full and trimmed)."""
    if len(first_values) == 0 or len(max_values) == 0:
        print("Insufficient data for t-test.")
        return

    print("\n=== Paired t-test (First MI vs Max MI) ===")
    t_stat, p_val = ttest_rel(first_values, max_values, nan_policy='omit')
    print(f"t-statistic = {t_stat:.4f}, p-value = {p_val:.9f}")
    if p_val < 0.05:
        print("→ Statistically significant difference (p < 0.05)")
    else:
        print("→ No statistically significant difference (p ≥ 0.05)")


def main():
    parser = argparse.ArgumentParser(description="Plot maximum MutualInfo from CSV directory.")
    parser.add_argument("--inputdir", required=True, help="Path to directory containing CSV files.")
    args = parser.parse_args()

    if not os.path.isdir(args.inputdir):
        print(f"Error: Directory '{args.inputdir}' not found.")
        return

    first_values, max_values = plot_mutual_info(args.inputdir)
    run_paired_ttest(first_values, max_values)


if __name__ == "__main__":
    main()
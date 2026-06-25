# Title: combine_regtraces.py

# Description:
# Concatenate multiple .csv registration traces within a single directory into one .csv table.
#

# Created on: August 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import os
import pandas as pd
import argparse


def combine_csv_files(input_directory, output_filepath):
    # List to hold the dataframes
    df_list = []

    # Iterate over all files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith('.csv'):
            file_path = os.path.join(input_directory, filename)
            # Read the CSV file
            df = pd.read_csv(file_path)
            # Append the dataframe to the list
            df_list.append(df)

    # Concatenate all dataframes in the list
    combined_df = pd.concat(df_list, ignore_index=True)

    # Write the combined dataframe to the output CSV file
    combined_df.to_csv(output_filepath, index=False)


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Combine multiple CSV files into one master CSV file.")
    parser.add_argument('--inputdir', type=str, help='The directory containing the input CSV files.', required=True)
    parser.add_argument('--outputfilename', type=str, help='The path to the output master CSV file.', required=True)
    args = parser.parse_args()

    outputfilepath = os.path.join(args.inputdir, args.outputfilename)

    # Combine CSV files using the provided arguments
    combine_csv_files(args.inputdir, outputfilepath)


#!/bin/bash

# RUN INSTRUCTIONS:
# Copy this bash script (run-vol2vol-registration.sh) into the parent directory that contains a sub-folder of image volume files
# cd to the location of the bash script
# Execute the following run command: bash run-vol2vol-registration.sh ./input_sub-directory ./input_sub-directory/reference_volume.nrrd

# Note: if you do not specify a reference volume file, the script will default to the first volume file in the sub-directory

# Create log directory and log file
log_dir="./logs"
mkdir -p "$log_dir"
current_date_time=$(date +"%Y%m%d_%H%M%S")
log_file="$log_dir/vol2vol_registration_terminal_log_$current_date_time.txt"
exec > >(tee -a "$log_file") 2>&1
echo "Executing vol-to-vol registration."
echo "Current time stamp: $current_date_time"

# Parse command line input arguments
input_dir=$1
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi
echo "Input directory of target volumes: $input_dir"

if [ -z "$2" ]; then
    ref_volume=$(ls "$input_dir"/* | head -n 1)
    echo "No reference volume specified. Using first file in input directory as reference: $ref_volume"
else
    ref_volume=$2
    if [ ! -f "$ref_volume" ]; then
      echo "Error: Reference volume '$ref_volume' does not exist."
      exit 1
    fi
fi
echo "Input reference volume: $ref_volume"
echo

# Create identity transform about the reference volume center
echo "Creating identity transform at reference volume center"
docker run --rm -it -v "$(pwd)":/data crl/sms-mi-reg /usr/src/moco/pyhelpers/identity-transform-at-volume-center.py --refvolume $ref_volume --transformfile "./identity_centered.tfm"

# Run sms-mi-reg on each volume file in the input directory
for input_file in "$input_dir"/*; do
    if [ -f "$input_file" ]; then
        echo
        echo "Input target volume: $input_file"

        # Extract volume number for output filenames
        numbers=$(echo "$input_file" | grep -oP '(?<=-volume-)\d{4}')
        echo "Output file identifier: _$numbers"

        # Check for a prior transform file
        prior_numbers=$(printf "%04d" $(($(echo "$numbers" | sed 's/^0*//') - 1)))
        prior_transform_file="sliceTransform_${prior_numbers}.tfm"
        if [ -f "./$prior_transform_file" ]; then
            initial_transform="./$prior_transform_file"
            echo "Found prior transform file: $initial_transform"
        else
            initial_transform="./identity_centered.tfm"
            echo "No prior transform file found! Using default input transform: $initial_transform"
        fi

        # Call sms-mi-reg program with selected inputs
        echo
        echo "Running sms-mi-reg vol-to-vol registration."
        docker run --rm -it --init -v `pwd`:/data --user $(id -u):$(id -g) crl/sms-mi-reg sms-mi-reg $ref_volume "$initial_transform" "_$numbers" "$input_file"
        echo "Completed sms-mi-reg vol-to-vol registration."
        echo
        echo

        sleep 1s
    else
        echo "Skipping non-file entry: $input_file"
    fi
done

echo "Registration process completed. Log saved to $log_file"

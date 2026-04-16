#!/bin/bash

# RUN INSTRUCTIONS:
# Copy this bash script (run-slice2vol-registration.sh) into the parent directory that contains a sub-folder of slice image files
# cd to the location of the bash script
# Execute the following run command: bash run-slice2vol-registration.sh ./input_sub-directory_to_slices ./different_sub-directory/reference_volume.nii ./different_sub-directory/metadata.json


# Create log directory and log file
log_dir="./logs"
mkdir -p "$log_dir"
current_date_time=$(date +"%Y%m%d_%H%M%S")
log_file="$log_dir/slice2vol_registration_terminal_log_$current_date_time.txt"
exec > >(tee -a "$log_file") 2>&1
echo "Executing slice-to-vol registration."
echo "Current time stamp: $current_date_time"

# Parse command line input arguments
input_dir=$1
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi
echo "Input directory of target slices: $input_dir"

ref_volume=$2
if [ ! -f "$ref_volume" ]; then
    echo "Error: Reference volume '$ref_volume' does not exist."
    exit 1
fi
echo "Input reference volume: $ref_volume"
echo

json_metadata=$3
if [ ! -f "$json_metadata" ]; then
    echo "Error: JSON metadata file '$json_metadata' does not exist."
    exit 1
fi
echo "JSON metadata file: $json_metadata"
echo

# Create identity transform about the reference volume center
echo "Creating identity transform at reference volume center"
docker run --rm -it -v "$(pwd)":/data crl/sms-mi-reg /usr/src/moco/pyhelpers/crl-identity-transform-at-volume-center.py --refvolume $ref_volume --transformfile "./identity_centered.tfm"
echo

# Read JSON file of scan metadata to pull sms factor, slice timing groups, and number of volumes
echo "Parsing metadata from $json_metadata"
slice_timing_array=$(jq '.SliceTiming' "$json_metadata")
echo "Slice timings: $slice_timing_array"
unique_slice_timings=$(echo "$slice_timing_array" | jq 'unique | sort')
echo "Unique slice timings (acquisition groups): $unique_slice_timings"
sms_factor=$(echo "$slice_timing_array" | jq -r 'group_by(.) | map(length) | max')
echo "SMS factor: $sms_factor"
max_acquisition_number=$(echo "$unique_slice_timings" | jq 'length')
echo "Total acquisitions per volume: $max_acquisition_number"
volume_numbers=$(find "$input_dir" -name "*-volume-*.nii" | grep -oP '(?<=-volume-)\d{4}' | sort -u)
#echo "Detected volume numbers: $volume_numbers"

# Iterate over each volume number
for volume in $volume_numbers; do
    echo
    echo "Processing volume: $volume"

    # Process slice timing groups for this volume
    acquisition_number=1  # Initialize acquisition number for each volume
    for timing in $(echo "$unique_slice_timings" | jq -r '.[]'); do
        echo "Processing acquisition $(printf "%04d" $acquisition_number) (slices timing: $timing) of volume: $volume"

        # Find the indices of slices that match the current timing
        matching_slices=$(echo "$slice_timing_array" | jq --arg timing "$timing" 'map(select(. == ($timing|tonumber))) | to_entries | map(.key)')
        # Convert matching_slices from JSON to a shell-readable list
        slice_indices=$(echo "$matching_slices" | jq -r '.[]')

        # Check if valid slices were found
        if [ -z "$slice_indices" ]; then
            echo "No valid slices found for acquisition $acq_str of volume $volume."
            continue
        fi
        echo "Matching slices for acquisition $acq_str: $slice_indices"

        # Build an array of corresponding slice files for this timing value
        input_slices=()
        for idx in $slice_indices; do
            # Assuming the slice files follow the pattern: *volume-XXXX-slice-XXX.nii
            slice_file=$(printf "$input_dir/*-volume-%04d-slice-%03d.nii" "$volume" "$idx")
            matching_files=($(ls $slice_file 2>/dev/null))
            if [ ${#matching_files[@]} -gt 0 ]; then
                input_slices+=("${matching_files[@]}")
            else
                echo "Warning: No matching slice files found for index $idx of volume $volume. Skipping."
            fi
        done

        # If there are any slices with this timing value, run sms-mi-reg
        if [ ${#input_slices[@]} -gt 0 ]; then
            echo "Running sms-mi-reg for acquisition $acq_str of volume $volume"
            echo "Input slices: ${input_slices[*]}"
            input_slices_str=$(IFS=,; echo "${input_slices[*]}")

            # Check for a prior transform file
            if [ "$acquisition_number" -gt 1 ]; then
                # Check previous acquisition within the same volume
                prior_acq=$(printf "%04d" $(("$acquisition_number" - 1)))
                prior_transform_file="sliceTransform_${volume_number}-${prior_acq}.tfm"
            else
                # If first acquisition, check the last acquisition of the previous volume
                prior_vol=$(printf "%04d" $(("$volume" - 1)))
                prior_acq=$max_acquisition_number
                prior_transform_file="sliceTransform_${prior_vol}-${prior_acq}.tfm"
            fi

            if [ -f "$prior_transform_file" ]; then
                initial_transform="./$prior_transform_file"
                echo "Found prior transform file: $initial_transform"
            else
                initial_transform="./identity_centered.tfm"
                echo "No prior transform file found! Using default input transform: $initial_transform"
            fi

            # Run sms-mi-reg for all slices with this timing value in the current volume
            docker run --rm -it --init -v `pwd`:/data --user $(id -u):$(id -g) crl/sms-mi-reg sms-mi-reg $ref_volume "$initial_transform" "_$volume-$acquisition_number" "$input_slices_str"
            echo "Completed sms-mi-reg for acquisition $acq_str of volume $volume"
            echo
        else
            echo "No valid slices found for acquisition $acq_str of volume $volume."
        fi

        acquisition_number=$((acquisition_number + 1))  # Increment acquisition number for the next timing group
        sleep 1s
    done
done

echo "Registration process completed. Log saved to $log_file"
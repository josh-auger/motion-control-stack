#!/bin/bash

# RUN INSTRUCTIONS:
# Copy this bash script (run-process-enhanced-dicoms.sh) into the parent directory that contains a sub-folder of enhanced dicom files
# cd to the location of the bash script
# Execute the following run command: bash run-process-enhanced-dicoms.sh ./input_sub-directory output_filename_base
#
# Note: if NO output_filename_base is provided in the run command, the default will be "output"

echo "Processing directory of enhanced DICOMs for offline registration."

# Parse command line input argument
input_dir=$1
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi
echo "Input directory of target volumes: $input_dir"

output_filename_base=${2:-output}
echo "Output filename base: $output_filename_base"

# Docker run command to uncompress enhanced dicoms
docker run --rm -u $(id -u):$(id -g) -v "`pwd`":/data crl/dicom-tools \
  uncompress_dicoms.py "./$input_dir" "./$output_filename_base-1_uncompressed"

# Docker run command to sort uncompressed dicoms
docker run --rm -u $(id -u):$(id -g) -v "`pwd`":/data crl/dicom-tools \
  sort_dicoms.py "./$output_filename_base-1_uncompressed" "./$output_filename_base-2_sorted"

# Docker run command to convert sorted uncompressed dicoms into 4D nifti
docker run --rm -u $(id -u):$(id -g) -v "`pwd`":/data crl/dicom-tools \
  dicom_tree_to_nifti.py "./$output_filename_base-2_sorted" "./$output_filename_base-3_converted"

# Docker run command to parse 3D nifti volumes from 4D nifti
input_4d_filepath=$(find "./$output_filename_base-3_converted" -type f -name "*.nii.gz")
if [ -z "$input_4d_filepath" ]; then
    echo "Error: No .nii.gz file found in './$output_filename_base-3_converted' directory."
    exit 1
fi
echo "Found 4D NIFTI file: $input_4d_filepath"
docker run --rm -it -v "$(pwd)":/data crl/sms-mi-reg /usr/src/moco/pyhelpers/crl-extract-volumes-from-fmri.py $input_4d_filepath "$output_filename_base-volume.nii"
output_vol_dir="./$output_filename_base-4_nifti_volumes"
mkdir -p "$output_vol_dir"
mv *.nii "$output_vol_dir"
echo "ITK-compatible volumes were saved to: $output_vol_dir"

# Loop through the generated volumes and extract slices
output_slice_dir="./$output_filename_base-5_nifti_slices"
mkdir -p "$output_slice_dir"

for volume in "$output_vol_dir"/*.nii; do
    volume_basename=$(basename "$volume" .nii)
    echo "Extracting slices from $volume..."
    docker run --rm -it -v "$(pwd)":/data crl/sms-mi-reg /usr/src/moco/pyhelpers/crl-fmri-volume-to-slices.py "$volume" "$volume_basename-slice.nii"
    mv *.nii "$output_slice_dir"
done

echo "ITK-compatible slice images were saved to: $output_slice_dir"
echo "Processing of enhanced DICOMs is complete. Proceed to run offline registration."

exit 0



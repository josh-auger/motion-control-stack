#!/bin/bash

# Example usage: ./shared/myhelpers_jauger/retrieve_from_research_pacs.sh 71451428 20260921

usage() {
  echo "Usage: $0 <MRN> <StudyDate>" >&2
  echo "  StudyDate must use YYYYMMDD format." >&2
}

if [[ $# -ne 2 ]]; then
  usage
  exit 2
fi

MRN="$1"
StudyDate="$2"

if [[ ! ${StudyDate} =~ ^[0-9]{8}$ ]]; then
  echo "Error: StudyDate must be an 8-digit date in YYYYMMDD format." >&2
  usage
  exit 2
fi

# Pull data from Research Synapse:
sudo docker run --rm -it -u $(id -u):$(id -g) \
--volume `pwd`:/data crl/dicom-tools:latest retrieve_dicoms.py \
--outputDir /data/${MRN} \
--aec SYNAPSERESEARCH --aet PACSDCM --namednode 10.20.2.28 --modality MR \
--subjectID ${MRN} --studyDate ${StudyDate}

# Step 1: Uncompress the 4D enhanced DICOMs to separate 3D dicoms
sudo docker run --rm -u $(id -u):$(id -g) -v "`pwd`":/data crl/dicom-tools \
  uncompress_dicoms.py ./${MRN} ./${MRN}-1-uncompressed

# Step 2: Sort the uncompressed 3D DICOMs into each sequence
sudo docker run --rm -u $(id -u):$(id -g) -v "`pwd`":/data crl/dicom-tools \
  sort_dicoms.py ./${MRN}-1-uncompressed ./${MRN}-2-sorted

# Step 3: Convert 3D DICOMs of each sequence into single 4D NIFTI files (w/ JSON metadata file)
sudo docker run --rm -u $(id -u):$(id -g) -v "`pwd`":/data crl/dicom-tools \
  dicom_tree_to_nifti.py ./${MRN}-2-sorted ${MRN}-3-converted

exit 0


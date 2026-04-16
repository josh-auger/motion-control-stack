#!/bin/bash

MRN=4162674
StudyDate=20260123

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



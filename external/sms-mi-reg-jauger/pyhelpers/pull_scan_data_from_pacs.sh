
MRN=5993633
StudyDate=20240725

# Pull the data from Research Synapse:
#sudo docker run --rm -it -u $(id -u):$(id -g) \
#--volume `pwd`:/data crl/dicom-tools:latest retrieve_dicoms.py \
#--outputDir /data/${MRN} \
#--aec SYNAPSERESEARCH --aet PACSDCM --namednode 10.20.2.28 --modality MR \
#--subjectID ${MRN} --studyDate ${StudyDate}

sudo docker run --rm -u $(id -u):$(id -g) -v "`pwd`":/data crl/dicom-tools \
  uncompress_dicoms.py ./${MRN} ./${MRN}-uncompressed
sudo docker run --rm -u $(id -u):$(id -g) -v "`pwd`":/data crl/dicom-tools \
  sort_dicoms.py ./${MRN}-uncompressed ./${MRN}-sorted
sudo docker run --rm -u $(id -u):$(id -g) -v "`pwd`":/data crl/dicom-tools \
  dicom_tree_to_nifti.py ./${MRN}-sorted ${MRN}-converted

exit 0



# python-fire-server
Docker container to monitor a TCP port for scan data communications, parse out image data, and process those images.
Developed for use as part of a multi-container asynchronous processing of FIRE stream scan data for real-time 
volume-to-volume registration.

# Asynchronous multi-container real-time motion monitoring

## General Description
This setup uses two separate containers, running simultaneously, to ensure that two separate processes can run 
independent of each other (i.e. not accidentally block each other from doing their respective job). The two containers 
and their jobs are:

1. "Python-fire-server" is responsible for saving image volumes. It will listen to a TCP port, ingest any sent data, 
compile volumes from the image slices, and save those volume files to a local file directory.

2. "Local-queue-processor" is responsible for executing volume-to-volume registration (VVR). It will monitor a local 
file directory and whenever a new volume file appears, read in that volume data, run VVR (using sms-mi-reg), and save 
the final alignment transform files to a different local file directory.

## Requirements
1. Docker container: jauger/python-fire-server: https://github.com/josh-auger/python-fire-server-jauger
2. Docker container: jauger/local-queue-processor: https://github.com/josh-auger/local-queue-processor

## Build Instructions
1. Git clone both required repositories
2. Within each repo directory, build the respective docker container
    - cd python-fire-server-jauger
      - sh build_fire_server_image.sh
    - cd local-queue-processor
      - sh build_local_queue_processor.sh
3. Check for the necessary script(s) and sub-directoy(s) within each repo directory
    - Within /python-fire-server-jauger:
      - start_fire_server.sh : this bash script will launch the python-fire-server container to receive image data from 
      a TCP port.
      - /received_data sub-directory : this is the sub-directory where all incoming image data files will be saved. A 
      - different output directory can be specified in the start_fire_server.sh bash script.
    - Within /local-queue-processor:
      - start_local_queue_processor.sh : this bash script will launch the local-queue-processor container to monitor a 
      directory for new image data and executing image registration on those files.
      - /working sub-directory : this is the sub-directory where all alignment transforms from registration will be 
      saved. A different output directory can be specified in the start_local_queue_processor.sh bash script.

## Run Instructions
1. In one terminal window, start the python-fire-server docker container
   - cd python-fire-server-jauger
   - sh start_fire_server.sh
     - After a few initialization lines, the terminal window should state:\
     *Starting server and listening for data at 0.0.0.0:9002\
     Serving...*
     - The python-fire-server is now running and ready to receive image data from the specified TCP port. All volume 
     files will be saved to the /received_data sub-directory.
2. In a second terminal window, start the local-queue-processor docker container
   - cd local-queue-processor
   - sh start_local_queue_processor.sh
     - After some initialization lines, the terminal window should state:\
     *Logging will be saved to : /working/log_local_queue_processor_<datestamp>_<timestamp>.log\
     Monitoring directory : /data/*
     - The local-queue-processor is now running and ready to execute image registration on saved volume files. All 
     registration output files will be saved to teh /working sub-directory.
3. Leave both containers running while image data is being communicated through the TCP port. While image data is being 
acquired and sent, both docker containers should be actively logging the receipt and processing of image data.
   - Once image acquisition is complete, proceed with the stop instructions below.

NOTE: These steps are configured for real-time, at-scanner motion measurement and reporting. For offline, retrospective 
motion characterization, please see the "Additional Notes" section below for the "retrospective motion monitoring" 
sub-section. The multi-container setup described here is the same, but requires the additional configuration of 
MRD-formatted image data to be sent to the specified TCP port. Effectively, simulating at-scanner image data acquisition 
for offline replay.

## Stop Instructions
1. At the conclusion of image data acquisition, first check that both containers are done saving/processing data
   - The python-fire-server will receive an MRD_MESSAGE_CLOSE at the conclusion of the sequence, prompting a send 
   MRD_MESSAGE_CLOSE by the container. The container has received all image data available and can be stopped.
   - The local-queue-processor is agnostic to what sequence is being acquired or how much image data is generated, but 
   it will log how many volume files have been processed. If this number corresponds to the total number of volumes 
   acquired in the scan sequence, then processing is complete and the container can be stopped.

IMPORTANT: IF either container is still working (i.e. still saving image data files or still running registrations), 
then let the containers continue running to completion. You can proceed with the next sequence on the scanner side, but 
do NOT send these new data to the TCP port!

2. IF both containers are done, then manually stop each container. The containers will not automatically stop.
   - In each terminal window, press ctrl + c

## Expected Outputs
There are two separate locations where output files get saved:
1. /python-fire-server-jauger/received_data/ : This directory contains all of the saved image data files (.nhdr header 
file for each volume, .raw file for each slice image, and .log logfile of the container).
2. /local-queue-processor/working/ : This directory contains all of the final alignment transforms from running 
registrations (.tfm transform files, .log logfile of the container).


## Additional Notes
### Evaluating motion parameters from registration transforms
The local-queue-processor will save the final alignment transform files from registration into the /working/ 
sub-directory. The motion-monitor container (https://github.com/josh-auger/motion-monitor) can be pointed to the same 
/working sub-directory to evaluate and save the motion characterization data from these transform files. The linked 
github repository includes setup and run instructions.

### Retrospective motion monitoring
The aforementioned multi-container setup remains the same, but offline scan data (in enhanced DICOM format) must be 
converted to MRD format to be sent to the specified TCP port, thus simulating at-scanner image data communications.
1. Converting enhanced DICOM scan data to MRD format: 
https://github.com/josh-auger/python-fire-server-jauger/blob/main/myhelpers-jauger/dicom2mrd_df.py
   - Light revisions likely required to accommodate different scan sequences
2. With the aforementioned python-fire-server and local-queue-processor containers running, open a third terminal 
window to send the offline MRD data
3. cd /python-fire-server-jauger
4. python3 client.py [path-to-MRD-data-file.h5] -a [IP address] -p [TCP port]
   - Specify the IP address of the running python-fire-server container
   - Specify the TCP port to send the MRD data
5. At the conclusion of the simulated image data acquisition, the final alignment transform files can be evaluated by 
the motion-monitor to generate the motion characterization outputs. See "Evaluating motion parameters from registration 
transforms" above.
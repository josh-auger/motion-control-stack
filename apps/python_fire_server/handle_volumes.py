# Title: handle_volumes.py

# Description:
# Receive MRD image data from a TCP port and compile image volumes.
# Using Kelvin's invertcontrast.py code as a base for processing incoming images (https://github.com/kspaceKelvin/python-ismrmrd-server/blob/master/invertcontrast.py)

# Created on: March 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import ismrmrd
import os
import itertools
import logging
import traceback
import numpy as np
import numpy.fft as fft
import xml.dom.minidom
import base64
import ctypes
import re
import mrdhelper
import constants
import h5py
import subprocess
import SimpleITK as sitk
from datetime import datetime
import time
from pprint import pprint
# from connection import Connection
from convert_transform_for_moco import *
import json
from math import floor

class handleVolumes:
    # Initiate an iterator to read each item in the connection
    def __init__(self, connection, datafolder, moco_enabled=False, registration_type=None):
        self.connection = connection
        self.is_exhausted = self.connection.is_exhausted
        self.datafolder = datafolder
        self.imageNo = 0
        self.sliceNo = 0
        self.volcount = 0
        self.groupcount = 0
        self.groupsize = 1
        self.delay = 0.15
        self.framenumber = 0
        self.frameNumberLookup = []
        self.protocol_name = f"protocol_name"
        self.moco_enabled = moco_enabled
        self.registration_type = registration_type

        # Log moco state
        logging.info(f"handleVolumes initialized (moco={'ON' if self.moco_enabled else 'OFF'})")

        try:
            now = datetime.now()
            mrdFilepath = os.path.join(self.datafolder, f"measurement-{now.strftime('%Y%m%dT%H%M%S')}.hdf5")
            self.hf = h5py.File(mrdFilepath,'w')
        except Exception as e:
            logging.exception(e)

    def __iter__(self):
        while not self.is_exhausted:
            yield self.next()

    def __next__(self):
        return self.next()

    def next(self):
        imgGroup = []
        start_time = time.time()

        for item in self.connection:
            if item is None:
                logging.info("\tThe connection will be closed since no data has been received.")
                logging.info(f"\tReceive elapsed time (sec) : {time.time() - start_time}")

                try:
                    logging.info(f"")
                    time.sleep(2)   # sympathetic pause to let queue processor catch up
                    new_subdir_name = self.consolidate_outputs_in_directory()  # Move all output files into a date-time stamped subdirectory
                    self.generate_close_file("closeQ")  # generate dummy close file to pass to local-queue-processor container
                    self.generate_close_file("closeM")  # generate dummy close file to pass to motion-monitor container
                    time.sleep(3)
                    self.consolidate_outputs_in_directory(new_subdir_name=new_subdir_name)  # Another consolidation sweep of final output files
                except Exception as e:
                    logging.exception("Error consolidating output files into subdirectory.")

                # self.connection.send_close()  # no need to send close message. If socket is open, keep listening
                self.is_exhausted = True
                return self.hf

            # Every connection item begins with the MRD message ID, item[0]
            elif item[0] == 1:  # config file, should be first
                self.save_config(item)
                logging.info(f"Received [{item[0]}] config file.")

            elif item[0] == 3:  # series metadata, should be second
                self.save_metadata(item)
                metadata_message, metadata_xml, metadata_bytes = item
                logging.info(f"Received [{item[0]}] metadata.")
                # logging.debug("XML Metadata: %s", metadata_xml)

                self.metadata_object = ismrmrd.xsd.CreateFromDocument(metadata_xml)
                self.protocol_name = self.get_protocol_name_from_metadata()
                self.sms_factor = self.get_smsfactor_from_metadata()
                self.nslices_per_volume = self.get_nslices_per_volume_from_metadata()

                self.get_slice_thickness_from_metadata()
                self.get_spacing_between_slices_from_metadata()

                self.set_groupsize()

            elif item[0] == 1022:   # image, should be all other items
                self.save_image(item)
                image_message, ismrmrd_image, header_bytes, attribute_bytes, data_bytes = item
                logging.info(f"\tReceived [{item[0]}] image number {self.imageNo}")
                img_framenumber = ismrmrd_image.getHead().user_float[7]
                logging.info(f"\tImage header frame number : {img_framenumber}")
                logging.info(f"\tImage position (center) : {[ismrmrd_image.position[x] for x in range(3)]}")

                # save raw image bytes to disk
                img_filename = self.save_raw_image_from_ismrmrd_data(data_bytes)

                img_entry = {'raw_filename': img_filename,
                 'image': ismrmrd_image,
                 'header_bytes': header_bytes,
                 'attribute_bytes': attribute_bytes,
                 'data_bytes': data_bytes}
                imgGroup.append(img_entry)

                # JDA: after first reference volume, all future image slices get individual NRRD headers
                if self.volcount > 0:
                    self.save_nrrd_slice_from_ismrmrd_data(ismrmrd_image, img_filename)
                    self.groupcount += 1
                    logging.info(f"\tVolume {self.volcount} group {floor(self.sliceNo / self.groupsize)} : {self.groupcount}/{self.groupsize}")
                    # JDA: when desired group size is reached (slice, slice group, volume), generate group pointer file
                    if self.groupcount == self.groupsize:
                        recent_group = imgGroup[-self.groupsize:]
                        self.save_group_pointer_file(recent_group)
                        self.groupcount = 0


                self.sliceNo += 1  # increment slice count for each image slice in a volume
                self.imageNo += 1  # increment image count for entire sequence acquisition

                # Check for and send moco feedback after receiving a new image slice
                # Only send motion correction feedback if enabled
                if self.moco_enabled and self.volcount > 0:
                    try:
                        self.package_and_send_moco_feedback()
                    except Exception as e:
                        logging.exception("Error during motion correction feedback.")

                if len(imgGroup) == self.nslices_per_volume:
                    # sort ismrmrd image items by image center z-position (low to high, neg to pos)
                    imgGroup = self.sort_image_slices_by_z_position(imgGroup)

                    if self.volcount == 0:
                        self.series_metadata_dict = self.save_series_metadata_to_json(imgGroup) # includes slice timings

                        # grab the image coordinate frame (direction cosines) from first image of first volume
                        logging.info(f"Establishing image coordinate frame from reference volume :")
                        self.refImgCoordFrame = self.get_image_direction_matrix(imgGroup[0]['image'])
                        self.refImgCoordFrame = self.refImgCoordFrame / np.linalg.norm(self.refImgCoordFrame, axis=0)   # make sure direction cosines are unit vectors
                        logging.info(f"\tRef image coord frame : {np.array2string(self.refImgCoordFrame, precision=4, separator=',', suppress_small=True)}")
                        self.refImgOrigin = self.get_first_voxel_center(imgGroup[0]['image'])

                        self.save_nrrd_volume_from_ismrmrd_data(imgGroup)  # Save as NRRD format, undo moco feedback when writing NHDR header
                    # else:
                    #     self.save_nrrd_volume_from_ismrmrd_data(imgGroup)

                    imgGroup = []   # make sure imgGroup gets cleared after each volume compilation
                    self.sliceNo = 0    # reset slice count to 0 for next volume
                    self.volcount += 1

            elif item[0] == 4:
                logging.info(f"Received [{item[0]}] close message.")
                logging.info(f"Receive elapsed time (sec) : {time.time() - start_time}")

                try:
                    time.sleep(2)
                    new_subdir_name = self.consolidate_outputs_in_directory()  # Move all output files into subdirectory
                    self.generate_close_file("closeQ")  # generate dummy close file to reset local-queue-processor container
                    self.generate_close_file("closeM")  # generate dummy close file to reset motion-monitor container
                    time.sleep(3)
                    self.consolidate_outputs_in_directory(new_subdir_name=new_subdir_name)  # Consolidate outputs from reset processes
                except Exception as e:
                    logging.exception("Error consolidating output files into subdirectory.")

                self.connection.send_close()
                self.is_exhausted = True
                return self.hf



    # ---------- CLASS FUNCTIONS ----------

    def consolidate_outputs_in_directory(self, new_subdir_name=None):
        """
        Move all files (except those with certain file extensions) into a newly created subdirectory. Optionally, the
        output directory can be specified.
        """
        if new_subdir_name is None:
            now = datetime.now()
            new_subdir_name = f"savedData_{now.strftime('%Y%m%dT%H%M%S')}_{self.protocol_name}"

        new_subdir_path = os.path.join(self.datafolder, new_subdir_name)
        os.makedirs(new_subdir_path, exist_ok=True)

        excluded_exts = ('.hdf5','.closeQ','.closeM')  # files to exclude when moving output files to subdirectory

        for fname in os.listdir(self.datafolder):
            fpath = os.path.join(self.datafolder, fname)
            if os.path.isfile(fpath) and not fname.endswith(excluded_exts):
                new_path = os.path.join(new_subdir_path, fname)
                os.rename(fpath, new_path)

        logging.info(f"\tMoved output files into subdirectory : {new_subdir_path}")
        return new_subdir_name


    def log_ismrmrd_image_metadata(self, ismrmrd_img):
        """
        Report all the image metadata fields of the ismrmrd image object
        """
        try:
            meta_dict = ismrmrd.Meta.deserialize(ismrmrd_img.attribute_string)
        except Exception as e:
            logging.error(f"Could not deserialize meta attributes: {e}")
            return

        logging.info("---- ISMRMRD Image Meta Attributes ----")
        for key, value in meta_dict.items():
            logging.info(f"{key}: {value}")
        logging.info("---------------------------------------")

        logging.info("---- ISMRMRD Image Header Attributes ----")
        hdr = ismrmrd_img.getHead()
        for field in dir(hdr):
            if not field.startswith("_"):
                try:
                    logging.info(f"{field}: {getattr(hdr, field)}")
                except Exception:
                    pass
        logging.info("---------------------------------------")


    def get_image_direction_matrix(self, ismrmrd_img):
        """
        Grab the image direction cosines from the first image slice of the reference volume.
        The direction vectors define the image coordinate frame w.r.t. the global device coordinate frame.
        """
        read_dir = [ismrmrd_img.read_dir[x] for x in range(3)]
        phase_dir = [ismrmrd_img.phase_dir[x] for x in range(3)]
        slice_dir = [ismrmrd_img.slice_dir[x] for x in range(3)]
        logging.info(f"\tismrmrd image read_dir : {read_dir}")
        logging.info(f"\tismrmrd image phase_dir : {phase_dir}")
        logging.info(f"\tismrmrd image slice_dir : {slice_dir}")

        imgRowDir = [float(ismrmrd_img.meta["ImageRowDir"][x]) for x in range(3)]
        imgColDir = [float(ismrmrd_img.meta["ImageColumnDir"][x]) for x in range(3)]
        imgSliceDir = [float(ismrmrd_img.meta["ImageSliceNormDir"][x]) for x in range(3)]
        logging.info(f"\tismrmrd image ImageRowDir : {imgRowDir}")
        logging.info(f"\tismrmrd image ImageColumnDir : {imgColDir}")
        logging.info(f"\tismrmrd image ImageSliceNormDir : {imgSliceDir}")
        logging.info(f"\tFlipping z-axis to account for ismrmrd.position in LPS and sorting image slices in ascending order.")
        imgSliceDir = [-x for x in imgSliceDir]
        logging.info(f"\tflipped image ImageSliceDir : {imgSliceDir}")

        # Column stack direction vectors: col 1 = x dir, col 2 = y dir, col 3 = z dir
        D = np.column_stack([imgRowDir, imgColDir, imgSliceDir])

        return D


    def get_latest_transform_file(self):
        """
        Search data directory for all transform files with extension .tfm and grab the most recently modified file.
        """
        VALID_EXTENSIONS = {'.tfm'}
        transform_files = []
        for f in os.listdir(self.datafolder):
            name, ext = os.path.splitext(f)
            if ext not in VALID_EXTENSIONS:
                continue
            transform_files.append(f)

        if not transform_files:
            logging.info("No transform file(s) found...")
            return None, None

        # Sort transform files by modification time
        transform_files.sort(key=lambda f: os.path.getmtime(os.path.join(self.datafolder, f)))
        latest_transform_filename = transform_files[-1]
        latest_transform_filepath = os.path.join(self.datafolder, latest_transform_filename)

        # Check last entry only
        if self.frameNumberLookup:
            last_entry = self.frameNumberLookup[-1]
            if last_entry['regTransformFilename'] == latest_transform_filename:
                return None, None

        # Read and return latest transform
        try:
            latest_transform = sitk.ReadTransform(latest_transform_filepath)
        except Exception as e:
            logging.warning(f"Skipping unreadable transform file {latest_transform_filename}: {e}")
            return None, None
        logging.info(f"Latest transform file to be sent as moco feedback : {latest_transform_filename}")
        logging.info(f"\t{latest_transform}")
        return latest_transform, latest_transform_filename


    def get_moco_transform_from_image_framenumber(self, ismrmrd_image):
        """
        Find the moco transform corresponding to the image header frame number. This is the transform that was applied
        to the imaging FOV during that image slice acquisition.
        """
        # Pull the frame number from the image header
        try:
            img_framenumber = ismrmrd_image.getHead().user_float[7]
            logging.info(f"\timage header frame number : {img_framenumber}")
        except (AttributeError, IndexError) as e:
            logging.warning(f"\tNo frame number found in image header user_float[7]. Returning identity transform.")
            return sitk.AffineTransform(3), sitk.AffineTransform(3)

        # Grab the moco lookup table entry corresponding to the image header frame number
        matched_entry = next(
            (entry for entry in self.frameNumberLookup if entry['frameNumber'] == int(img_framenumber)),
            None
        )
        if matched_entry is None:
            logging.info(f"\t\tNo matching transform found for frame number {img_framenumber}. Returning identity transform.")
            return sitk.AffineTransform(3), sitk.AffineTransform(3)

        moco_transform = matched_entry['mocoTransformObject']
        reg_transform = matched_entry['regTransformObject']
        reg_transform_filename = matched_entry['regTransformFilename']

        logging.info(f"\t\tMoco transform for frame number {img_framenumber} (from {reg_transform_filename}) :")
        logging.info(f"\t\t{moco_transform}")

        logging.info(f"\t\tRegistration transform for frame number {img_framenumber} :")
        logging.info(f"\t\t{reg_transform}")

        return moco_transform, reg_transform


    def update_image_frame_with_moco_transform(self, origin, direction_matrix, transform):
        """
        Update the image header metadata to account for the moco feedback applied by the scanner, placing the image
        data in its true position in global space.
        Moco rotation gets applied about the device isocenter, altering the image origin location and image directions.
        Moco translation gets applied to the image origin.
        Use the alignment transform (in image coord frame) corresponding to the generated moco feedback (in device
        frame) because the image header metadata is written in the image frame.
        """
        logging.info(f"\tUpdating header of volume {self.volcount:04d} to account for applied moco feedback :")
        R = np.array(transform.GetMatrix()).reshape((3, 3))
        t = np.array(transform.GetTranslation())
        # Try to get center
        if hasattr(transform, "GetCenter"):
            c = np.array(transform.GetCenter())
        else:
            fixed_params = transform.GetFixedParameters()
            if len(fixed_params) >= 3:
                c = np.array(fixed_params[:3])
            else:
                c = np.zeros(3)

        # Apply moco rotation (R) and translation (t) to the image header direction cosine matrix and origin, respectively
        # R_inv = np.array((moco_transform.GetInverse()).GetMatrix()).reshape((3,3))   # Get inverse rotation matrix directly from transform object
        R_inv = np.linalg.inv(R)
        updated_direction_matrix = R_inv @ direction_matrix
        updated_origin = (origin - c) @ R + c - t

        logging.info(f"\t\tmoco updated origin = {updated_origin}")
        logging.info(f"\t\tmoco updated dir_x = {updated_direction_matrix[:, 0]}")
        logging.info(f"\t\tmoco updated dir_y = {updated_direction_matrix[:, 1]}")
        logging.info(f"\t\tmoco updated dir_z = {updated_direction_matrix[:, 2]}")

        return updated_origin, updated_direction_matrix


    def package_and_send_moco_feedback(self):
        """
        Convert the most recent registration transform into the device coordinate frame and then package as a moco
        feedback struct and send as a message back to the scanner.
        """
        versor_transform, versor_transform_filename = self.get_latest_transform_file()
        if versor_transform is None:
            logging.info("No new transform file found. Skipping moco feedback.")
            return

        # Convert sms-mi-reg alignment transform to device coordinate frame for MOCO (import functions from convert_transform_for_moco.py)
        moco_transform = convert_transform_for_moco(versor_transform, self.refImgCoordFrame)
        logging.info(f"Moco transform to be sent to scanner :")
        logging.info(f"\t{moco_transform}")

        # Reformat transform for FIRE send and track frame numbers
        self.framenumber += 1
        timestamp = int(time.time())
        moco_struct = package_transform_as_datastruct(moco_transform, timestamp, self.framenumber)

        # Update global frame number lookup table with current frame number, moco data struct, and affine transform object
        entry = {'frameNumber': self.framenumber,
                 'mocoStruct': moco_struct,
                 'mocoTransformObject': moco_transform,
                 'regTransformFilename': versor_transform_filename,
                 'regTransformObject': versor_transform }
        self.frameNumberLookup.append(entry)

        # Write moco feedback entry to cumulative log file
        feedback_log_path = os.path.join(self.datafolder, "log_moco_feedback_sent.log")
        with open(feedback_log_path, 'a') as f:
            f.write(format_moco_struct(moco_struct) + '\n')
            f.write(f"from transform: {versor_transform_filename}\n\n")
        logging.info(f"MOCO return struct : \n{format_moco_struct(moco_struct)}")

        # Send latest moco structure as feedback
        self.connection.send_feedback("MyMocoFeedbackFIRE", moco_struct)
        logging.info(f"Sent MOCO feedback frame number : {self.framenumber} at {timestamp} (ms)")
        return


    def get_smsfactor_from_metadata(self):
        """
        Extract the simultaneous multi-slice acceleration factor from the sequence metadata.
        SMS factor details how many image slices were excited simultaneously (i.e. have the same acquisition time).
        """
        params = getattr(self.metadata_object.userParameters, "userParameterLong", [])
        if not params:
            logging.info("No userParameters found. Returning SMS factor = 1.")
            return 1

        sms_factor = [p.value for p in params if p.name == "MultiBandFactor"]
        if not sms_factor:
            logging.info("\tNo 'MultiBandFactor' field found in series metadata. Using SMS factor = 1.")
            return 1

        logging.info(f"\tSMS factor = {sms_factor[0]}")
        return sms_factor[0]


    def get_nslices_per_volume_from_metadata(self):
        """
        Extract the number of slices per volume from the sequence metadata.
        """
        encoding = self.metadata_object.encoding[0]
        nslices = encoding.encodingLimits.slice.maximum + 1

        logging.info(f"\tnSlices per volume = {nslices}")
        return nslices


    def get_protocol_name_from_metadata(self):
        """
        Extract the protocol name from the sequence metadata.
        """
        protocol_name = self.metadata_object.measurementInformation.protocolName
        protocol_name = protocol_name.replace(" ", "_")  # Replace spaces with underscores

        logging.info(f"\tSequence name = {protocol_name}")
        return protocol_name


    def get_slice_thickness_from_metadata(self):
        sp = getattr(self.metadata_object, "sequenceParameters", None)
        if sp and hasattr(sp, "SliceThickness"):
            logging.info(f"\tSlice thickness = {sp.sliceThickness} mm")
            return sp.sliceThickness

        logging.info(f"\tNo 'SliceThickness' field found in series metadata. Will calculate from first full volume.")
        return None


    def get_spacing_between_slices_from_metadata(self):
        up = getattr(self.metadata_object, "userParameters", None)
        if not up:
            return None

        for p in getattr(up, "userParameterDouble", []):
            if p.name.lower() == "SpacingBetweenSlices":
                logging.info(f"\tSpacing between slices = {p.value} mm")
                return p.value

        logging.info(f"\tNo 'SpacingBetweenSlices' field found in series metadata. Will calculate from first full volume.")
        return None


    def set_groupsize(self):
        """
        Set the group size according to the user-set input flag for registration type. This determines how many image
        slice files are listed in each group pointer file.
        """
        if self.registration_type is not None:
            # Image data collected into groups based on desired registration type
            if self.registration_type == "slice":
                self.groupsize = 1
            elif self.registration_type == "smsgroup":
                self.groupsize = self.sms_factor
                # self.groupsize = 2
            elif self.registration_type == "volume":
                self.groupsize = self.nslices_per_volume
            else:
                logging.info(
                    f"Unknown registration type '{self.registration_type}'. Using default 'slice', group size = 1.")

            logging.info(f"\tRegistration type '{self.registration_type}', image slice group size : {self.groupsize}")
        else:
            logging.info(f"Registration type not specified. Using default 'slice', group size = 1.")
            self.groupsize = 1

        return


    def save_config(self, item):
        config_message, config, config_bytes = item
        if config_message == 1:
            self.hf.create_dataset("Config File", data=bytearray(config_bytes))


    def save_metadata(self, item):
        metadata_message, metadata_xml, metadata_bytes = item

        try:
            metadata = ismrmrd.xsd.CreateFromDocument(metadata_xml)
            if metadata.acquisitionSystemInformation.systemFieldStrength_T is not None:
                logging.info("\tData is from a %s %s at %1.1fT",
                             metadata.acquisitionSystemInformation.systemVendor,
                             metadata.acquisitionSystemInformation.systemModel,
                             metadata.acquisitionSystemInformation.systemFieldStrength_T)
                logging.info("\tIncoming dataset contains %d encoding(s)", len(metadata.encoding))
                logging.info(
                    "\tEncoding type: '%s', FOV: (%s x %s x %s)mm^3, Matrix Size: (%s x %s x %s)",
                    metadata.encoding[0].trajectory,
                    metadata.encoding[0].encodedSpace.matrixSize.x,
                    metadata.encoding[0].encodedSpace.matrixSize.y,
                    metadata.encoding[0].encodedSpace.matrixSize.z,
                    metadata.encoding[0].encodedSpace.fieldOfView_mm.x,
                    metadata.encoding[0].encodedSpace.fieldOfView_mm.y,
                    metadata.encoding[0].encodedSpace.fieldOfView_mm.z)
                self.hf.create_dataset("Metadata XML", data=bytearray(metadata_bytes))
        except:
            logging.warning("Metadata is not a valid MRD XML structure. Passing on metadata as text")


    def save_image(self, item):
        image_message, image, header_bytes, attribute_bytes, data_bytes = item
        self.hf.create_dataset("image_" + str(self.imageNo) + "/header", data=bytearray(header_bytes))
        self.hf.create_dataset("image_" + str(self.imageNo) + "/attribute", data=bytearray(attribute_bytes))
        self.hf.create_dataset("image_" + str(self.imageNo) + "/data", data=bytearray(data_bytes))


    def sort_image_slices_by_z_position(self, item_list):
        """
        Sort list of image slices by z-position of ismrmrd.position (center of acquired image FoV), from lowest (most
        negative) to highest (most positive) value.
        IMPORTANT: ismrmrd.position field is recorded in LPS orientation already! IF ImageSliceNormDir is negative, as
        would be the case this for a head first supine patient orientation in the scanner, then this sort process now
        accounts for the z-axis flip.
        """
        logging.info(f"Sorting image slices by z-position of ismrmrd.position :")
        item_list.sort(key=lambda x: x['image'].position[2])
        first_img = item_list[0]['image']
        last_img = item_list[-1]['image']
        logging.info(f"\tFirst image {item_list[0]['raw_filename']} : Position : {[first_img.position[x] for x in range(3)]}")
        logging.info(f"\tLast image {item_list[-1]['raw_filename']} : Position : {[last_img.position[x] for x in range(3)]}")

        return item_list


    def get_slice_timings(self, item_list):
        """
        Read the acquisition time from every image slice in the input list of images. Convert the acquisition time into
        units of seconds.
        """
        logging.info("Extracting slice acquisition timings of each image slice in the volume :")
        slice_times = []
        for i, item in enumerate(item_list):
            img = item['image']
            acq_time = img.acquisition_time_stamp  # JDA: unsigned 32-bit integer given as "number of samples since midnight"
            slice_times.append(acq_time)
            logging.debug(f"\tSlice {item['raw_filename']} : Acquisition Time = {acq_time}")

        # Compute relative timing in seconds
        PHYSIO_CLOCK_HZ = 400  # JDA: fixed internal clock rate (Hz)
        earliest_time = min(slice_times)
        slice_timing = {
            i: round((t - earliest_time) / PHYSIO_CLOCK_HZ, 4)  # JDA: delta time, in seconds
            for i, t in enumerate(slice_times)
        }

        # Confirm number of slice groups and size of slice groups
        nslices = len(slice_timing)
        nslicegroups = len(np.unique(list(slice_timing.values())))
        nslices_per_group = nslices / nslicegroups
        logging.info(f"\tNumber of slices : {nslices}")
        logging.info(f"\tNumber of slice groups (unique acquisition times) : {nslicegroups}")
        logging.info(f"\tNumber of slices per group : {nslices_per_group}")
        if nslices_per_group > self.sms_factor:
            logging.info(f"\tUpdating SMS factor from {self.sms_factor} to {nslices_per_group}")
            self.sms_factor = int(nslices_per_group)
            logging.info(f"\tResetting group size.")
            self.set_groupsize()

        # Log slice order by acquisition time
        sorted_indices = sorted(range(len(slice_timing)), key=lambda i: slice_timing[i])
        logging.info(f"\tSlice acquisition order : {sorted_indices}")

        return slice_timing


    def save_series_metadata_to_json(self, item_list=None):
        """
        Save series metadata to json file.
        Option: Include slice acquisition timings from image group, if provided as input.
        """
        now = datetime.now()
        metadata_filename = f"{self.protocol_name}_metadata_{now.strftime('%Y%m%dT%H%M%S')}.json"
        metadata_filepath = os.path.join(self.datafolder, metadata_filename)

        # Recursive helper to convert ISMRMRD object to dict
        def to_dict(obj):
            """Recursively convert ismrmrd.xsd objects to dictionaries."""
            if isinstance(obj, list):
                return [to_dict(item) for item in obj]
            elif hasattr(obj, "__dict__"):
                return {
                    k: to_dict(v)
                    for k, v in obj.__dict__.items()
                    if not k.startswith("_")
                }
            else:
                return obj

        # Convert metadata object to a serializable dictionary
        metadata_dict = to_dict(self.metadata_object)

        # Optionally compute and add SliceTiming if item_list is provided
        if item_list is not None:
            slice_timing = self.get_slice_timings(item_list)
            metadata_dict["SliceTiming"] = slice_timing

        # Write formatted JSON
        with open(metadata_filepath, "w") as f:
            json.dump(metadata_dict, f, indent=4, default=str)

        logging.info(f"\tSequence metadata saved as : {metadata_filename}")

        return metadata_dict


    def save_group_pointer_file(self, item_list):
        """
        Save a separate pointer file (.txt) with the filenames of the image slices in the current group.
        Group size is determined by the user-set "registration_type" (slice, slice group, volume).
        """
        groupNo = floor(self.sliceNo / self.groupsize)
        logging.info(f"\t\tCompiling {len(item_list)} images as volume {self.volcount} group {groupNo}")

        # grab image slice filenames from current item list
        slice_filenames = []
        for i, img_item in enumerate(item_list):
            img_filename = item_list[i]['raw_filename']
            img_header_filename = os.path.splitext(img_filename)[0] + ".nhdr"   # get equivalent image slice header
            slice_filenames.append(img_header_filename)

        # generate a .txt file that points to the group of image slices
        group_filename = f"{self.protocol_name}_volume_{self.volcount:04d}_group_{groupNo:04d}.txt"
        group_filepath = os.path.join(self.datafolder, group_filename)

        # write all slice header filenames to group pointer file, one per line
        try:
            with open(group_filepath, 'w') as f:
                for fname in slice_filenames:
                    f.write(fname + '\n')
            logging.info(f"\t\tGroup pointer saved as : {group_filepath}")
        except Exception as e:
            logging.error(f"Failed to save group pointer file: {e}")
            raise

        return


    def save_nrrd_slice_from_ismrmrd_data(self, ismrmrd_img, img_filename):
        """
        Save a single image slice .nhdr header file corresponding to one .raw image slice data.
        Update the image slice header to place the image data in the "true position" in global device coordinate frame.
        """
        header_filename = os.path.splitext(img_filename)[0] + ".nhdr"
        header_filepath = os.path.join(self.datafolder, header_filename)
        logging.info(f"\tImage slice header saved as : {header_filename}")

        # Determine image data type
        image_type = str(ismrmrd.get_dtype_from_data_type(ismrmrd_img.data_type))
        nrrd_type_overrides = {
            "float32": "float",
            "float64": "double"
        }
        if image_type in nrrd_type_overrides:
            logging.warning(
                f"\timage type {image_type} unsupported by NRRD format. Set image type to '{nrrd_type_overrides[image_type]}'.")
            image_type = nrrd_type_overrides[image_type]

        # Size of single slice
        size = [ismrmrd_img.matrix_size[x] for x in range(3)]  # (z, y, x)

        # Compute voxel spacing (x, y, z)
        spacing = [
            float(ismrmrd_img.field_of_view[0]) / size[2],
            float(ismrmrd_img.field_of_view[1]) / size[1],
            float(ismrmrd_img.field_of_view[2]) / size[0],
        ]

        # Compute the slice’s origin in scanner space
        slice_origin = self.get_first_voxel_center(ismrmrd_img)

        # Orientation and scaling
        scaled_direction_matrix = self.refImgCoordFrame * spacing

        # Apply any motion correction if available
        moco_transform, reg_transform = self.get_moco_transform_from_image_framenumber(ismrmrd_img)
        updated_origin, updated_scaled_direction_matrix = self.update_image_frame_with_moco_transform(
            slice_origin, scaled_direction_matrix, reg_transform)

        dir_x = updated_scaled_direction_matrix[:, 0]
        dir_y = updated_scaled_direction_matrix[:, 1]
        dir_z = updated_scaled_direction_matrix[:, 2]

        # Write NRRD header for the single slice
        with open(header_filepath, "w") as hdr:
            hdr.write("NRRD0004\n")
            hdr.write("# Created by slice header writer\n")
            hdr.write(f"type: {image_type}\n")
            hdr.write("dimension: 3\n")
            hdr.write(f"sizes: {size[2]} {size[1]} {size[0]}\n")  # x y z
            hdr.write("space: left-posterior-superior\n")
            hdr.write(f"space origin: ({updated_origin[0]}, {updated_origin[1]}, {updated_origin[2]})\n")
            hdr.write(
                f"space directions: "
                f"({dir_x[0]},{dir_x[1]},{dir_x[2]}) "
                f"({dir_y[0]},{dir_y[1]},{dir_y[2]}) "
                f"({dir_z[0]},{dir_z[1]},{dir_z[2]})\n"
            )
            hdr.write("encoding: raw\n")
            hdr.write("endian: little\n")
            hdr.write(f"data file: {img_filename}\n")

        return


    def save_nrrd_volume_from_ismrmrd_data(self, item_list):
        """
        Save a list of image slices as a single NRRD image volume defined by a detached header file (.nhdr) and the
        corresponding image slice data (.raw).
        """
        logging.info(f"Processing {len(item_list)} images into NRRD volume {self.volcount}")
        slice_filenames = []

        # save each image data_bytes to a .raw file
        for i, img_item in enumerate(item_list):
            img_filename = item_list[i]['raw_filename']
            slice_filenames.append(img_filename)

        # generate a .nhdr file for the volume
        now = datetime.now()
        header_filename = f"{self.protocol_name}_volume_{self.volcount:04d}_{now.strftime('%Y%m%dT%H%M%S')}.nhdr"
        header_filepath = os.path.join(self.datafolder, header_filename)
        self.generate_nrrd_header(item_list[0]['image'], header_filepath, slice_filenames, item_list[1]['image'])

        # generate pointer file to volume
        groupNo = floor(self.sliceNo / self.groupsize)
        group_filename = f"{self.protocol_name}_volume_{self.volcount:04d}_group_{groupNo:04d}.txt"
        group_filepath = os.path.join(self.datafolder, group_filename)
        try:
            with open(group_filepath, 'w') as f:
                f.write(header_filename + '\n')
            logging.info(f"Saved group pointer file: {group_filepath}")
        except Exception as e:
            logging.error(f"Failed to save group pointer file: {e}")
            raise

        return


    def save_raw_image_from_ismrmrd_data(self, data_bytes):
        """
        Save a byte string of image data directly to a .raw file.
        """
        img_filename = f"{self.protocol_name}_volume_{self.volcount:04d}_slice_{self.sliceNo:04d}.raw"
        img_filepath = os.path.join(self.datafolder, img_filename)
        logging.info(f"\tRaw image data saved as : {img_filename}")
        with open(img_filepath, "wb") as f:
            f.write(data_bytes)

        return img_filename


    def get_first_voxel_center(self, ismrmrd_img):
        """
        Calculate the spatial location of the center of the first (top left) voxel of the image bounding box (FoV).
        The ismrmrd image header "position" field is the center of the excited volume in LPS orientation, relative to
        the device isocenter.
        NRRD format requires that "space origin" = center of the first sample in the array (that is, the location of
        the first voxel).
        """
        logging.info(f"\tIdentifying image origin as center of first (top left) voxel :")
        center = np.array(ismrmrd_img.position)
        logging.info(f"\t\tCenter of first image slice (ismrmrd.position) : {center}")
        size = np.array(ismrmrd_img.matrix_size)
        # logging.info(f"\t\tSize of image slice (z y x) : {size}")
        spacing = np.array([
            float(ismrmrd_img.field_of_view[0]) / size[2],
            float(ismrmrd_img.field_of_view[1]) / size[1],
            float(ismrmrd_img.field_of_view[2]) / size[0]])
        # logging.info(f"\t\tFoV (mm) of image slice (x y z) : {[ismrmrd_img.field_of_view[x] for x in range(3)]}")
        # logging.info(f"\t\tSpacing of image slice (x y z) : {spacing}")

        D = self.refImgCoordFrame * spacing
        # logging.info(f"\t\tDirection matrix of image (scaled) : {np.array2string(D, precision=4, separator=',', suppress_small=True)}")

        offset_indices = (size - 1) / 2.0   # indices of center voxel (z y x)
        offset_indices = np.array([offset_indices[2], offset_indices[1], offset_indices[0]])    # indices of center (x y z)
        logging.info(f"\t\tIndices (in pixels) of center voxel (x y z) : {offset_indices}")
        first_voxel_center = center - D @ offset_indices
        logging.info(f"\t\tOrigin as center (in mm) of first voxel (x y z) : {first_voxel_center}")

        return first_voxel_center


    def check_for_slice_gap(self, ismrmrd_img, next_ismrmrd_img):
        """
        Compare image slice FOV thickness with difference in z-position between two image slices to determine if there
        is any slice gap.
        The two image slices must be consecutive in global space (not necessarily consecutive acquisitions).
        Slice gap matters when compiling a single detached header file of multiple image slices (i.e. saving volumes)
        """
        z0 = ismrmrd_img.position[2]
        z1 = next_ismrmrd_img.position[2]
        dz = abs(z1 - z0)   # spacing between images slices
        fov_z = ismrmrd_img.field_of_view[2]    # image slice thickness
        slice_gap = (dz - fov_z)
        logging.info(f"\tMeasuring slice thickness and spacing between slices from full image volume:")
        logging.info(f"\tSlice thickness = {fov_z}")
        logging.info(f"\tSpacing between slices = {dz}")
        if dz > fov_z:
            logging.info(f"\tSlice gap = {slice_gap}")
            return slice_gap
        else:
            logging.info(f"\tNo slice gap to accommodate")
            return None


    def generate_nrrd_header(self, ismrmrd_img, header_filepath, slice_filenames, next_ismrmrd_img):
        """
        Write a detached header file (.nhdr) for the image volume defined by the list of image slice filenames.
        The origin and direction vectors are updated to place the image data in the "true position" in global device
        coordinate space.
        """
        logging.info(f"Volume header file saved as : {header_filepath}")

        # Pull header information from first ismrmrd image item
        image_type = str(ismrmrd.get_dtype_from_data_type(ismrmrd_img.data_type))     # cast NumPy dtype as a string to check for overrides
        nrrd_type_overrides = {
            "float32": "float",
            "float64": "double"
        }
        if image_type in nrrd_type_overrides:
            logging.warning(f"\timage type {image_type} unsupported by NRRD format. Set image type to '{nrrd_type_overrides[image_type]}'.")
            image_type = nrrd_type_overrides[image_type]
        logging.info(f"\timage type : {image_type}")

        size = [ismrmrd_img.matrix_size[x] for x in range(3)]   # for a single slice image = 1 x 90 x 90 (z y x)
        size[0] = len(slice_filenames)  # z dim = num slices per volume, sizes order = z y x
        spacing = [float(ismrmrd_img.field_of_view[0]) / size[2], \
                   float(ismrmrd_img.field_of_view[1]) / size[1], \
                   float(ismrmrd_img.field_of_view[2])]
        logging.info(f"\torigin (mm, x y z) : {self.refImgOrigin}")
        logging.info(f"\tsize (pixels, x y z) : {size[2]} {size[1]} {size[0]}")
        logging.info(f"\tspacing (mm, x y z) : {spacing}")

        # IF spacing between slices > slice thickness (slice gap exists!), then use "spacing between slices" to properly space image slices in volume header
        slice_gap = self.check_for_slice_gap(ismrmrd_img, next_ismrmrd_img)
        if slice_gap is not None:
            spacing[2] += slice_gap
            logging.info(f"\tUpdated spacing with slice gap : {spacing}")

        # update the image coordinate frame with the applied moco feedback to get the true image position/orientation
        scaled_direction_matrix = self.refImgCoordFrame * spacing
        moco_transform, reg_transform = self.get_moco_transform_from_image_framenumber(ismrmrd_img)
        updated_origin, updated_scaled_direction_matrix = self.update_image_frame_with_moco_transform(self.refImgOrigin, scaled_direction_matrix, reg_transform)
        # updated_origin, updated_scaled_direction_matrix = self.refImgOrigin, scaled_direction_matrix  # Dev: keep original image header fields in order to "spoof motion" in the acquired images

        dir_x = updated_scaled_direction_matrix[:,0]    # column 1 = x direction
        dir_y = updated_scaled_direction_matrix[:,1]    # column 2 = y direction
        dir_z = updated_scaled_direction_matrix[:,2]    # column 3 = z direction

        with open(header_filepath, "w") as hdr:
            hdr.write("NRRD0004\n")
            hdr.write("# Created by manual header writer\n")
            hdr.write(f"type: {image_type}\n")
            hdr.write("dimension: 3\n")
            hdr.write(f"sizes: {size[2]} {size[1]} {size[0]}\n")
            hdr.write("space: left-posterior-superior\n")
            hdr.write(f"space origin: ({updated_origin[0]}, {updated_origin[1]}, {updated_origin[2]})\n")
            hdr.write(f"space directions: ({dir_x[0]},{dir_x[1]},{dir_x[2]}) ({dir_y[0]},{dir_y[1]},{dir_y[2]}) ({dir_z[0]},{dir_z[1]},{dir_z[2]})\n")
            hdr.write("encoding: raw\n")
            hdr.write("endian: little\n")
            hdr.write("data file: LIST\n")
            for i in range(size[0]):
                hdr.write(f"{slice_filenames[i]}\n")

        logging.info(f"\n")
        return


    def generate_close_file(self, close_extension):
        """
        Write a dummy file with custom extension (close_extension) to trigger "close" operations in another process.
        """
        now = datetime.now()
        filename = f"{self.protocol_name}_CLOSE_{now.strftime('%Y%m%d_%H%M%S')}.{close_extension}"
        filepath = os.path.join(self.datafolder, filename)
        tmp_filepath = filepath + ".tmp"  # Atomic rename to avoid partial reads
        logging.info(f"\tWriting dummy close file : {filepath}")
        try:
            with open(tmp_filepath, "w", encoding="utf-8") as f:
                f.write(f"# Dummy CLOSE trigger file : Timestamp {now.strftime('%Y%m%d_%H%M%S')}\n")
            os.replace(tmp_filepath, filepath)
        except Exception as e:
            logging.error(f"Failed to save CLOSE trigger file: {e}")    # Clean up temp file if something went wrong
            try:
                if os.path.exists(tmp_filepath):
                    os.remove(tmp_filepath)
            except Exception:
                pass
            raise
        return


# Title: mrd2tuncH5.py

# Description:
# Convert standard MRD files into a modified "Tunc h5" format in which all data is in byte string format in acquisition order.
# The resulting "Tunc h5" file is compatible with Tunc's FIRE send code to replay scan data over a TCP port offline.

# Created on: March 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import logging
import h5py
from datetime import datetime
import ismrmrd
import os

# Use Kelvin's send_image() code in connection.py as reference for parsing out each piece of the MRD object

# Kelvin's standard MRD format saves all images within /dataset/images_0
# Read MRD object
file_path = "/home/jauger/Radiology_Research/SLIMM_data/Saved_scan_data/server_saved_data/func-bold_task-rest960_run-01slimmon_2024-03-14-183816_29.h5"
dataIn = ismrmrd.Dataset(file_path, '/dataset', True)

# Read config file

# Read series metadata
metadataIn = data.read_xml_header()

# loop through all images within /images_0 and grab header, attribute, and data for each image
num_images = data.number_of_images("images_0")

# ----- MRD_MESSAGE_ISMRMRD_IMAGE (1022) -----------------------------------
    # This message contains a single [x y z cha] image.
    # Message consists of:
    #   ID               (   2 bytes, unsigned short)
    #   Fixed header     ( 198 bytes, mixed         )
    #   Attribute length (   8 bytes, uint64_t      )
    #   Attribute data   (  variable, char          )
    #   Image data       (  variable, variable      )

# Convert each item to byte string format
# self.socket.send(image.getHead())
# self.socket.send(constants.MrdMessageAttribLength.pack(len(image.attribute_string)))
# self.socket.send(bytes(image.attribute_string, 'utf-8'))
# self.socket.send(bytes(image.data))
img_head = dataIn.read_image("images_0",400).meta
img_attr = dataIn.read_image("images_0", 400).attribute_string
img_data = dataIn.read_image("images_0", 400).data

# Need to re-order all images according to acquisition order for Tunc h5 format for scan replay
# Store all images in list array and then sort by instance number?


# Use Tunc's saveData.py code to generate the Tunc h5 formatted file
class SaveData:
    def __init__(self, connection, datafolder):
        self.connection = connection
        self.is_exhausted = self.connection.is_exhausted
        self.datafolder = datafolder
        self.imageNo = 0
        try:
            now = datetime.now()
            self.hf = h5py.File(self.datafolder + 'measurement-' + now.strftime("%Y%m%dT%H%M%S") + '.hdf5','w')
            while not os.path.exists(self.datafolder + 'measurement-' + now.strftime("%Y%m%dT%H%M%S") + '.hdf5'):
                logging.info("Waiting OS for file creation...")
        except Exception as e:
            logging.exception(e)

    def __iter__(self):
        while not self.is_exhausted:
            yield self.next()

    def __next__(self):
        return self.next()

    def next(self):
        for item in self.connection:
            # Break out if a connection was established but no data was received
            if item is None:
                logging.info("\tThe connection will be closed since no data has been received.")
                self.connection.send_close()
                self.is_exhausted = True
                return self.hf
            elif item[0] == 1:
                # If message ID = 1, then message is the config (file or text), should be first
                self.save_config(item)
            elif item[0] == 3:
                # IF message ID = 3, then message is the metadata (text), should be second
                self.save_metadata(item)
            elif item[0] == 1022:
                # If message ID = 1022, then message is an image
                # Rest of  messages are images
                self.save_image(item)
            else:
                self.connection.send_close()
                self.is_exhausted = True
                return self.hf

    def save_config(self, item):
        config_message, config, config_bytes = item
        # set config_message = 1
        if config_message == 1:
            self.hf.create_dataset("Config File", data=bytearray(config_bytes))

    def save_metadata(self, item):
        xml_message, metadata_xml, metadata_bytes = item
        # set xml_message = 3
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
        image_message, item, header_bytes, attribute_bytes, data_bytes = item
        # set image_message = 1022
        self.hf.create_dataset("image_" + str(self.imageNo) + "/header", data=bytearray(header_bytes))
        self.hf.create_dataset("image_" + str(self.imageNo) + "/attribute", data=bytearray(attribute_bytes))
        self.hf.create_dataset("image_" + str(self.imageNo) + "/data", data=bytearray(data_bytes))
        self.imageNo += 1
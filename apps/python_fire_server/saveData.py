import logging
import h5py
from datetime import datetime
import ismrmrd
import os

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
                # First message is the config (file or text)
                self.save_config(item)
            elif item[0] == 3:
                # Second messages is the metadata (text)
                self.save_metadata(item)
            elif item[0] == 1022:
                # Rest of the messages are images
                self.save_image(item)
            else:
                self.connection.send_close()
                self.is_exhausted = True
                return self.hf

    def save_config(self, item):
        config_message, config, config_bytes = item
        if config_message == 1:
            self.hf.create_dataset("Config File", data=bytearray(config_bytes))

    def save_metadata(self, item):
        xml_message, metadata_xml, metadata_bytes = item
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
        self.hf.create_dataset("image_" + str(self.imageNo) + "/header", data=bytearray(header_bytes))
        self.hf.create_dataset("image_" + str(self.imageNo) + "/attribute", data=bytearray(attribute_bytes))
        self.hf.create_dataset("image_" + str(self.imageNo) + "/data", data=bytearray(data_bytes))
        self.imageNo += 1
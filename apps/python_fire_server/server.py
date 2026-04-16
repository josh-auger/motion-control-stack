
import constants
import socket
import logging
import ismrmrd.xsd
import importlib
import os
import json
from datetime import datetime

from connection import Connection
from handle_volumes import handleVolumes
from saveData import SaveData

class Server:
    """
    Something something docstring.
    """

    def __init__(self, address, port, savedataFolder, registration_type=None, moco_enabled=False):
        self.savedataFolder = savedataFolder
        self.moco_enabled = moco_enabled
        self.registration_type = registration_type
        self.log_handler = None
        self.log_filename = None
        self.setup_logging()

        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.socket.bind((address, port))
        logging.info(
            "Starting server and listening for data at %s:%d (moco=%s)",
            address, port, "ON" if self.moco_enabled else "OFF"
        )

    def setup_logging(self):
        """Attach (or replace) a FileHandler for logging to a new file."""
        if not os.path.exists(self.savedataFolder):
            os.makedirs(self.savedataFolder)

        self.log_filename = os.path.join(self.savedataFolder,f"python-fire-server-jauger_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

        # Remove old file handler if present
        if self.log_handler:
            logging.getLogger().removeHandler(self.log_handler)
            self.log_handler.close()

        self.log_handler = logging.FileHandler(self.log_filename)
        self.log_handler.setFormatter(logging.Formatter("%(asctime)s - %(message)s"))
        logging.getLogger().addHandler(self.log_handler)
        logging.info(f"New connection — logging reset to: {self.log_filename}")

    def set_moco_enabled(self, enable: bool):
        """Enable/disable moco after initialization."""
        self.moco_enabled = enable
        logging.info("Motion correction flag set to: %s", "ON" if enable else "OFF")

    def serve(self):
        logging.debug("Serving... ")
        self.socket.listen(0)

        while True:
            sock, (remote_addr, remote_port) = self.socket.accept()
            self.setup_logging()    # Always create a new logfile for each new established connection
            logging.info("Accepting connection from: %s:%d", remote_addr, remote_port)
            self.handle(sock)

    def handle(self, sock):
        try:
            connection = Connection(sock)
            outdata = handleVolumes(connection, self.savedataFolder, self.moco_enabled, self.registration_type)    # Run analysis code "handle_volumes.py"
            for item in outdata:
                hfile = item
        except Exception as e:
            logging.exception(e)
        finally:
            # Encapsulate shutdown in a try block because the socket may have already been closed on the other side
            try:
                sock.shutdown(socket.SHUT_RDWR)
            except:
                pass
            sock.close()
            logging.info("\tSocket closed")

            # Dataset may not be closed properly if a close message is not received
            if hfile:  # if hfile was returned by SaveData(), then close the file and end
                try:
                    hfile.close()
                    logging.info("\tIncoming data stream saved to %s", self.savedataFolder)
                    logging.info(f"\n\n")
                    logging.info("---- Python-fire-server reset ----\n")
                except Exception as e:
                    logging.exception(e)
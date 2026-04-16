#!/usr/bin/python3

from server import Server

import argparse
import logging
import sys
import os
import signal
import subprocess
from datetime import datetime

defaults = {
    'host':           '0.0.0.0',
    'port':           9002,
    'savedataFolder': '/tmp/share/saved_data',
    'moco':           'off'  # default to no moco feedback send
}

def main(args):
    # Print motion correction status
    if args.moco.lower() == 'on':
        print("Motion correction ENABLED.")
        moco_enabled = True
    else:
        print("Motion correction DISABLED.")
        moco_enabled = False

    # Create a multi-threaded dispatcher to handle incoming connections
    server = Server(args.host, args.port, args.savedataFolder, args.regtype, moco_enabled=(args.moco.lower() == 'on'))

    # Trap signal interrupts (e.g. ctrl+c, SIGTERM) and gracefully stop
    def handle_signals(signum, frame):
        print("Received signal interrupt -- stopping server")
        server.socket.close()
        sys.exit(0)

    signal.signal(signal.SIGTERM, handle_signals)
    signal.signal(signal.SIGINT,  handle_signals)

    # Start server
    server.serve()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Example server for MRD streaming format',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--port',            type=int,            help='Port')
    parser.add_argument('-H', '--host',            type=str,            help='Host')
    parser.add_argument('-v', '--verbose',         action='store_true', help='Verbose output.')
    parser.add_argument('-S', '--savedataFolder',  type=str,            help='Folder to save incoming data')
    parser.add_argument('-r', '--crlf',            action='store_true', help='Use Windows (CRLF) line endings')
    parser.add_argument('--moco', type=str, choices=['on', 'off'], help='Enable/disable motion correction')
    parser.add_argument('--regtype', type=str, choices=['slice', 'smsgroup', 'volume'], help='Registration type')

    parser.set_defaults(**defaults)
    args = parser.parse_args()

    # Logging format
    fmt = '%(asctime)s - %(message)s\r' if args.crlf else '%(asctime)s - %(message)s'
    # Set up root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    # Remove any pre-existing logger handlers
    for h in root_logger.handlers[:]:
        root_logger.removeHandler(h)
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter(fmt))
    root_logger.addHandler(console_handler)

    if args.savedataFolder:
        print("Saving to folder: ", args.savedataFolder)
        if not os.path.exists(args.savedataFolder):
            os.makedirs(args.savedataFolder)

    main(args)

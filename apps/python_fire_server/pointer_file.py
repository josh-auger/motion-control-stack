"""Helpers for publishing queue pointer files."""

import os
import uuid


def write_pointer_file_atomically(filepath, filenames):
    """Write pointer entries before atomically publishing ``filepath``."""
    temporary_filepath = f"{filepath}.{uuid.uuid4().hex}.tmp"
    descriptor = None

    try:
        descriptor = os.open(
            temporary_filepath,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL,
            0o666,
        )
        pointer_file = os.fdopen(descriptor, "w")
        descriptor = None
        with pointer_file:
            for filename in filenames:
                pointer_file.write(filename + "\n")

        os.replace(temporary_filepath, filepath)
    except Exception:
        if descriptor is not None:
            try:
                os.close(descriptor)
            except OSError:
                pass
        try:
            os.remove(temporary_filepath)
        except OSError:
            pass
        raise

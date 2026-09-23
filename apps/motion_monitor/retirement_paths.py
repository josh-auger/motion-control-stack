"""Shared path construction for motion-monitor retirement outputs."""

from __future__ import annotations

import os
from os import PathLike


def processed_acquisition_directory(
    input_dir: str | PathLike[str],
    protocol_name: str,
) -> str:
    """Return the existing per-protocol processed-directory convention."""
    return os.path.join(
        os.path.abspath(os.fspath(input_dir)),
        f"processed_{protocol_name}",
    )

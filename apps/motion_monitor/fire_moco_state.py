"""Validation for FIRE's atomically published MOCO release state."""

from __future__ import annotations

from dataclasses import dataclass
import json
import logging
import os
import re
from os import PathLike


# JSON content uses a non-.json suffix so existing queue and motion-monitor
# acquisition-file filters cannot mistake coordination state for metadata.
FIRE_MOCO_STATE_FILENAME = "fire_moco_state.moco"
FIRE_MOCO_STATE_SCHEMA_VERSION = 1
_REGISTRATION_TRANSFORM_RE = re.compile(
    r"^alignTransform_(?P<index>\d+)_(?P<volume>\d+)-(?P<group>\d+)\.tfm$"
)


@dataclass(frozen=True)
class FireMocoState:
    committed_registration_index: int
    committed_transform_filename: str
    release_before_registration_index: int


def parse_registration_transform_index(filename: str) -> int | None:
    """Return a strict production registration index, excluding identity files."""
    match = _REGISTRATION_TRANSFORM_RE.fullmatch(filename)
    return int(match.group("index")) if match is not None else None


def load_fire_moco_state(
    input_dir: str | PathLike[str],
) -> FireMocoState:
    """Load and validate one complete FIRE state publication.

    FIRE publishes with ``os.replace``, so no retries are needed: readers see
    either the previous complete state or the next complete state.
    """
    path = os.path.join(os.fspath(input_dir), FIRE_MOCO_STATE_FILENAME)
    with open(path, "r") as stream:
        payload = json.load(stream)

    if not isinstance(payload, dict):
        raise ValueError("FIRE MOCO state is not an object")
    if payload.get("schema_version") != FIRE_MOCO_STATE_SCHEMA_VERSION:
        raise ValueError("unsupported FIRE MOCO state schema")

    committed = payload.get("committed_registration_index")
    release_before = payload.get("release_before_registration_index")
    filename = payload.get("committed_transform_filename")
    if isinstance(committed, bool) or not isinstance(committed, int) or committed < 0:
        raise ValueError("invalid committed_registration_index")
    if (
        isinstance(release_before, bool)
        or not isinstance(release_before, int)
        or release_before != committed
    ):
        raise ValueError("release-before index does not match committed index")
    if not isinstance(filename, str):
        raise ValueError("invalid committed_transform_filename")
    filename_index = parse_registration_transform_index(filename)
    if filename_index != committed:
        raise ValueError("committed filename/index mismatch")

    return FireMocoState(
        committed_registration_index=committed,
        committed_transform_filename=filename,
        release_before_registration_index=release_before,
    )


def read_fire_moco_state(
    input_dir: str | PathLike[str],
    *,
    logger: logging.Logger | None = None,
) -> FireMocoState | None:
    """Return validated state or conservatively report it unavailable."""
    try:
        return load_fire_moco_state(input_dir)
    except (OSError, json.JSONDecodeError, ValueError, TypeError) as error:
        target_logger = logger if logger is not None else logging.getLogger(__name__)
        target_logger.warning("FIRE MOCO state unavailable; retirement delayed: %s", error)
        return None

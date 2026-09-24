"""Atomic publication of FIRE's committed MOCO transform advancement."""

from __future__ import annotations

import json
import logging
import os
import time
import uuid
from os import PathLike

from moco_transform_discovery import parse_registration_transform_filename


FIRE_MOCO_STATE_FILENAME = "fire_moco_state.moco"
FIRE_MOCO_STATE_SCHEMA_VERSION = 1


class FireMocoStatePublisher:
    """Publish FIRE's current committed transform anchor without blocking MOCO."""

    def __init__(
        self,
        directory: str | PathLike[str],
        *,
        logger: logging.Logger | None = None,
    ) -> None:
        self.directory = os.path.abspath(os.fspath(directory))
        self.path = os.path.join(self.directory, FIRE_MOCO_STATE_FILENAME)
        self._last_published_index: int | None = None
        self._logger = logger if logger is not None else logging.getLogger(__name__)

    def publish(
        self,
        *,
        committed_registration_index: int,
        committed_transform_filename: str,
        feedback_frame_number: int,
        trigger_image_identifier: int,
        trigger_volume: int,
        trigger_slice: int,
        trigger_group: int,
    ) -> bool:
        """Atomically publish a non-regressing release-before watermark.

        Publication is best effort. ``False`` means the prior complete state, if
        any, remains authoritative and retirement must stay conservative.
        """
        parsed = parse_registration_transform_filename(committed_transform_filename)
        if (
            parsed is None
            or isinstance(committed_registration_index, bool)
            or parsed.registration_index != committed_registration_index
        ):
            self._logger.warning(
                "FIRE MOCO state not published: invalid committed transform %s "
                "for index %r",
                committed_transform_filename,
                committed_registration_index,
            )
            return False

        # FIRE is the sole writer. Validate an existing publication only until
        # this publisher has successfully established its in-memory watermark;
        # subsequent commits avoid an unnecessary read on every image.
        existing_index = (
            self._existing_committed_index()
            if self._last_published_index is None
            else None
        )
        known_indices = tuple(
            index
            for index in (self._last_published_index, existing_index)
            if index is not None
        )
        newest_known = max(known_indices) if known_indices else None
        if newest_known is not None and committed_registration_index < newest_known:
            self._logger.warning(
                "FIRE MOCO state regression refused: attempted %d, published %d",
                committed_registration_index,
                newest_known,
            )
            return False
        payload = {
            "schema_version": FIRE_MOCO_STATE_SCHEMA_VERSION,
            "committed_registration_index": committed_registration_index,
            "committed_transform_filename": committed_transform_filename,
            "release_before_registration_index": committed_registration_index,
            "committed_transform_volume": parsed.volume,
            "committed_transform_group": parsed.group,
            "feedback_frame_number": feedback_frame_number,
            "trigger_image_identifier": trigger_image_identifier,
            "trigger_volume": trigger_volume,
            "trigger_slice": trigger_slice,
            "trigger_group": trigger_group,
            "published_at_unix_ns": time.time_ns(),
        }
        temporary_path = f"{self.path}.{uuid.uuid4().hex}.tmp"
        descriptor = None
        try:
            descriptor = os.open(
                temporary_path,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                0o666,
            )
            stream = os.fdopen(descriptor, "w")
            descriptor = None
            with stream:
                json.dump(payload, stream, sort_keys=True)
                stream.write("\n")
                stream.flush()
            os.replace(temporary_path, self.path)
        except Exception as error:
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except OSError:
                    pass
            try:
                os.remove(temporary_path)
            except OSError:
                pass
            self._logger.warning(
                "FIRE MOCO state publication failed for index %d: %s",
                committed_registration_index,
                error,
            )
            return False

        self._last_published_index = committed_registration_index
        return True

    def _existing_committed_index(self) -> int | None:
        """Return a valid existing index solely to prevent publication regression."""
        try:
            with open(self.path, "r") as stream:
                payload = json.load(stream)
        except FileNotFoundError:
            return None
        except (OSError, json.JSONDecodeError) as error:
            self._logger.warning(
                "FIRE MOCO state could not validate existing publication; "
                "replacing it from committed in-memory state: %s",
                error,
            )
            return None

        if not isinstance(payload, dict):
            self._logger.warning(
                "FIRE MOCO state existing publication is invalid; replacing it "
                "from committed in-memory state"
            )
            return None

        committed = payload.get("committed_registration_index")
        release_before = payload.get("release_before_registration_index")
        filename = payload.get("committed_transform_filename")
        parsed = (
            parse_registration_transform_filename(filename)
            if isinstance(filename, str)
            else None
        )
        if (
            payload.get("schema_version") != FIRE_MOCO_STATE_SCHEMA_VERSION
            or isinstance(committed, bool)
            or not isinstance(committed, int)
            or committed < 0
            or release_before != committed
            or parsed is None
            or parsed.registration_index != committed
        ):
            self._logger.warning(
                "FIRE MOCO state existing publication is invalid; replacing it "
                "from committed in-memory state"
            )
            return None
        return committed

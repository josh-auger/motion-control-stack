"""Buffered, acquisition-local profiling for FIRE MOCO transform handling."""

from __future__ import annotations

import csv
from datetime import datetime
import logging
import os


PROFILE_FIELDS = (
    "lookup_id",
    "incoming_image_identifier",
    "volume",
    "slice",
    "group",
    "root_entry_count",
    "valid_transform_count",
    "newer_candidate_count",
    "selected_transform_filename",
    "selected_registration_index",
    "last_consumed_registration_index_before",
    "last_consumed_transform_filename_before",
    "last_consumed_registration_index_after",
    "last_consumed_transform_filename_after",
    "discovery_selection_ms",
    "transform_read_ms",
    "conversion_ms",
    "packaging_ms",
    "conversion_package_ms",
    "feedback_log_ms",
    "feedback_send_ms",
    "outcome",
)


class FireMocoProfiler:
    """Write one low-overhead CSV row for each FIRE MOCO lookup."""

    def __init__(self, directory: str) -> None:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        self.path = os.path.join(directory, f"fire_moco_profile_{stamp}.csv")
        self._stream = open(
            self.path,
            "w",
            newline="",
            buffering=1024 * 1024,
        )
        self._writer = csv.DictWriter(
            self._stream,
            fieldnames=PROFILE_FIELDS,
            extrasaction="ignore",
        )
        self._writer.writeheader()
        self._lookup_id = 0
        self._written = 0
        self.disabled = False

    def record(self, row: dict) -> None:
        """Record a lookup without allowing profiling failures to affect MOCO."""
        if self.disabled:
            return
        self._lookup_id += 1
        output = {field: "" for field in PROFILE_FIELDS}
        output.update(row)
        output["lookup_id"] = self._lookup_id
        try:
            self._writer.writerow(output)
            self._written += 1
            if self._written % 256 == 0:
                self._stream.flush()
        except Exception as error:
            self.disabled = True
            logging.warning("FIRE MOCO profiling disabled after write failure: %s", error)

    def close(self) -> None:
        try:
            self._stream.close()
        except Exception:
            pass

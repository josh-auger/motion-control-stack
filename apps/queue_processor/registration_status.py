"""Durable registration outcome tracking for one acquisition."""

from __future__ import annotations

import json
import logging
import os
import tempfile
from collections.abc import Iterable
from os import PathLike


SCHEMA_VERSION = 1
TERMINAL_GROUP_STATUSES = frozenset({"registered", "skipped", "failed"})


class RegistrationStatusTracker:
    """Track terminal group outcomes and atomically publish finalized volumes."""

    def __init__(
        self,
        output_path: str | PathLike[str],
        *,
        logger: logging.Logger | None = None,
    ) -> None:
        self.output_path = os.path.abspath(os.fspath(output_path))
        self._logger = logger if logger is not None else logging.getLogger(__name__)
        self.reset_memory()
        self._load_existing_snapshot()

    def reset_memory(self) -> None:
        """Clear in-memory acquisition state without replacing published history."""
        self._volumes: dict[int, dict[str, object]] = {}
        self._revision = 0

    def set_expected_groups(
        self,
        volume_number: int,
        group_ids: Iterable[int],
    ) -> None:
        """Set the exact registration-group membership for a volume."""
        volume = self._volume_state(volume_number)
        expected_groups = {int(group_id) for group_id in group_ids}
        if not expected_groups or min(expected_groups) < 0:
            raise ValueError("Expected registration groups must be non-negative and nonempty.")

        existing = volume["expected_groups"]
        if existing is not None and existing != expected_groups:
            raise ValueError(
                f"Volume {volume_number} expected groups changed from "
                f"{sorted(existing)} to {sorted(expected_groups)}."
            )
        if volume["registration_finalized"] and existing != expected_groups:
            raise RuntimeError(f"Volume {volume_number} is already finalized.")
        volume["expected_groups"] = expected_groups

    def record_group(
        self,
        volume_number: int,
        group_id: int,
        status: str,
        *,
        transform_filename: str | None = None,
    ) -> None:
        """Record one terminal group result, rejecting changes after finalization."""
        if status not in TERMINAL_GROUP_STATUSES:
            raise ValueError(f"Unsupported registration group status: {status}")
        if status == "registered" and not transform_filename:
            raise ValueError("Registered groups require a transform filename.")
        if status != "registered" and transform_filename is not None:
            raise ValueError("Only registered groups may include a transform filename.")

        volume = self._volume_state(volume_number)
        if volume["registration_finalized"]:
            raise RuntimeError(
                f"Volume {volume_number} is finalized; no later group result is allowed."
            )

        expected_groups = volume["expected_groups"]
        group_id = int(group_id)
        if expected_groups is not None and group_id not in expected_groups:
            raise ValueError(
                f"Group {group_id} is not expected for volume {volume_number}."
            )

        group_record = {"status": status}
        if transform_filename is not None:
            group_record["transform_filename"] = os.path.basename(transform_filename)

        groups = volume["groups"]
        existing = groups.get(group_id)
        if existing is not None and existing != group_record:
            raise ValueError(
                f"Conflicting result for volume {volume_number} group {group_id}: "
                f"{existing} versus {group_record}."
            )
        groups[group_id] = group_record

    def is_finalized(self, volume_number: int) -> bool:
        volume = self._volumes.get(int(volume_number))
        return bool(volume and volume["registration_finalized"])

    def finalize_volume_if_ready(self, volume_number: int) -> bool:
        """Finalize and publish a volume once every expected group is terminal."""
        volume_number = int(volume_number)
        volume = self._volume_state(volume_number)
        if volume["registration_finalized"]:
            return False

        expected_groups = volume["expected_groups"]
        groups = volume["groups"]
        if expected_groups is None or set(groups) != expected_groups:
            return False

        volume["registration_finalized"] = True
        self._revision += 1
        self.publish()

        counts = {
            status: sum(record["status"] == status for record in groups.values())
            for status in TERMINAL_GROUP_STATUSES
        }
        self._logger.info(
            "REG STATUS: volume %d finalized (%d registered, %d skipped, %d failed)",
            volume_number,
            counts["registered"],
            counts["skipped"],
            counts["failed"],
        )
        return True

    def finalize_all_ready(self) -> list[int]:
        """Finalize every complete volume and retry durable publication at close."""
        finalized_now = []
        for volume_number in sorted(self._volumes):
            if self.finalize_volume_if_ready(volume_number):
                finalized_now.append(volume_number)
        if any(volume["registration_finalized"] for volume in self._volumes.values()):
            self.publish()
        return finalized_now

    def unfinalized_volumes(self) -> dict[int, list[int]]:
        """Return missing terminal group IDs for diagnostic close handling."""
        missing: dict[int, list[int]] = {}
        for volume_number, volume in self._volumes.items():
            if volume["registration_finalized"]:
                continue
            expected_groups = volume["expected_groups"]
            groups = volume["groups"]
            if expected_groups is None:
                missing[volume_number] = []
            else:
                missing[volume_number] = sorted(expected_groups - set(groups))
        return missing

    def publish(self) -> None:
        """Atomically publish all finalized volume records."""
        document = self.as_document()
        output_dir = os.path.dirname(self.output_path)
        basename = os.path.basename(self.output_path)
        temporary_path: str | None = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=output_dir,
                prefix=f".{basename}.",
                suffix=".tmp",
                delete=False,
            ) as temporary_file:
                temporary_path = temporary_file.name
                json.dump(document, temporary_file, indent=2, sort_keys=True)
                temporary_file.write("\n")
                temporary_file.flush()
                os.fsync(temporary_file.fileno())
            os.replace(temporary_path, self.output_path)
            temporary_path = None
        finally:
            if temporary_path is not None:
                try:
                    os.remove(temporary_path)
                except FileNotFoundError:
                    pass

    def as_document(self) -> dict[str, object]:
        """Return the compact public representation of finalized history."""
        volumes: dict[str, object] = {}
        for volume_number in sorted(self._volumes):
            volume = self._volumes[volume_number]
            if not volume["registration_finalized"]:
                continue
            groups = volume["groups"]
            volumes[str(volume_number)] = {
                "registration_finalized": True,
                "groups": {
                    str(group_id): groups[group_id]
                    for group_id in sorted(groups)
                },
            }
        return {
            "schema_version": SCHEMA_VERSION,
            "revision": self._revision,
            "volumes": volumes,
        }

    def _volume_state(self, volume_number: int) -> dict[str, object]:
        volume_number = int(volume_number)
        if volume_number < 0:
            raise ValueError("Volume numbers must be non-negative.")
        return self._volumes.setdefault(
            volume_number,
            {
                "expected_groups": None,
                "groups": {},
                "registration_finalized": False,
            },
        )

    def _load_existing_snapshot(self) -> None:
        """Restore finalized history so a queue restart cannot schedule it again."""
        if not os.path.isfile(self.output_path):
            return
        try:
            with open(self.output_path, "r", encoding="utf-8") as status_file:
                document = json.load(status_file)
            if document.get("schema_version") != SCHEMA_VERSION:
                raise ValueError("unsupported or missing schema_version")
            volumes = document.get("volumes")
            if not isinstance(volumes, dict):
                raise ValueError("volumes must be an object")

            restored: dict[int, dict[str, object]] = {}
            for volume_text, public_volume in volumes.items():
                if not str(volume_text).isdigit() or not isinstance(public_volume, dict):
                    raise ValueError("invalid volume record")
                if public_volume.get("registration_finalized") is not True:
                    raise ValueError("snapshot contains a non-finalized volume")
                public_groups = public_volume.get("groups")
                if not isinstance(public_groups, dict) or not public_groups:
                    raise ValueError("finalized volume has no groups")
                groups: dict[int, dict[str, str]] = {}
                for group_text, group_record in public_groups.items():
                    if not str(group_text).isdigit() or not isinstance(group_record, dict):
                        raise ValueError("invalid group record")
                    status = group_record.get("status")
                    if status not in TERMINAL_GROUP_STATUSES:
                        raise ValueError("invalid terminal group status")
                    if status == "registered" and not isinstance(
                        group_record.get("transform_filename"), str
                    ):
                        raise ValueError("registered group lacks transform filename")
                    groups[int(group_text)] = dict(group_record)
                restored[int(volume_text)] = {
                    "expected_groups": set(groups),
                    "groups": groups,
                    "registration_finalized": True,
                }
            self._volumes = restored
            self._revision = int(document.get("revision", 0))
        except Exception as error:
            self._logger.error(
                "REG STATUS: failed to restore existing snapshot %s: %s",
                self.output_path,
                error,
            )

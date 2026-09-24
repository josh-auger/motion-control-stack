"""Registration-index based transform discovery for FIRE MOCO feedback."""

from __future__ import annotations

from dataclasses import dataclass
import os
import re
from os import PathLike


_REGISTRATION_TRANSFORM_RE = re.compile(
    r"^alignTransform_(?P<index>\d+)_(?P<volume>\d+)-(?P<group>\d+)\.tfm$"
)


@dataclass(frozen=True)
class RegistrationTransform:
    """The sequence fields encoded in a production registration filename."""

    filename: str
    registration_index: int
    volume: int
    group: int


@dataclass(frozen=True)
class TransformDiscovery:
    """One root-directory snapshot and its highest eligible transform."""

    root_entry_count: int
    valid_transform_count: int
    newer_candidate_count: int
    selected: RegistrationTransform | None


def parse_registration_transform_filename(
    filename: str,
) -> RegistrationTransform | None:
    """Parse a production transform name, excluding identity and other TFM files."""
    match = _REGISTRATION_TRANSFORM_RE.fullmatch(filename)
    if match is None:
        return None
    return RegistrationTransform(
        filename=filename,
        registration_index=int(match.group("index")),
        volume=int(match.group("volume")),
        group=int(match.group("group")),
    )


def discover_newest_registration_transform(
    directory: str | PathLike[str],
    last_consumed_registration_index: int | None,
) -> TransformDiscovery:
    """Select the highest valid index newer than FIRE's consumed index.

    Discovery intentionally performs one ``os.listdir`` and one pass over its
    names. It does not inspect mtimes, sort candidates, or open transforms.
    """
    filenames = os.listdir(directory)
    valid_count = 0
    newer_count = 0
    selected = None

    for filename in filenames:
        candidate = parse_registration_transform_filename(filename)
        if candidate is None:
            continue
        valid_count += 1
        if (
            last_consumed_registration_index is not None
            and candidate.registration_index <= last_consumed_registration_index
        ):
            continue
        newer_count += 1
        if (
            selected is None
            or candidate.registration_index > selected.registration_index
        ):
            selected = candidate

    return TransformDiscovery(
        root_entry_count=len(filenames),
        valid_transform_count=valid_count,
        newer_candidate_count=newer_count,
        selected=selected,
    )

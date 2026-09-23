"""Archive transforms after motion-monitor consumes their successors."""

from __future__ import annotations

import logging
import os
import re
from os import PathLike

try:
    from .retirement_paths import processed_acquisition_directory
except ImportError:  # Support execution from the motion_monitor application directory.
    from retirement_paths import processed_acquisition_directory


_REGISTRATION_TRANSFORM_RE = re.compile(
    r"^alignTransform_\d+_\d+-\d+\.tfm$"
)
_IDENTITY_TRANSFORM_RE = re.compile(
    r"^alignTransform_\d+_\d+-\d+_identity\.tfm$"
)


def retirement_enabled_for_moco_flag(moco_flag: str | None = None) -> bool:
    """Enable retirement only for the project's explicit MOCO-off mode."""
    value = os.environ.get("MOCO_FLAG", "off") if moco_flag is None else moco_flag
    return value.strip().lower() == "off"


class TransformRetirer:
    """Move a consumed predecessor while retaining identity and current files.

    This is intentionally an uninterrupted-runtime optimization. It does not
    rebuild motion history from the archive after a motion-monitor restart.
    """

    def __init__(
        self,
        input_dir: str | PathLike[str],
        *,
        enabled: bool,
        logger: logging.Logger | None = None,
    ) -> None:
        self.input_dir = os.path.abspath(os.fspath(input_dir))
        self.enabled = bool(enabled)
        self._logger = logger if logger is not None else logging.getLogger(__name__)

    def retire_after_success(
        self,
        processing_succeeded: bool,
        predecessor_path: str | PathLike[str],
        successor_path: str | PathLike[str],
        protocol_name: str | None,
    ) -> bool:
        """Archive the actual predecessor only after its successor succeeded."""
        if not processing_succeeded or not self.enabled:
            return False

        predecessor = os.path.abspath(os.fspath(predecessor_path))
        successor = os.path.abspath(os.fspath(successor_path))
        predecessor_name = os.path.basename(predecessor)
        successor_name = os.path.basename(successor)

        # The first observed transform is processed against itself. Never move
        # that current transform, and independently protect the queue identity.
        if predecessor == successor:
            return False
        if _IDENTITY_TRANSFORM_RE.fullmatch(predecessor_name):
            return False
        if not protocol_name:
            self._logger.warning(
                "TRANSFORM RETIRE: predecessor left active (protocol unavailable): %s",
                predecessor_name,
            )
            return False
        if not _REGISTRATION_TRANSFORM_RE.fullmatch(predecessor_name):
            self._logger.warning(
                "TRANSFORM RETIRE: predecessor left active (unexpected name): %s",
                predecessor_name,
            )
            return False
        if os.path.dirname(predecessor) != self.input_dir:
            self._logger.warning(
                "TRANSFORM RETIRE: predecessor left active (not in acquisition root): %s",
                predecessor,
            )
            return False
        if not os.path.isfile(predecessor):
            self._logger.warning(
                "TRANSFORM RETIRE: predecessor unavailable: %s",
                predecessor,
            )
            return False

        processed_dir = processed_acquisition_directory(
            self.input_dir,
            protocol_name,
        )
        transform_dir = os.path.join(processed_dir, "transforms")
        destination = os.path.join(transform_dir, predecessor_name)
        try:
            os.makedirs(transform_dir, exist_ok=True)
        except OSError as error:
            self._logger.warning(
                "TRANSFORM RETIRE: archive directory unavailable for %s -> %s: %s",
                predecessor,
                destination,
                error,
            )
            return False

        if os.path.lexists(destination):
            self._logger.warning(
                "TRANSFORM RETIRE: collision for %s -> %s; source left active",
                predecessor,
                destination,
            )
            return False

        try:
            os.rename(predecessor, destination)
        except OSError as error:
            self._logger.warning(
                "TRANSFORM RETIRE: failed to archive %s -> %s: %s",
                predecessor,
                destination,
                error,
            )
            return False

        self._logger.info(
            "TRANSFORM RETIRE: archived %s after successor %s",
            predecessor_name,
            successor_name,
        )
        return True

"""Archive transforms after motion-monitor consumes their successors."""

from __future__ import annotations

import json
import logging
import os
import re
from dataclasses import dataclass
from os import PathLike

try:
    from .fire_moco_state import (
        load_fire_moco_state,
        parse_registration_transform_index,
    )
    from .retirement_paths import processed_acquisition_directory
except ImportError:  # Support execution from the motion_monitor application directory.
    from fire_moco_state import (
        load_fire_moco_state,
        parse_registration_transform_index,
    )
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


def retirement_mode_for_moco_flag(moco_flag: str | None = None) -> str:
    """Return the conservative transform-retirement mode for ``MOCO_FLAG``."""
    value = os.environ.get("MOCO_FLAG", "off") if moco_flag is None else moco_flag
    normalized = value.strip().lower()
    return normalized if normalized in {"off", "on"} else "disabled"


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


@dataclass(frozen=True)
class _PendingTransform:
    registration_index: int
    predecessor_path: str
    successor_path: str
    protocol_name: str | None


class TransformRetirementCoordinator:
    """Coordinate monitor last-use with FIRE release-through state.

    Pending candidates are intentionally in-memory only. This experimental
    optimization assumes uninterrupted motion-monitor operation and does not
    reconstruct pending retirement after restart.
    """

    def __init__(
        self,
        input_dir: str | PathLike[str],
        *,
        moco_flag: str | None,
        logger: logging.Logger | None = None,
    ) -> None:
        self.input_dir = os.path.abspath(os.fspath(input_dir))
        self.mode = retirement_mode_for_moco_flag(moco_flag)
        self._logger = logger if logger is not None else logging.getLogger(__name__)
        self._retirer = TransformRetirer(
            self.input_dir,
            enabled=self.mode in {"off", "on"},
            logger=self._logger,
        )
        self._pending: dict[int, _PendingTransform] = {}
        self._highest_fire_release_before: int | None = None
        self._last_state_error: str | None = None

    @property
    def pending_registration_indices(self) -> tuple[int, ...]:
        return tuple(sorted(self._pending))

    def reset(self) -> None:
        """Drop acquisition-local coordination state at the existing reset."""
        self._pending.clear()
        self._highest_fire_release_before = None
        self._last_state_error = None

    def retire_pending_released(self) -> int:
        """Recheck monitor-consumed candidates at the shutdown boundary."""
        if self.mode != "on" or not self._pending:
            return 0
        release_before = self._read_release_before()
        if release_before is None:
            return 0
        pending_before = len(self._pending)
        self._retire_released_pending(release_before)
        return pending_before - len(self._pending)

    def retire_after_success(
        self,
        processing_succeeded: bool,
        predecessor_path: str | PathLike[str],
        successor_path: str | PathLike[str],
        protocol_name: str | None,
    ) -> bool:
        """Retire a monitor-consumed predecessor under the configured mode."""
        if self.mode == "off":
            return self._retirer.retire_after_success(
                processing_succeeded,
                predecessor_path,
                successor_path,
                protocol_name,
            )
        if self.mode != "on" or not processing_succeeded:
            return False

        predecessor = os.path.abspath(os.fspath(predecessor_path))
        successor = os.path.abspath(os.fspath(successor_path))
        predecessor_name = os.path.basename(predecessor)

        # The first observed transform is processed against itself. Identity is
        # protected here as well as in the physical mover.
        if predecessor == successor or _IDENTITY_TRANSFORM_RE.fullmatch(predecessor_name):
            return False
        registration_index = parse_registration_transform_index(predecessor_name)
        if registration_index is None:
            self._logger.warning(
                "TRANSFORM RETIRE DECISION: transform=%s decision=pending_rejected "
                "reason=unexpected_name",
                predecessor_name,
            )
            return False

        candidate = _PendingTransform(
            registration_index=registration_index,
            predecessor_path=predecessor,
            successor_path=successor,
            protocol_name=protocol_name,
        )
        self._pending[registration_index] = candidate

        release_before = self._read_release_before()
        if release_before is None:
            self._log_decision(candidate, None, "retained_no_state")
            return False

        self._retire_released_pending(release_before)
        if registration_index in self._pending:
            # Eligible move failures are already recorded by the helper below.
            if registration_index >= release_before:
                self._log_decision(candidate, release_before, "pending")
            return False
        return True

    def _read_release_before(self) -> int | None:
        try:
            state = load_fire_moco_state(self.input_dir)
        except (OSError, json.JSONDecodeError, ValueError, TypeError) as error:
            message = str(error)
            if message != self._last_state_error:
                self._logger.warning(
                    "FIRE MOCO state unavailable; transform retirement delayed: %s",
                    error,
                )
                self._last_state_error = message
            return None

        release_before = state.release_before_registration_index
        if (
            self._highest_fire_release_before is not None
            and release_before < self._highest_fire_release_before
        ):
            message = (
                f"regressed from {self._highest_fire_release_before} to {release_before}"
            )
            if message != self._last_state_error:
                self._logger.warning(
                    "FIRE MOCO state regression ignored; transform retirement delayed: %s",
                    message,
                )
                self._last_state_error = message
            return None

        self._highest_fire_release_before = release_before
        self._last_state_error = None
        return release_before

    def _retire_released_pending(self, release_before: int) -> None:
        for registration_index in sorted(self._pending):
            if registration_index >= release_before:
                continue
            candidate = self._pending[registration_index]
            retired = self._retirer.retire_after_success(
                True,
                candidate.predecessor_path,
                candidate.successor_path,
                candidate.protocol_name,
            )
            if retired:
                del self._pending[registration_index]
                self._log_decision(candidate, release_before, "retired")
            else:
                self._log_decision(candidate, release_before, "retirement_failed")

    def _log_decision(
        self,
        candidate: _PendingTransform,
        release_before: int | None,
        decision: str,
    ) -> None:
        self._logger.info(
            "TRANSFORM RETIRE DECISION: transform=%s index=%d "
            "fire_release_before=%s decision=%s",
            os.path.basename(candidate.predecessor_path),
            candidate.registration_index,
            "unavailable" if release_before is None else release_before,
            decision,
        )

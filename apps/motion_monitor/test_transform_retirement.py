from __future__ import annotations

import json
import logging
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from apps.motion_monitor.retirement_paths import processed_acquisition_directory
from apps.motion_monitor.fire_moco_state import FIRE_MOCO_STATE_FILENAME
from apps.motion_monitor.transform_retirement import (
    TransformRetirementCoordinator,
    TransformRetirer,
    retirement_enabled_for_moco_flag,
    retirement_mode_for_moco_flag,
)


class TransformRetirementTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.input_dir = Path(self.temporary_directory.name)
        self.protocol = "func-bold_task-rest_run-01"
        self.logger = logging.getLogger(f"transform-retirement-test-{id(self)}")

    def _transform(
        self,
        index: int,
        volume: int,
        group: int,
        *,
        identity: bool = False,
        contents: bytes | None = None,
    ) -> Path:
        suffix = "_identity" if identity else ""
        path = self.input_dir / (
            f"alignTransform_{index:04d}_{volume:04d}-{group:04d}{suffix}.tfm"
        )
        path.write_bytes(contents if contents is not None else f"transform-{index}".encode())
        return path

    def _retirer(self, enabled: bool = True) -> TransformRetirer:
        return TransformRetirer(
            self.input_dir,
            enabled=enabled,
            logger=self.logger,
        )

    def _archive(self) -> Path:
        return Path(
            processed_acquisition_directory(self.input_dir, self.protocol)
        ) / "transforms"

    def _coordinator(self, moco_flag: str = "on") -> TransformRetirementCoordinator:
        return TransformRetirementCoordinator(
            self.input_dir,
            moco_flag=moco_flag,
            logger=self.logger,
        )

    def _write_fire_state(self, index: int, *, filename_index: int | None = None) -> None:
        filename_index = index if filename_index is None else filename_index
        payload = {
            "schema_version": 1,
            "committed_registration_index": index,
            "committed_transform_filename": (
                f"alignTransform_{filename_index:04d}_0001-0000.tfm"
            ),
            "release_before_registration_index": index,
        }
        (self.input_dir / FIRE_MOCO_STATE_FILENAME).write_text(json.dumps(payload))

    def test_successful_successor_retires_only_actual_predecessor(self) -> None:
        identity = self._transform(0, 0, 44, identity=True)
        predecessor = self._transform(1, 1, 0)
        successor = self._transform(2, 1, 7)
        self.assertFalse(self._archive().exists())

        retired = self._retirer().retire_after_success(
            True, predecessor, successor, self.protocol
        )

        self.assertTrue(retired)
        self.assertTrue((self._archive() / predecessor.name).is_file())
        self.assertFalse(predecessor.exists())
        self.assertTrue(successor.is_file())
        self.assertTrue(identity.is_file())

    def test_successive_retirement_leaves_identity_and_current_active(self) -> None:
        identity = self._transform(0, 0, 44, identity=True)
        first = self._transform(10, 5, 3)
        second = self._transform(11, 5, 9)
        third = self._transform(12, 6, 1)
        retirer = self._retirer()

        self.assertTrue(
            retirer.retire_after_success(True, first, second, self.protocol)
        )
        self.assertTrue(
            retirer.retire_after_success(True, second, third, self.protocol)
        )

        self.assertTrue((self._archive() / first.name).is_file())
        self.assertTrue((self._archive() / second.name).is_file())
        self.assertTrue(third.is_file())
        self.assertTrue(identity.is_file())

    def test_moco_enabled_disables_retirement(self) -> None:
        predecessor = self._transform(1, 1, 0)
        successor = self._transform(2, 1, 1)
        self.assertTrue(retirement_enabled_for_moco_flag("off"))
        self.assertTrue(retirement_enabled_for_moco_flag(" OFF "))
        self.assertFalse(retirement_enabled_for_moco_flag("on"))
        self.assertFalse(retirement_enabled_for_moco_flag("unexpected"))

        retired = self._retirer(enabled=False).retire_after_success(
            True, predecessor, successor, self.protocol
        )

        self.assertFalse(retired)
        self.assertTrue(predecessor.is_file())
        self.assertTrue(successor.is_file())
        self.assertFalse(self._archive().exists())

    def test_identity_and_current_transform_are_never_retired(self) -> None:
        identity = self._transform(0, 0, 44, identity=True)
        current = self._transform(1, 1, 0)
        retirer = self._retirer()

        self.assertFalse(
            retirer.retire_after_success(True, identity, current, self.protocol)
        )
        self.assertFalse(
            retirer.retire_after_success(True, current, current, self.protocol)
        )

        self.assertTrue(identity.is_file())
        self.assertTrue(current.is_file())
        self.assertFalse(self._archive().exists())

    def test_processing_failure_leaves_predecessor_active(self) -> None:
        predecessor = self._transform(1, 1, 0)
        successor = self._transform(2, 1, 1)

        retired = self._retirer().retire_after_success(
            False, predecessor, successor, self.protocol
        )

        self.assertFalse(retired)
        self.assertTrue(predecessor.is_file())
        self.assertFalse(self._archive().exists())

    def test_archive_collision_does_not_overwrite_or_move_source(self) -> None:
        predecessor = self._transform(1, 1, 0, contents=b"current acquisition")
        successor = self._transform(2, 1, 1)
        self._archive().mkdir(parents=True)
        destination = self._archive() / predecessor.name
        destination.write_bytes(b"earlier acquisition")

        retired = self._retirer().retire_after_success(
            True, predecessor, successor, self.protocol
        )

        self.assertFalse(retired)
        self.assertEqual(destination.read_bytes(), b"earlier acquisition")
        self.assertEqual(predecessor.read_bytes(), b"current acquisition")

    def test_move_failure_is_nonfatal_and_leaves_source_active(self) -> None:
        predecessor = self._transform(1, 1, 0)
        successor = self._transform(2, 1, 1)

        with patch(
            "apps.motion_monitor.transform_retirement.os.rename",
            side_effect=OSError("injected move failure"),
        ):
            retired = self._retirer().retire_after_success(
                True, predecessor, successor, self.protocol
            )

        self.assertFalse(retired)
        self.assertTrue(predecessor.is_file())
        self.assertTrue(successor.is_file())

    def test_directory_creation_failure_is_nonfatal(self) -> None:
        predecessor = self._transform(1, 1, 0)
        successor = self._transform(2, 1, 1)

        with patch(
            "apps.motion_monitor.transform_retirement.os.makedirs",
            side_effect=OSError("injected directory failure"),
        ):
            retired = self._retirer().retire_after_success(
                True, predecessor, successor, self.protocol
            )

        self.assertFalse(retired)
        self.assertTrue(predecessor.is_file())
        self.assertFalse(self._archive().exists())

    def test_noncontiguous_lifo_sequence_uses_passed_predecessor(self) -> None:
        predecessor = self._transform(4, 2, 3)
        unrelated = self._transform(87, 40, 20)
        successor = self._transform(103, 41, 7)

        retired = self._retirer().retire_after_success(
            True, predecessor, successor, self.protocol
        )

        self.assertTrue(retired)
        self.assertTrue((self._archive() / predecessor.name).is_file())
        self.assertTrue(unrelated.is_file())
        self.assertTrue(successor.is_file())

    def test_transform_archive_coexists_with_image_retirement_layout(self) -> None:
        processed_dir = Path(
            processed_acquisition_directory(self.input_dir, self.protocol)
        )
        processed_dir.mkdir()
        image_names = {
            f"{self.protocol}_volume_0001_group_0000.txt",
            f"{self.protocol}_volume_0001_slice_0000.nhdr",
            f"{self.protocol}_volume_0001_slice_0000.raw",
        }
        for name in image_names:
            (processed_dir / name).write_bytes(b"image work")
        predecessor = self._transform(1, 1, 0)
        successor = self._transform(2, 1, 7)

        self.assertTrue(
            self._retirer().retire_after_success(
                True, predecessor, successor, self.protocol
            )
        )

        self.assertTrue(image_names.issubset({path.name for path in processed_dir.iterdir()}))
        self.assertEqual(
            {path.name for path in self._archive().iterdir()},
            {predecessor.name},
        )

    def test_moco_on_missing_or_malformed_state_keeps_candidate_pending(self) -> None:
        predecessor = self._transform(100, 4, 0)
        successor = self._transform(101, 4, 1)
        coordinator = self._coordinator()

        self.assertFalse(
            coordinator.retire_after_success(
                True, predecessor, successor, self.protocol
            )
        )
        self.assertEqual(coordinator.pending_registration_indices, (100,))
        self.assertTrue(predecessor.exists())

        (self.input_dir / FIRE_MOCO_STATE_FILENAME).write_text("{malformed")
        self.assertFalse(
            coordinator.retire_after_success(
                True, predecessor, successor, self.protocol
            )
        )
        self.assertTrue(predecessor.exists())

    def test_fire_anchor_and_newer_candidates_remain_pending(self) -> None:
        coordinator = self._coordinator()
        self._write_fire_state(105)
        anchor = self._transform(105, 5, 0)
        newer = self._transform(106, 5, 1)
        successor = self._transform(109, 5, 2)

        self.assertFalse(
            coordinator.retire_after_success(True, anchor, successor, self.protocol)
        )
        self.assertFalse(
            coordinator.retire_after_success(True, newer, successor, self.protocol)
        )
        self.assertEqual(coordinator.pending_registration_indices, (105, 106))
        self.assertTrue(anchor.exists())
        self.assertTrue(newer.exists())

    def test_lower_candidate_retires_but_current_transform_stays_active(self) -> None:
        coordinator = self._coordinator()
        self._write_fire_state(105)
        predecessor = self._transform(104, 5, 0)
        current = self._transform(105, 5, 1)

        self.assertTrue(
            coordinator.retire_after_success(
                True, predecessor, current, self.protocol
            )
        )
        self.assertTrue((self._archive() / predecessor.name).exists())
        self.assertTrue(current.exists())

    def test_pending_candidates_retire_after_later_fire_jump_with_gaps(self) -> None:
        coordinator = self._coordinator()
        self._write_fire_state(100)
        first = self._transform(101, 10, 0)
        second = self._transform(103, 10, 2)
        third = self._transform(105, 10, 4)
        current = self._transform(106, 10, 5)

        for predecessor in (first, second, third):
            self.assertFalse(
                coordinator.retire_after_success(
                    True, predecessor, current, self.protocol
                )
            )
        self.assertEqual(coordinator.pending_registration_indices, (101, 103, 105))

        self._write_fire_state(106)
        # A later successful monitor event is the deliberate recheck trigger.
        self.assertFalse(
            coordinator.retire_after_success(True, current, current, self.protocol)
        )
        # The first-observation self-pair is not a consumption event, so force a
        # real successor and verify all previously consumed gaps are released.
        later = self._transform(109, 11, 0)
        self.assertFalse(
            coordinator.retire_after_success(True, current, later, self.protocol)
        )

        self.assertEqual(coordinator.pending_registration_indices, (106,))
        for predecessor in (first, second, third):
            self.assertTrue((self._archive() / predecessor.name).exists())
        self.assertTrue(current.exists())
        self.assertTrue(later.exists())

    def test_shutdown_recheck_retires_only_fire_released_pending(self) -> None:
        coordinator = self._coordinator()
        self._write_fire_state(100)
        released = self._transform(101, 10, 0)
        anchor = self._transform(105, 10, 4)
        current = self._transform(106, 10, 5)
        coordinator.retire_after_success(True, released, current, self.protocol)
        coordinator.retire_after_success(True, anchor, current, self.protocol)
        self.assertEqual(coordinator.pending_registration_indices, (101, 105))

        self._write_fire_state(105)
        self.assertEqual(coordinator.retire_pending_released(), 1)

        self.assertTrue((self._archive() / released.name).is_file())
        self.assertTrue(anchor.is_file())
        self.assertEqual(coordinator.pending_registration_indices, (105,))

    def test_identity_processing_failure_and_unexpected_flag_are_conservative(self) -> None:
        identity = self._transform(0, 0, 44, identity=True)
        predecessor = self._transform(10, 2, 0)
        successor = self._transform(11, 2, 1)
        self._write_fire_state(100)
        coordinator = self._coordinator()

        self.assertFalse(
            coordinator.retire_after_success(True, identity, successor, self.protocol)
        )
        self.assertFalse(
            coordinator.retire_after_success(False, predecessor, successor, self.protocol)
        )
        self.assertTrue(identity.exists())
        self.assertTrue(predecessor.exists())
        self.assertEqual(coordinator.pending_registration_indices, ())

        disabled = self._coordinator("unexpected")
        self.assertEqual(retirement_mode_for_moco_flag("unexpected"), "disabled")
        self.assertFalse(
            disabled.retire_after_success(True, predecessor, successor, self.protocol)
        )
        self.assertTrue(predecessor.exists())

    def test_moco_off_path_does_not_require_fire_state(self) -> None:
        predecessor = self._transform(10, 2, 0)
        successor = self._transform(11, 2, 1)
        coordinator = self._coordinator(" OFF ")

        self.assertTrue(
            coordinator.retire_after_success(
                True, predecessor, successor, self.protocol
            )
        )
        self.assertTrue((self._archive() / predecessor.name).exists())

    def test_move_failure_and_collision_remain_pending_without_overwrite(self) -> None:
        coordinator = self._coordinator()
        self._write_fire_state(100)
        failed = self._transform(10, 2, 0)
        current = self._transform(11, 2, 1)
        with patch(
            "apps.motion_monitor.transform_retirement.os.rename",
            side_effect=OSError("injected move failure"),
        ):
            self.assertFalse(
                coordinator.retire_after_success(
                    True, failed, current, self.protocol
                )
            )
        self.assertEqual(coordinator.pending_registration_indices, (10,))
        self.assertTrue(failed.exists())

        collision = self._transform(12, 2, 2, contents=b"active")
        self._archive().mkdir(parents=True, exist_ok=True)
        destination = self._archive() / collision.name
        destination.write_bytes(b"archived")
        self.assertFalse(
            coordinator.retire_after_success(
                True, collision, current, self.protocol
            )
        )
        self.assertEqual(destination.read_bytes(), b"archived")
        self.assertEqual(collision.read_bytes(), b"active")
        # The next successful monitor event retries the earlier nonfatal move.
        self.assertTrue((self._archive() / failed.name).exists())
        self.assertEqual(coordinator.pending_registration_indices, (12,))

    def test_highest_generated_transform_and_mtime_do_not_release_candidates(self) -> None:
        coordinator = self._coordinator()
        self._write_fire_state(50)
        candidate = self._transform(60, 3, 0)
        current = self._transform(999, 99, 9)
        os.utime(current, (999999, 999999))
        (self.input_dir / "registration_status.json").write_text(
            json.dumps({"highest_successful_registration_index": 999})
        )

        self.assertFalse(
            coordinator.retire_after_success(
                True, candidate, current, self.protocol
            )
        )
        self.assertTrue(candidate.exists())
        self.assertEqual(coordinator.pending_registration_indices, (60,))

    def test_filename_mismatch_and_regressing_state_delay_retirement(self) -> None:
        coordinator = self._coordinator()
        candidate = self._transform(10, 3, 0)
        current = self._transform(11, 3, 1)
        self._write_fire_state(100, filename_index=99)
        self.assertFalse(
            coordinator.retire_after_success(
                True, candidate, current, self.protocol
            )
        )
        self.assertTrue(candidate.exists())

        self._write_fire_state(105)
        released = self._transform(104, 4, 0)
        later = self._transform(105, 4, 1)
        self.assertTrue(
            coordinator.retire_after_success(True, released, later, self.protocol)
        )

        self._write_fire_state(104)
        regressed_candidate = self._transform(103, 4, 2)
        self.assertFalse(
            coordinator.retire_after_success(
                True, regressed_candidate, later, self.protocol
            )
        )
        self.assertTrue(regressed_candidate.exists())
        self.assertIn(103, coordinator.pending_registration_indices)


if __name__ == "__main__":
    unittest.main()

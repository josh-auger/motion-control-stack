from __future__ import annotations

import logging
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from apps.motion_monitor.retirement_paths import processed_acquisition_directory
from apps.motion_monitor.transform_retirement import (
    TransformRetirer,
    retirement_enabled_for_moco_flag,
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


if __name__ == "__main__":
    unittest.main()

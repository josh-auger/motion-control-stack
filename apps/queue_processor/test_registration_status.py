"""Focused tests for the queue-to-monitor registration status contract."""

import json
import logging
import os
import tempfile
import unittest

from apps.queue_processor.registration_status import RegistrationStatusTracker


class RegistrationStatusTrackerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.status_path = os.path.join(
            self.temporary_directory.name, "registration_status.json"
        )
        logger = logging.Logger(f"registration-status-test-{id(self)}")
        logger.addHandler(logging.NullHandler())
        self.tracker = RegistrationStatusTracker(self.status_path, logger=logger)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _read_status(self) -> dict:
        with open(self.status_path, "r", encoding="utf-8") as status_file:
            return json.load(status_file)

    def test_registered_skipped_and_failed_groups_are_published(self) -> None:
        self.tracker.set_expected_groups(28, {0, 1, 2})
        self.tracker.record_group(
            28,
            0,
            "registered",
            transform_filename="alignTransform_0001_0028-0000.tfm",
        )
        self.tracker.record_group(28, 1, "skipped")
        self.assertFalse(self.tracker.finalize_volume_if_ready(28))
        self.assertFalse(os.path.exists(self.status_path))

        self.tracker.record_group(28, 2, "failed")
        self.assertTrue(self.tracker.finalize_volume_if_ready(28))
        groups = self._read_status()["volumes"]["28"]["groups"]
        self.assertEqual(groups["0"]["status"], "registered")
        self.assertEqual(
            groups["0"]["transform_filename"],
            "alignTransform_0001_0028-0000.tfm",
        )
        self.assertEqual(groups["1"], {"status": "skipped"})
        self.assertEqual(groups["2"], {"status": "failed"})

    def test_no_group_can_be_recorded_after_volume_finalization(self) -> None:
        self.tracker.set_expected_groups(1, {0})
        self.tracker.record_group(1, 0, "skipped")
        self.assertTrue(self.tracker.finalize_volume_if_ready(1))

        with self.assertRaisesRegex(RuntimeError, "finalized"):
            self.tracker.record_group(1, 0, "failed")
        with self.assertRaisesRegex(RuntimeError, "finalized"):
            self.tracker.record_group(1, 1, "skipped")

        restored_tracker = RegistrationStatusTracker(self.status_path)
        self.assertTrue(restored_tracker.is_finalized(1))
        with self.assertRaisesRegex(RuntimeError, "finalized"):
            restored_tracker.record_group(1, 0, "failed")

    def test_finalized_history_is_retained_after_later_updates(self) -> None:
        for volume_number in (1, 2):
            self.tracker.set_expected_groups(volume_number, {0})
            self.tracker.record_group(volume_number, 0, "skipped")
            self.assertTrue(self.tracker.finalize_volume_if_ready(volume_number))

        document = self._read_status()
        self.assertEqual(document["schema_version"], 1)
        self.assertEqual(set(document["volumes"]), {"1", "2"})

    def test_close_finalization_publishes_the_final_acquisition_volume(self) -> None:
        self.tracker.set_expected_groups(192, {0, 1})
        self.tracker.record_group(192, 0, "skipped")
        self.tracker.record_group(192, 1, "failed")

        self.assertEqual(self.tracker.finalize_all_ready(), [192])
        record = self._read_status()["volumes"]["192"]
        self.assertTrue(record["registration_finalized"])
        self.assertEqual(set(record["groups"]), {"0", "1"})

    def test_atomic_publication_always_leaves_valid_canonical_json(self) -> None:
        self.tracker.set_expected_groups(3, {0})
        self.tracker.record_group(3, 0, "skipped")
        self.tracker.finalize_volume_if_ready(3)
        first_document = self._read_status()

        self.tracker.set_expected_groups(4, {0})
        self.tracker.record_group(4, 0, "failed")
        self.tracker.finalize_volume_if_ready(4)
        second_document = self._read_status()

        self.assertEqual(set(first_document["volumes"]), {"3"})
        self.assertEqual(set(second_document["volumes"]), {"3", "4"})
        temporary_files = [
            name
            for name in os.listdir(self.temporary_directory.name)
            if name != "registration_status.json"
        ]
        self.assertEqual(temporary_files, [])


if __name__ == "__main__":
    unittest.main()

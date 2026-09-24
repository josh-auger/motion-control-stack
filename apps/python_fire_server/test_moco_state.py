"""Focused tests for atomic FIRE MOCO release-state publication."""

import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import moco_state
from moco_state import FIRE_MOCO_STATE_FILENAME, FireMocoStatePublisher


def transform_name(index, volume=1, group=0):
    return f"alignTransform_{index:04d}_{volume:04d}-{group:04d}.tfm"


def publish(publisher, index):
    return publisher.publish(
        committed_registration_index=index,
        committed_transform_filename=transform_name(index, volume=2, group=3),
        feedback_frame_number=4,
        trigger_image_identifier=88,
        trigger_volume=2,
        trigger_slice=3,
        trigger_group=3,
    )


class FireMocoStatePublisherTests(unittest.TestCase):
    def test_state_absent_until_first_publish_and_schema_matches_commit(self):
        with tempfile.TemporaryDirectory() as directory:
            publisher = FireMocoStatePublisher(directory)
            state_path = Path(directory, FIRE_MOCO_STATE_FILENAME)
            self.assertFalse(state_path.exists())

            self.assertTrue(publish(publisher, 105))

            payload = json.loads(state_path.read_text())
            self.assertEqual(payload["schema_version"], 1)
            self.assertEqual(payload["committed_registration_index"], 105)
            self.assertEqual(payload["release_before_registration_index"], 105)
            self.assertEqual(
                payload["committed_transform_filename"],
                transform_name(105, volume=2, group=3),
            )
            self.assertEqual(payload["feedback_frame_number"], 4)

    def test_publication_uses_atomic_replace_and_leaves_no_temporary_file(self):
        with tempfile.TemporaryDirectory() as directory:
            publisher = FireMocoStatePublisher(directory)
            real_replace = os.replace
            with patch.object(moco_state.os, "replace", wraps=real_replace) as replace:
                self.assertTrue(publish(publisher, 7))

            replace.assert_called_once()
            temporary, final = replace.call_args.args
            self.assertTrue(temporary.endswith(".tmp"))
            self.assertEqual(final, str(Path(directory, FIRE_MOCO_STATE_FILENAME)))
            self.assertFalse(Path(temporary).exists())
            self.assertEqual(
                sorted(path.name for path in Path(directory).iterdir()),
                [FIRE_MOCO_STATE_FILENAME],
            )
            self.assertNotEqual(Path(final).suffix, ".json")

    def test_regression_is_refused_and_existing_state_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            publisher = FireMocoStatePublisher(directory)
            self.assertTrue(publish(publisher, 105))
            state_path = Path(directory, FIRE_MOCO_STATE_FILENAME)
            original = state_path.read_bytes()

            # A fresh publisher also validates a pre-existing final state.
            restarted_publisher = FireMocoStatePublisher(directory)
            self.assertFalse(publish(restarted_publisher, 104))
            self.assertEqual(state_path.read_bytes(), original)

    def test_publication_failure_is_nonfatal_and_later_state_catches_up(self):
        with tempfile.TemporaryDirectory() as directory:
            publisher = FireMocoStatePublisher(directory)
            with patch.object(
                moco_state.os,
                "replace",
                side_effect=OSError("simulated replace failure"),
            ):
                self.assertFalse(publish(publisher, 105))
            self.assertFalse(Path(directory, FIRE_MOCO_STATE_FILENAME).exists())
            self.assertFalse(any(path.suffix == ".tmp" for path in Path(directory).iterdir()))

            self.assertTrue(publish(publisher, 109))
            payload = json.loads(Path(directory, FIRE_MOCO_STATE_FILENAME).read_text())
            self.assertEqual(payload["committed_registration_index"], 109)

    def test_invalid_or_identity_filename_cannot_be_published(self):
        with tempfile.TemporaryDirectory() as directory:
            publisher = FireMocoStatePublisher(directory)
            for filename in (
                "alignTransform_0000_0000-0000_identity.tfm",
                transform_name(6),
            ):
                with self.subTest(filename=filename):
                    self.assertFalse(
                        publisher.publish(
                            committed_registration_index=5,
                            committed_transform_filename=filename,
                            feedback_frame_number=1,
                            trigger_image_identifier=1,
                            trigger_volume=1,
                            trigger_slice=1,
                            trigger_group=1,
                        )
                    )
            self.assertFalse(Path(directory, FIRE_MOCO_STATE_FILENAME).exists())


if __name__ == "__main__":
    unittest.main()

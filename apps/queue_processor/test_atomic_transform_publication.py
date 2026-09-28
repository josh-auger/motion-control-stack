"""Focused tests for atomic identity-transform publication."""

from __future__ import annotations

import multiprocessing
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import SimpleITK as sitk

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))

from apps.queue_processor import process_queue_directory as queue


def _observe_final(final_path, ready, done, observations):
    """Observe a publication from a separate process on the local test FS."""
    ready.set()
    seen = []
    while not done.is_set():
        try:
            seen.append(Path(final_path).read_text())
        except FileNotFoundError:
            pass
        time.sleep(0.001)
    try:
        seen.append(Path(final_path).read_text())
    except FileNotFoundError:
        pass
    observations.put(seen)


class AtomicTransformPublicationTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.final = self.root / "alignTransform_0001_0001-0000.tfm"
        self.transform = sitk.VersorRigid3DTransform()
        self.transform.SetTranslation((1.25, -2.5, 3.75))

    def test_success_is_parseable_and_staging_is_same_directory(self):
        real_write = sitk.WriteTransform
        observed = {}

        def recording_write(transform, path):
            observed["stage"] = Path(path)
            self.assertEqual(Path(path).parent, self.final.parent)
            self.assertEqual(Path(path).suffix, ".tfm")
            self.assertTrue(Path(path).name.startswith("."))
            self.assertTrue(Path(path).name.endswith(".partial.tfm"))
            self.assertFalse(self.final.exists())
            real_write(transform, path)

        with patch.object(queue.sitk, "WriteTransform", side_effect=recording_write):
            queue.write_transform_atomically(self.transform, self.final)

        self.assertTrue(self.final.is_file())
        self.assertFalse(observed["stage"].exists())
        loaded = sitk.ReadTransform(str(self.final))
        self.assertEqual(tuple(loaded.GetParameters()), tuple(self.transform.GetParameters()))

    def test_existing_final_is_preserved_and_writer_is_not_called(self):
        self.final.write_bytes(b"existing-transform")
        with patch.object(queue.sitk, "WriteTransform") as writer:
            with self.assertRaises(FileExistsError):
                queue.write_transform_atomically(self.transform, self.final)
        writer.assert_not_called()
        self.assertEqual(self.final.read_bytes(), b"existing-transform")

    def test_serialization_failure_cleans_staging_and_propagates(self):
        observed = {}

        def failing_write(_transform, path):
            observed["stage"] = Path(path)
            Path(path).write_bytes(b"partial")
            raise RuntimeError("injected serialization failure")

        with patch.object(queue.sitk, "WriteTransform", side_effect=failing_write):
            with self.assertRaisesRegex(RuntimeError, "serialization failure"):
                queue.write_transform_atomically(self.transform, self.final)

        self.assertFalse(self.final.exists())
        self.assertFalse(observed["stage"].exists())

    def test_rename_failure_cleans_staging_and_propagates(self):
        observed = {}
        real_write = sitk.WriteTransform

        def recording_write(transform, path):
            observed["stage"] = Path(path)
            real_write(transform, path)

        with patch.object(queue.sitk, "WriteTransform", side_effect=recording_write), \
             patch.object(queue.os, "replace", side_effect=OSError("injected rename failure")):
            with self.assertRaisesRegex(OSError, "rename failure"):
                queue.write_transform_atomically(self.transform, self.final)

        self.assertFalse(self.final.exists())
        self.assertFalse(observed["stage"].exists())

    def test_cleanup_failure_does_not_hide_serialization_failure(self):
        def failing_write(_transform, path):
            Path(path).write_bytes(b"partial")
            raise RuntimeError("primary serialization failure")

        with patch.object(queue.sitk, "WriteTransform", side_effect=failing_write), \
             patch.object(queue.os, "remove", side_effect=PermissionError("cleanup failure")), \
             self.assertLogs(level="WARNING") as captured:
            with self.assertRaisesRegex(RuntimeError, "primary serialization failure"):
                queue.write_transform_atomically(self.transform, self.final)

        self.assertTrue(any("cleanup failure" in line for line in captured.output))
        self.assertFalse(self.final.exists())

    @unittest.skipUnless(hasattr(os, "fork"), "requires fork-based local observer")
    def test_separate_process_observer_never_sees_partial_final(self):
        context = multiprocessing.get_context("fork")
        ready = context.Event()
        done = context.Event()
        observations = context.Queue()
        observer = context.Process(
            target=_observe_final,
            args=(str(self.final), ready, done, observations),
        )
        observer.start()
        self.assertTrue(ready.wait(2.0))

        def slow_write(_transform, path):
            with open(path, "w") as output:
                output.write("partial")
                output.flush()
                os.fsync(output.fileno())
                time.sleep(0.05)
                output.write("-complete")
                output.flush()
                os.fsync(output.fileno())

        try:
            with patch.object(queue.sitk, "WriteTransform", side_effect=slow_write):
                queue.write_transform_atomically(self.transform, self.final)
        finally:
            done.set()
            observer.join(2.0)
            if observer.is_alive():
                observer.terminate()
                observer.join()

        self.assertEqual(observer.exitcode, 0)
        seen = observations.get(timeout=2.0)
        self.assertTrue(seen)
        self.assertEqual(set(seen), {"partial-complete"})

    def test_identity_helper_publishes_parseable_transform_in_requested_root(self):
        reference = self.root / "reference.nrrd"
        sitk.WriteImage(sitk.Image([2, 2, 2], sitk.sitkFloat32), str(reference))

        result = queue.create_identityTransformFile(
            self.root, reference, 0, 0, 0
        )

        result_path = Path(result)
        self.assertEqual(result_path.parent, self.root)
        self.assertEqual(
            result_path.name,
            "alignTransform_0000_0000-0000_identity.tfm",
        )
        self.assertEqual(sitk.ReadTransform(str(result_path)).GetName(),
                         "VersorRigid3DTransform")

    def test_predecessor_glob_cannot_select_hidden_staging_name(self):
        identity = self.root / "alignTransform_0000_0000-0000_identity.tfm"
        identity.write_text("identity")
        stage = self.root / ".alignTransform_0001_0001-0000.token.partial.tfm"
        stage.write_text("partial")

        chosen = queue.select_input_transform(
            self.root, str(identity), 2, reference_volume_flag=1
        )

        self.assertEqual(chosen, str(identity))
        self.assertTrue(stage.exists())


class IdentityStatusBoundaryTests(unittest.TestCase):
    def test_identity_publication_failure_is_not_recorded_as_registered(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            (root / "protocol_volume_0000_group_0000.txt").write_text(
                "reference.nrrd\n"
            )

            with patch.object(queue, "resample_nrrd_volume", return_value="reference.nrrd"), \
                 patch.object(queue, "create_identityTransformFile",
                              side_effect=OSError("publication failed")), \
                 patch.object(queue, "read_slice_timings_from_json", return_value=[0.0]), \
                 patch.object(queue, "reset_logging", return_value=None), \
                 patch.object(
                     queue.RegistrationStatusTracker,
                     "record_group",
                     autospec=True,
                 ) as record_group:
                with self.assertRaisesRegex(OSError, "publication failed"):
                    queue.monitor_directory(root, "on", "cuda")

            record_group.assert_not_called()

            status_path = root / queue.REGISTRATION_STATUS_FILENAME
            if status_path.exists():
                status = json.loads(status_path.read_text())
                for volume in status.get("volumes", {}).values():
                    for group in volume.get("groups", {}).values():
                        self.assertNotEqual(group.get("status"), "registered")


if __name__ == "__main__":
    unittest.main()

"""Tests for transform publication names at the monitor entry point."""

from __future__ import annotations

import inspect
import os
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import MagicMock, patch


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import monitor_directory as monitor
from apps.python_fire_server.pointer_file import write_pointer_file_atomically
from monitor_directory import is_published_transform_filename


class StopAfterPointerProcessing(Exception):
    pass


class StopAfterTransformProcessing(Exception):
    pass


class StopAfterIdlePoll(Exception):
    pass


class TransformPublicationFilterTests(unittest.TestCase):
    def test_final_registration_and_identity_names_are_accepted(self):
        published_names = (
            ("cuda", "alignTransform_0001_0002-0003.tfm"),
            ("cpu", "alignTransform_0042_0002-0004.tfm"),
            ("identity", "alignTransform_0000_0000-0000_identity.tfm"),
        )
        for producer, filename in published_names:
            with self.subTest(producer=producer):
                self.assertTrue(is_published_transform_filename(filename))

    def test_staging_and_malformed_tfm_names_are_rejected(self):
        invalid = (
            ".alignTransform_0001_0002-0003.token.partial.tfm",
            "alignTransform_0001_0002-0003.partial.tfm",
            "alignTransform_0001.tfm",
            "unrelated.tfm",
        )
        for filename in invalid:
            with self.subTest(filename=filename):
                self.assertFalse(is_published_transform_filename(filename))

    def _assert_published_transform_bypasses_wait(self, filename, transform):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.json"
            transform_path = root / filename
            metadata.write_text("{}")
            monitor.sitk.WriteTransform(transform, transform_path)
            os.utime(metadata, (1, 1))
            os.utime(transform_path, (2, 2))
            metadata_object = {
                "SliceTiming": {"0": 0.0},
                "measurementInformation": {"protocolName": "test"},
                "encoding": [
                    {"encodingLimits": {"repetition": {"maximum": 1}}}
                ],
            }
            retirer = MagicMock()
            retirer.mode = "test"
            retirer.retire_after_success.side_effect = StopAfterTransformProcessing

            with patch.object(monitor, "reset_logging", return_value=None), \
                 patch.object(
                     monitor,
                     "load_metadata_from_json",
                     return_value=metadata_object,
                 ), \
                 patch.object(
                     monitor,
                     "TransformRetirementCoordinator",
                     return_value=retirer,
                 ), \
                 patch.object(
                     monitor,
                     "dashboard_profiling_enabled",
                     return_value=False,
                 ), \
                 patch.object(monitor, "load_dashboard_interval", return_value=3600), \
                 patch.object(
                     monitor.time,
                     "sleep",
                     side_effect=AssertionError("published transform entered wait"),
                 ) as sleep:
                with self.assertRaises(StopAfterTransformProcessing):
                    monitor.monitor_directory(root, 50, 1.0, 5000, "off")

            sleep.assert_not_called()
            retirer.retire_after_success.assert_called_once_with(
                True,
                os.fspath(transform_path),
                os.fspath(transform_path),
                "test",
            )

    def test_published_final_transforms_bypass_wait_and_reach_fd_state(self):
        cases = (
            (
                "cuda",
                "alignTransform_0001_0001-0000.tfm",
                monitor.sitk.Euler3DTransform(),
            ),
            (
                "cpu",
                "alignTransform_0042_0001-0001.tfm",
                monitor.sitk.VersorRigid3DTransform(),
            ),
            (
                "identity",
                "alignTransform_0000_0000-0000_identity.tfm",
                monitor.sitk.VersorRigid3DTransform(),
            ),
        )
        for producer, filename, transform in cases:
            with self.subTest(producer=producer):
                self._assert_published_transform_bypasses_wait(filename, transform)

    def test_atomically_published_pointer_bypasses_wait_and_reaches_tsnr(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "test_volume_0001_slice_0000.raw"
            header = root / "test_volume_0001_slice_0000.nhdr"
            pointer = root / "test_volume_0001_group_0000.txt"
            raw.write_bytes(b"complete image payload")
            header.write_text(
                "NRRD0004\n"
                "type: uchar\n"
                "dimension: 3\n"
                "sizes: 1 1 1\n"
                "encoding: raw\n"
                f"data file: {raw.name}\n"
            )
            write_pointer_file_atomically(pointer, [header.name])
            temporary_pointer = root / f"{pointer.name}.deadbeef.tmp"
            temporary_pointer.write_text("not-yet-published.nhdr\n")
            tsnr_processor = MagicMock()
            tsnr_processor.handle_pointer.side_effect = StopAfterPointerProcessing

            with patch.object(monitor, "reset_logging", return_value=None), \
                 patch.object(
                     monitor,
                     "TSNRVolumeProcessor",
                     return_value=tsnr_processor,
                 ), \
                 patch.object(
                     monitor.time,
                     "sleep",
                     side_effect=AssertionError("published pointer entered wait"),
                 ) as sleep:
                with self.assertRaises(StopAfterPointerProcessing):
                    monitor.monitor_directory(root, 50, 1.0, 5000, "off")

            sleep.assert_not_called()
            tsnr_processor.handle_pointer.assert_called_once_with(
                os.fspath(pointer),
                None,
            )
            self.assertEqual(pointer.read_text(), f"{header.name}\n")
            self.assertEqual(raw.read_bytes(), b"complete image payload")
            self.assertTrue(temporary_pointer.is_file())

    def test_rejected_tfm_files_never_enter_transform_processing(self):
        rejected_names = (
            ".alignTransform_0001_0002-0003.token.partial.tfm",
            "unrelated.tfm",
        )
        for filename in rejected_names:
            with self.subTest(filename=filename), \
                 tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                rejected = root / filename
                rejected.write_text("not published")
                retirer = MagicMock()
                retirer.mode = "test"

                def stop_after_idle_poll(_delay):
                    stack_functions = {
                        frame.function for frame in inspect.stack()
                    }
                    self.assertNotIn("wait_for_complete_write", stack_functions)
                    raise StopAfterIdlePoll

                with patch.object(
                         monitor,
                         "TransformRetirementCoordinator",
                         return_value=retirer,
                     ), \
                     patch.object(
                         monitor,
                         "read_transform_as_euler",
                     ) as read_transform, \
                     patch.object(monitor, "compose_transform_pair") as compose, \
                     patch.object(monitor, "calculate_displacement") as displacement, \
                     patch.object(
                         monitor,
                         "load_metadata_from_json",
                     ) as load_metadata, \
                     patch.object(
                         monitor.time,
                         "sleep",
                         side_effect=stop_after_idle_poll,
                     ) as sleep:
                    with self.assertRaises(StopAfterIdlePoll):
                        monitor.monitor_directory(root, 50, 1.0, 5000, "off")

                sleep.assert_called_once_with(0.005)
                load_metadata.assert_not_called()
                read_transform.assert_not_called()
                compose.assert_not_called()
                displacement.assert_not_called()
                retirer.retire_after_success.assert_not_called()


if __name__ == "__main__":
    unittest.main()

"""Focused tests for the active motion-monitor TSNR orchestration."""

import logging
import json
import os
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import SimpleITK as sitk

from apps.motion_monitor.running_tsnr import RunningTSNR
from apps.motion_monitor.tsnr_integration import TSNRVolumeProcessor


class FailingAccumulator(RunningTSNR):
    """Accumulator double that accepts the reference and rejects later data."""

    def update(self, volume: np.ndarray) -> None:
        if self.count:
            raise RuntimeError("injected accumulator failure")
        super().update(volume)


class TSNRIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.input_dir = self.temporary_directory.name
        self.output_path = os.path.join(self.input_dir, "tsnr_dashboard.jpg")
        self.protocol = "test_protocol"
        self.logger = logging.Logger(f"tsnr-test-{id(self)}")
        self.logger.addHandler(logging.NullHandler())
        self.logger.propagate = False

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _processor(self, **kwargs) -> TSNRVolumeProcessor:
        return TSNRVolumeProcessor(
            self.input_dir,
            self.output_path,
            logger=self.logger,
            **kwargs,
        )

    def _pointer(self, volume: int, group: int = 0, contents: str = "") -> str:
        path = os.path.join(
            self.input_dir,
            f"{self.protocol}_volume_{volume:04d}_group_{group:04d}.txt",
        )
        with open(path, "w", encoding="utf-8") as pointer_file:
            pointer_file.write(contents)
        return path

    def _write_reference(self, value: float = 1.0) -> str:
        name = f"{self.protocol}_volume_0000_20260916T120000.nhdr"
        path = os.path.join(self.input_dir, name)
        array = np.full((3, 3, 4), value, dtype=np.float32)
        image = sitk.GetImageFromArray(array)
        image.SetSpacing((2.0, 2.0, 3.0))
        sitk.WriteImage(image, path)
        return self._pointer(0, contents=f"{name}\n")

    def _write_slice(
        self,
        volume: int,
        index: int,
        *,
        value: float | None = None,
        center_z: float | None = None,
    ) -> str:
        path = os.path.join(
            self.input_dir,
            f"{self.protocol}_volume_{volume:04d}_slice_{index:04d}.nhdr",
        )
        plane_value = float(index if value is None else value)
        array = np.full((1, 3, 4), plane_value, dtype=np.float32)
        image = sitk.GetImageFromArray(array)
        image.SetSpacing((2.0, 2.0, 3.0))
        image.SetOrigin((0.0, 0.0, float(index if center_z is None else center_z)))
        sitk.WriteImage(image, path)
        return path

    def _write_complete_volume(self, volume: int, value: float | None = None) -> str:
        for index in range(3):
            self._write_slice(volume, index, value=value)
        return self._pointer(volume)

    def _motion_ready(self, processor: TSNRVolumeProcessor, *volumes: int) -> None:
        status_path = os.path.join(self.input_dir, "registration_status.json")
        document = {
            "schema_version": 1,
            "revision": len(volumes),
            "volumes": {
                str(volume): {
                    "registration_finalized": True,
                    "groups": {
                        "0": {
                            "status": "registered",
                            "transform_filename": f"alignTransform_{volume:04d}.tfm",
                        }
                    },
                }
                for volume in volumes
            },
        }
        with open(status_path, "w", encoding="utf-8") as status_file:
            json.dump(document, status_file)
        self.assertTrue(processor.handle_registration_status(status_path))

    def test_incomplete_volume_remains_pending_then_processes_exactly_once(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(0.0), expected_slice_count=None)
        pointer = self._pointer(1)
        self._write_slice(1, 0, center_z=30.0)
        self._write_slice(1, 1, center_z=10.0)
        self._motion_ready(processor, 1)

        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.pending_volumes, {1})
        self.assertEqual(processor.processed_volumes, {0})
        self.assertEqual(processor.accumulator.count, 1)

        self._write_slice(1, 2, center_z=20.0)
        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.pending_volumes, set())
        self.assertEqual(processor.processed_volumes, {0, 1})
        self.assertEqual(processor.accumulator.count, 2)
        np.testing.assert_array_equal(
            processor.accumulator.mean[:, 0, 0], np.array([0.5, 1.0, 0.0])
        )

        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.accumulator.count, 2)

    def test_exact_slice_membership_is_required(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer = self._pointer(2)
        self._motion_ready(processor, 2)
        for index in (0, 1, 3):
            self._write_slice(2, index)

        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.pending_volumes, {2})
        self.assertEqual(processor.accumulator.count, 1)

    def test_read_failure_does_not_mark_volume_processed(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer = self._pointer(3)
        self._motion_ready(processor, 3)
        for index in (0, 1):
            self._write_slice(3, index)
        bad_path = os.path.join(
            self.input_dir,
            f"{self.protocol}_volume_0003_slice_0002.nhdr",
        )
        with open(bad_path, "w", encoding="utf-8") as bad_slice:
            bad_slice.write("not a valid NRRD\n")

        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.pending_volumes, {3})
        self.assertEqual(processor.processed_volumes, {0})
        self.assertEqual(processor.accumulator.count, 1)

    def test_accumulator_failure_retains_pending_volume(self) -> None:
        processor = self._processor(accumulator=FailingAccumulator())
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer = self._write_complete_volume(4)
        self._motion_ready(processor, 4)

        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.pending_volumes, {4})
        self.assertEqual(processor.processed_volumes, {0})

    def test_display_failure_does_not_make_volume_eligible_again(self) -> None:
        def fail_display(*args, **kwargs):
            raise RuntimeError("injected display failure")

        processor = self._processor(min_samples=2, mosaic_factory=fail_display)
        reference_pointer = self._write_reference(1.0)
        processor.handle_pointer(reference_pointer, expected_slice_count=None)
        volume_pointer = self._write_complete_volume(1, value=2.0)
        self._motion_ready(processor, 1)
        processor.handle_pointer(volume_pointer, expected_slice_count=3)

        self.assertEqual(processor.accumulator.count, 2)
        self.assertEqual(processor.processed_volumes, {0, 1})
        self.assertEqual(processor.pending_volumes, set())
        processor.handle_pointer(volume_pointer, expected_slice_count=3)
        self.assertEqual(processor.accumulator.count, 2)

    def test_display_starts_at_sample_twenty_with_slice_dimensions(self) -> None:
        saved_mosaics: list[np.ndarray] = []

        def capture_image(mosaic, output_path):
            saved_mosaics.append(mosaic.copy())

        processor = self._processor(min_samples=20, image_saver=capture_image)
        processor.handle_pointer(self._write_reference(1.0), expected_slice_count=None)
        for volume in range(1, 19):
            pointer = self._write_complete_volume(volume, value=float(volume + 1))
            self._motion_ready(processor, volume)
            processor.handle_pointer(pointer, expected_slice_count=3)

        self.assertEqual(processor.accumulator.count, 19)
        self.assertEqual(saved_mosaics, [])

        twentieth_pointer = self._write_complete_volume(19, value=20.0)
        self._motion_ready(processor, 19)
        processor.handle_pointer(twentieth_pointer, expected_slice_count=3)
        self.assertEqual(processor.accumulator.count, 20)
        self.assertEqual(len(saved_mosaics), 1)
        self.assertEqual(saved_mosaics[0].shape, (3, 4))

        twenty_first_pointer = self._write_complete_volume(20, value=21.0)
        self._motion_ready(processor, 20)
        processor.handle_pointer(twenty_first_pointer, expected_slice_count=3)
        self.assertEqual(processor.accumulator.count, 21)
        self.assertEqual(len(saved_mosaics), 2)
        self.assertEqual(saved_mosaics[1].shape, (3, 4))

    def test_reference_is_accumulated_exactly_once(self) -> None:
        processor = self._processor()
        pointer = self._write_reference()
        processor.handle_pointer(pointer, expected_slice_count=None)
        processor.handle_pointer(pointer, expected_slice_count=None)
        self.assertEqual(processor.accumulator.count, 1)
        self.assertEqual(processor.processed_volumes, {0})

    def test_idle_retry_is_throttled_and_reset_preserves_jpeg(self) -> None:
        processor = self._processor(idle_retry_interval=0.5)
        processor.pending_volumes.add(9)
        with open(self.output_path, "wb") as output_file:
            output_file.write(b"previous JPEG")

        self.assertFalse(processor.retry_pending_if_due(3, now=0.1))
        self.assertTrue(processor.retry_pending_if_due(3, now=0.5))
        self.assertFalse(processor.retry_pending_if_due(3, now=0.6))

        processor.accumulator.update(np.ones((3, 3, 4)))
        processor.processed_volumes.add(0)
        processor.reset()
        self.assertEqual(processor.accumulator.count, 0)
        self.assertEqual(processor.pending_volumes, set())
        self.assertEqual(processor.processed_volumes, set())
        self.assertEqual(processor.image_ready_volumes, {})
        self.assertEqual(processor.motion_ready_volumes, {})
        self.assertTrue(os.path.exists(self.output_path))

    def test_image_ready_waits_for_motion_ready(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer = self._write_complete_volume(1, value=2.0)

        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.accumulator.count, 1)
        self.assertIn(1, processor.image_ready_volumes)
        self.assertIn(1, processor.pending_volumes)

        self._motion_ready(processor, 1)
        self.assertEqual(processor.accumulator.count, 2)
        self.assertEqual(processor.processed_volumes, {0, 1})
        self.assertNotIn(1, processor.image_ready_volumes)

    def test_motion_ready_before_image_ready_and_duplicate_reads_are_idempotent(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        self._motion_ready(processor, 2)
        self._motion_ready(processor, 2)
        self.assertEqual(processor.accumulator.count, 1)

        pointer = self._write_complete_volume(2, value=3.0)
        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.accumulator.count, 2)
        self._motion_ready(processor, 2)
        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.accumulator.count, 2)

    def test_one_status_read_discovers_multiple_finalized_volumes(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer1 = self._write_complete_volume(1, value=2.0)
        pointer2 = self._write_complete_volume(2, value=3.0)
        processor.handle_pointer(pointer1, expected_slice_count=3)
        processor.handle_pointer(pointer2, expected_slice_count=3)

        self._motion_ready(processor, 1, 2)
        self.assertEqual(processor.accumulator.count, 3)
        self.assertEqual(processor.processed_volumes, {0, 1, 2})

    def test_malformed_registration_status_is_nonfatal(self) -> None:
        processor = self._processor()
        status_path = os.path.join(self.input_dir, "registration_status.json")
        with open(status_path, "w", encoding="utf-8") as status_file:
            status_file.write("{not valid JSON")

        self.assertFalse(processor.handle_registration_status(status_path))
        self.assertEqual(processor.accumulator.count, 0)
        self.assertEqual(processor.motion_ready_volumes, {})

    def test_successful_nonreference_volume_retires_only_image_work_files(self) -> None:
        processor = self._processor()
        reference_pointer = self._write_reference()
        with open(reference_pointer, "r", encoding="utf-8") as pointer_file:
            reference_header_name = pointer_file.read().strip()
        reference_raw_name = os.path.splitext(reference_header_name)[0] + ".raw"
        processor.handle_pointer(reference_pointer, expected_slice_count=None)

        volume_pointer = self._write_complete_volume(1, value=2.0)
        second_pointer = self._pointer(1, group=1)
        transform_name = "alignTransform_0001_0001-0000.tfm"
        unrelated_name = f"{self.protocol}_volume_0001_notes.txt"
        upsampled_header_name = (
            f"{self.protocol}_volume_0000_20260916T120000_upsampled.nhdr"
        )
        upsampled_raw_name = (
            f"{self.protocol}_volume_0000_20260916T120000_upsampled.raw"
        )
        for name in (
            transform_name,
            unrelated_name,
            upsampled_header_name,
            upsampled_raw_name,
        ):
            with open(os.path.join(self.input_dir, name), "wb") as output_file:
                output_file.write(b"leave in acquisition root")

        expected_retired = {
            os.path.basename(volume_pointer),
            os.path.basename(second_pointer),
        }
        for index in range(3):
            stem = f"{self.protocol}_volume_0001_slice_{index:04d}"
            expected_retired.update({f"{stem}.nhdr", f"{stem}.raw"})

        self._motion_ready(processor, 1)
        processor.handle_pointer(volume_pointer, expected_slice_count=3)

        processed_dir = os.path.join(
            self.input_dir, f"processed_{self.protocol}"
        )
        self.assertTrue(os.path.isdir(processed_dir))
        self.assertEqual(set(os.listdir(processed_dir)), expected_retired)
        retired_image = sitk.ReadImage(
            os.path.join(
                processed_dir,
                f"{self.protocol}_volume_0001_slice_0000.nhdr",
            )
        )
        self.assertEqual(retired_image.GetSize(), (4, 3, 1))
        for name in expected_retired:
            self.assertFalse(os.path.exists(os.path.join(self.input_dir, name)))
        for name in (
            transform_name,
            unrelated_name,
            upsampled_header_name,
            upsampled_raw_name,
            os.path.basename(reference_pointer),
            reference_header_name,
            reference_raw_name,
        ):
            self.assertTrue(os.path.exists(os.path.join(self.input_dir, name)), name)

    def test_retirement_waits_for_successful_tsnr_incorporation(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer = self._write_complete_volume(1, value=2.0)

        processor.handle_pointer(pointer, expected_slice_count=3)

        self.assertIn(1, processor.image_ready_volumes)
        self.assertIn(1, processor.pending_volumes)
        self.assertTrue(os.path.exists(pointer))
        self.assertFalse(
            os.path.exists(os.path.join(self.input_dir, f"processed_{self.protocol}"))
        )

    def test_pending_load_and_accumulator_failure_leave_files_active(self) -> None:
        pending_processor = self._processor()
        pending_processor.handle_pointer(
            self._write_reference(), expected_slice_count=None
        )
        pending_pointer = self._pointer(3)
        pending_header = self._write_slice(3, 0)
        pending_raw = os.path.splitext(pending_header)[0] + ".raw"
        self._motion_ready(pending_processor, 3)
        pending_processor.handle_pointer(pending_pointer, expected_slice_count=3)

        for path in (pending_pointer, pending_header, pending_raw):
            self.assertTrue(os.path.exists(path))
        self.assertFalse(
            os.path.exists(os.path.join(self.input_dir, f"processed_{self.protocol}"))
        )

        failing_processor = self._processor(accumulator=FailingAccumulator())
        failing_processor.handle_pointer(
            self._pointer(
                0,
                contents=f"{self.protocol}_volume_0000_20260916T120000.nhdr\n",
            ),
            expected_slice_count=None,
        )
        failing_pointer = self._write_complete_volume(4, value=4.0)
        self._motion_ready(failing_processor, 4)
        failing_processor.handle_pointer(failing_pointer, expected_slice_count=3)

        self.assertIn(4, failing_processor.pending_volumes)
        self.assertNotIn(4, failing_processor.processed_volumes)
        self.assertTrue(os.path.exists(failing_pointer))

    def test_destination_collision_preserves_existing_and_source_pair(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer = self._write_complete_volume(5, value=5.0)
        collided_header_name = f"{self.protocol}_volume_0005_slice_0001.nhdr"
        collided_raw_name = f"{self.protocol}_volume_0005_slice_0001.raw"
        processed_dir = os.path.join(
            self.input_dir, f"processed_{self.protocol}"
        )
        os.makedirs(processed_dir)
        collided_destination = os.path.join(processed_dir, collided_header_name)
        with open(collided_destination, "wb") as destination_file:
            destination_file.write(b"existing acquisition data")

        self._motion_ready(processor, 5)
        processor.handle_pointer(pointer, expected_slice_count=3)

        with open(collided_destination, "rb") as destination_file:
            self.assertEqual(destination_file.read(), b"existing acquisition data")
        self.assertTrue(os.path.exists(os.path.join(self.input_dir, collided_header_name)))
        self.assertTrue(os.path.exists(os.path.join(self.input_dir, collided_raw_name)))

    def test_retirement_failure_does_not_revert_successful_tsnr_state(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer = self._write_complete_volume(6, value=6.0)
        self._motion_ready(processor, 6)

        with patch(
            "apps.motion_monitor.tsnr_integration.os.rename",
            side_effect=OSError("injected rename failure"),
        ):
            processor.handle_pointer(pointer, expected_slice_count=3)

        self.assertEqual(processor.accumulator.count, 2)
        self.assertIn(6, processor.processed_volumes)
        self.assertNotIn(6, processor.pending_volumes)
        self.assertTrue(os.path.exists(pointer))

    def test_retirement_directory_failure_leaves_files_active(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer = self._write_complete_volume(7, value=7.0)
        processed_path = os.path.join(
            self.input_dir, f"processed_{self.protocol}"
        )
        with open(processed_path, "wb") as blocking_file:
            blocking_file.write(b"not a directory")
        self._motion_ready(processor, 7)

        processor.handle_pointer(pointer, expected_slice_count=3)

        self.assertIn(7, processor.processed_volumes)
        self.assertNotIn(7, processor.pending_volumes)
        self.assertTrue(os.path.isfile(processed_path))
        self.assertTrue(os.path.exists(pointer))

    def test_multiple_volumes_retire_sequentially(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)

        retired_names: set[str] = set()
        for volume in (1, 2):
            pointer = self._write_complete_volume(volume, value=float(volume + 1))
            self._motion_ready(processor, volume)
            processor.handle_pointer(pointer, expected_slice_count=3)
            retired_names.add(os.path.basename(pointer))
            for index in range(3):
                stem = f"{self.protocol}_volume_{volume:04d}_slice_{index:04d}"
                retired_names.update({f"{stem}.nhdr", f"{stem}.raw"})

        processed_dir = os.path.join(
            self.input_dir, f"processed_{self.protocol}"
        )
        self.assertTrue(retired_names.issubset(set(os.listdir(processed_dir))))
        self.assertEqual(processor.processed_volumes, {0, 1, 2})
        self.assertEqual(processor.accumulator.count, 3)


if __name__ == "__main__":
    unittest.main()

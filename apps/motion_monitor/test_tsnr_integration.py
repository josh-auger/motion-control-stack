"""Focused tests for the active motion-monitor TSNR orchestration."""

import logging
import os
import tempfile
import unittest

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

    def test_incomplete_volume_remains_pending_then_processes_exactly_once(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(0.0), expected_slice_count=None)
        pointer = self._pointer(1)
        self._write_slice(1, 0, center_z=30.0)
        self._write_slice(1, 1, center_z=10.0)

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
        for index in (0, 1, 3):
            self._write_slice(2, index)

        processor.handle_pointer(pointer, expected_slice_count=3)
        self.assertEqual(processor.pending_volumes, {2})
        self.assertEqual(processor.accumulator.count, 1)

    def test_read_failure_does_not_mark_volume_processed(self) -> None:
        processor = self._processor()
        processor.handle_pointer(self._write_reference(), expected_slice_count=None)
        pointer = self._pointer(3)
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
            processor.handle_pointer(pointer, expected_slice_count=3)

        self.assertEqual(processor.accumulator.count, 19)
        self.assertEqual(saved_mosaics, [])

        twentieth_pointer = self._write_complete_volume(19, value=20.0)
        processor.handle_pointer(twentieth_pointer, expected_slice_count=3)
        self.assertEqual(processor.accumulator.count, 20)
        self.assertEqual(len(saved_mosaics), 1)
        self.assertEqual(saved_mosaics[0].shape, (3, 4))

        twenty_first_pointer = self._write_complete_volume(20, value=21.0)
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
        self.assertTrue(os.path.exists(self.output_path))


if __name__ == "__main__":
    unittest.main()

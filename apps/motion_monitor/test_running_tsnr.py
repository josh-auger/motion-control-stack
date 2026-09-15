"""Lightweight tests for the real-time running TSNR helpers."""

import os
import tempfile
import unittest

import cv2
import numpy as np

from apps.motion_monitor.running_tsnr import (
    RunningTSNR,
    atomic_save_tsnr_image,
    create_tsnr_mosaic,
)


class RunningTSNRTests(unittest.TestCase):
    def setUp(self) -> None:
        self.stack = np.array(
            [
                np.full((3, 4, 5), 2.0),
                np.full((3, 4, 5), 4.0),
                np.full((3, 4, 5), 8.0),
                np.full((3, 4, 5), 10.0),
            ]
        )
        self.stack[:, 0, 0, 0] = 7.0

    def test_running_statistics_match_numpy(self) -> None:
        running = RunningTSNR()
        for volume in self.stack:
            running.update(volume)

        expected_mean = np.mean(self.stack, axis=0)
        expected_std = np.std(self.stack, axis=0, ddof=1)
        expected_tsnr = np.divide(
            expected_mean,
            expected_std,
            out=np.zeros_like(expected_mean),
            where=expected_std > 0,
        )

        np.testing.assert_allclose(running.mean, expected_mean)
        np.testing.assert_allclose(running.get_std(), expected_std)
        np.testing.assert_allclose(running.get_tsnr(), expected_tsnr)
        self.assertEqual(running.get_tsnr()[0, 0, 0], 0.0)

    def test_reset_clears_state(self) -> None:
        running = RunningTSNR()
        running.update(self.stack[0])
        running.reset()
        self.assertEqual(running.count, 0)
        self.assertIsNone(running.mean)
        self.assertIsNone(running.M2)

    def test_shape_mismatch_raises_clear_error(self) -> None:
        running = RunningTSNR()
        running.update(np.zeros((2, 3, 4)))
        with self.assertRaisesRegex(ValueError, "does not match initial shape"):
            running.update(np.zeros((2, 3, 5)))

    def test_first_volume_has_zero_tsnr(self) -> None:
        running = RunningTSNR()
        running.update(self.stack[0])
        np.testing.assert_array_equal(running.get_tsnr(), np.zeros_like(self.stack[0]))

    def test_invalid_statistics_and_display_values_become_zero(self) -> None:
        running = RunningTSNR()
        first = self.stack[0].copy()
        second = self.stack[1].copy()
        first[0, 0, 0] = np.nan
        second[0, 0, 1] = np.inf
        running.update(first)
        running.update(second)
        tsnr = running.get_tsnr()
        self.assertTrue(np.isfinite(tsnr).all())
        self.assertEqual(tsnr[0, 0, 0], 0.0)
        self.assertEqual(tsnr[0, 0, 1], 0.0)

        tsnr[0, 0, 0] = np.nan
        tsnr[1, 0, 0] = np.inf
        mosaic = create_tsnr_mosaic(tsnr, output_shape=(17, 19))
        self.assertEqual(mosaic.dtype, np.uint8)
        self.assertTrue(np.isfinite(mosaic).all())

    def test_mosaic_has_exact_requested_shape(self) -> None:
        volume = np.arange(11 * 12 * 13, dtype=np.float64).reshape(11, 12, 13)
        mosaic = create_tsnr_mosaic(volume, output_shape=(176, 175), display_max=100.0)
        self.assertEqual(mosaic.shape, (176, 175))

    def test_mosaic_uses_fixed_normalized_axial_locations(self) -> None:
        volume = np.broadcast_to(np.arange(10)[:, None, None], (10, 4, 4))
        mosaic = create_tsnr_mosaic(volume, output_shape=(4, 4), display_max=10.0)
        expected = np.array(
            [
                [51, 51, 102, 102],
                [51, 51, 102, 102],
                [153, 153, 204, 204],
                [153, 153, 204, 204],
            ],
            dtype=np.uint8,
        )
        np.testing.assert_array_equal(mosaic, expected)

    def test_atomic_save_produces_readable_jpeg(self) -> None:
        mosaic = create_tsnr_mosaic(self.stack[0], output_shape=(31, 29))
        with tempfile.TemporaryDirectory() as directory:
            output_path = os.path.join(directory, "tsnr_dashboard.jpg")
            atomic_save_tsnr_image(mosaic, output_path)
            loaded = cv2.imread(output_path, cv2.IMREAD_GRAYSCALE)
            self.assertIsNotNone(loaded)
            self.assertEqual(loaded.shape, mosaic.shape)
            self.assertFalse(any(name.endswith(".tmp.jpg") for name in os.listdir(directory)))


if __name__ == "__main__":
    unittest.main()

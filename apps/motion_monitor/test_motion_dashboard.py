"""Focused tests for the combined motion and running-TSNR dashboard."""

import os
import tempfile
import unittest

import cv2
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from apps.motion_monitor.generate_motion_plots import (
    DASHBOARD_PIXEL_HEIGHT,
    DASHBOARD_PIXEL_WIDTH,
    add_tsnr_dashboard_panel,
    plot_motion_dashboard,
)
from apps.motion_monitor.running_tsnr import RunningTSNR


class MotionDashboardTests(unittest.TestCase):
    def setUp(self) -> None:
        count = 24
        self.motion_df = pd.DataFrame(
            {
                "X_rotation(rad)": np.linspace(0.0, 0.02, count),
                "Y_rotation(rad)": np.linspace(0.0, -0.01, count),
                "Z_rotation(rad)": np.zeros(count),
                "X_translation(mm)": np.linspace(0.0, 1.0, count),
                "Y_translation(mm)": np.linspace(0.0, -0.5, count),
                "Z_translation(mm)": np.zeros(count),
                "Displacement(mm)": np.linspace(0.0, 0.6, count),
                "Cumulative_displacement(mm)": np.linspace(0.0, 4.0, count),
                "Volume_index": np.arange(1, count + 1),
                "Motion_flag": np.zeros(count, dtype=int),
            }
        )

    @staticmethod
    def _running_tsnr(sample_count: int) -> RunningTSNR:
        accumulator = RunningTSNR()
        base = np.arange(10 * 12 * 12, dtype=np.float64).reshape(10, 12, 12)
        for sample in range(sample_count):
            accumulator.update(base + sample)
        return accumulator

    def _render_dashboard(self, sample_count: int) -> np.ndarray:
        accumulator = self._running_tsnr(sample_count)
        tsnr_volume = (
            accumulator.get_tsnr()
            if accumulator.count >= 20
            else None
        )
        with tempfile.TemporaryDirectory() as directory:
            output_path = os.path.join(directory, "dashboard.jpg")
            plot_motion_dashboard(
                self.motion_df,
                output_filename=output_path,
                protocol_name="synthetic_test",
                threshold=0.3,
                num_expected_volumes=30,
                num_moved_volumes=0,
                livestream_enabled=False,
                tsnr_volume=tsnr_volume,
                tsnr_count=accumulator.count,
                tsnr_min_samples=20,
                tsnr_display_max=100.0,
            )
            image = cv2.imread(output_path)
            self.assertIsNotNone(image)
            return image

    def test_dashboard_before_threshold_is_valid_and_shows_waiting_state(self) -> None:
        figure = plt.figure()
        panel = add_tsnr_dashboard_panel(
            figure,
            figure.add_gridspec(1, 1)[0, 0],
            tsnr_count=19,
            min_samples=20,
        )
        self.assertEqual(panel["state"], "waiting")
        self.assertEqual(panel["text"].get_text(), "Waiting for TSNR\nN = 19 / 20")
        plt.close(figure)

        image = self._render_dashboard(19)
        self.assertEqual(
            image.shape[:2],
            (DASHBOARD_PIXEL_HEIGHT, DASHBOARD_PIXEL_WIDTH),
        )

    def test_dashboard_at_sample_twenty_uses_actual_count_and_fixed_scale(self) -> None:
        accumulator = self._running_tsnr(20)
        figure = plt.figure(figsize=(5, 4))
        panel = add_tsnr_dashboard_panel(
            figure,
            figure.add_gridspec(1, 1)[0, 0],
            tsnr_volume=accumulator.get_tsnr(),
            tsnr_count=accumulator.count,
            min_samples=20,
            display_max=100.0,
        )

        self.assertEqual(panel["state"], "ready")
        self.assertIn("N = 20", panel["overlay"].get_text())
        self.assertTrue(np.isfinite(panel["mean_tsnr"]))
        self.assertEqual(panel["image_artist"].get_clim(), (0.0, 100.0))
        self.assertEqual(len(panel["image_axes"]), 4)
        plt.close(figure)

        image = self._render_dashboard(20)
        self.assertEqual(
            image.shape[:2],
            (DASHBOARD_PIXEL_HEIGHT, DASHBOARD_PIXEL_WIDTH),
        )

    def test_dashboard_above_threshold_remains_valid(self) -> None:
        image = self._render_dashboard(21)
        self.assertEqual(
            image.shape[:2],
            (DASHBOARD_PIXEL_HEIGHT, DASHBOARD_PIXEL_WIDTH),
        )

    def test_tsnr_panel_failure_does_not_prevent_motion_dashboard(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output_path = os.path.join(directory, "dashboard.jpg")
            with self.assertLogs(level="ERROR"):
                plot_motion_dashboard(
                    self.motion_df,
                    output_filename=output_path,
                    protocol_name="synthetic_test",
                    threshold=0.3,
                    num_expected_volumes=30,
                    num_moved_volumes=0,
                    tsnr_volume=np.zeros((2, 2)),
                    tsnr_count=20,
                    tsnr_min_samples=20,
                )
            image = cv2.imread(output_path)
            self.assertIsNotNone(image)
            self.assertEqual(
                image.shape[:2],
                (DASHBOARD_PIXEL_HEIGHT, DASHBOARD_PIXEL_WIDTH),
            )


if __name__ == "__main__":
    unittest.main()

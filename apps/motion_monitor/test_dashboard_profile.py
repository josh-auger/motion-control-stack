"""Focused tests for observational motion-dashboard profiling."""

import csv
import os
import tempfile
import unittest
from unittest import mock

import cv2
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from apps.motion_monitor.dashboard_profile import (
    DashboardProfiler,
    dashboard_profiling_enabled,
)
from apps.motion_monitor.generate_motion_plots import (
    DASHBOARD_PIXEL_HEIGHT,
    DASHBOARD_PIXEL_WIDTH,
    plot_motion_dashboard,
)


class DashboardProfileTests(unittest.TestCase):
    @staticmethod
    def motion_dataframe():
        count = 40
        return pd.DataFrame(
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

    def test_profiled_render_preserves_output_and_closes_figure(self):
        motion_df = self.motion_dataframe()
        tsnr_volume = np.ones((10, 12, 12), dtype=float)
        with tempfile.TemporaryDirectory() as directory:
            baseline_path = os.path.join(directory, "baseline.jpg")
            profiled_path = os.path.join(directory, "profiled.jpg")
            baseline_figures = tuple(plt.get_fignums())

            plot_motion_dashboard(
                motion_df,
                output_filename=baseline_path,
                protocol_name="profile_test",
                threshold=0.3,
                num_expected_volumes=50,
                num_moved_volumes=0,
                tsnr_volume=tsnr_volume,
                tsnr_count=20,
                tsnr_min_samples=10,
            )
            timings = {}
            plot_motion_dashboard(
                motion_df,
                output_filename=profiled_path,
                protocol_name="profile_test",
                threshold=0.3,
                num_expected_volumes=50,
                num_moved_volumes=0,
                tsnr_volume=tsnr_volume,
                tsnr_count=20,
                tsnr_min_samples=10,
                profile_timings=timings,
            )

            baseline = cv2.imread(baseline_path)
            profiled = cv2.imread(profiled_path)
            self.assertIsNotNone(baseline)
            self.assertIsNotNone(profiled)
            self.assertEqual(
                profiled.shape[:2],
                (DASHBOARD_PIXEL_HEIGHT, DASHBOARD_PIXEL_WIDTH),
            )
            np.testing.assert_array_equal(profiled, baseline)
            self.assertEqual(tuple(plt.get_fignums()), baseline_figures)

        for field in (
            "plot_data_extract_ms",
            "figure_create_ms",
            "plot_population_ms",
            "annotation_ms",
            "tsnr_panel_ms",
            "layout_ms",
            "savefig_ms",
            "cleanup_ms",
            "plot_total_ms",
        ):
            self.assertGreaterEqual(timings[field], 0.0)

    def test_buffered_profiler_writes_component_accounting(self):
        with tempfile.TemporaryDirectory() as directory:
            profiler = DashboardProfiler(directory)
            profiler.record(
                {
                    "dashboard_index": 1,
                    "outcome": "published",
                    "monitor_data_prep_ms": 2.0,
                    "savefig_ms": 10.0,
                    "atomic_replace_ms": 0.5,
                    "plot_total_ms": 11.0,
                    "total_ms": 15.0,
                }
            )
            path = profiler.path
            profiler.close()

            with open(path, newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["outcome"], "published")
            self.assertAlmostEqual(float(rows[0]["component_sum_ms"]), 12.5)
            self.assertAlmostEqual(float(rows[0]["unattributed_ms"]), 2.5)

    def test_profile_flag_is_opt_in(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertFalse(dashboard_profiling_enabled(None))
        self.assertFalse(dashboard_profiling_enabled("off"))
        self.assertTrue(dashboard_profiling_enabled("ON"))


if __name__ == "__main__":
    unittest.main()

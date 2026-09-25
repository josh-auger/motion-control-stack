"""Opt-in, acquisition-local profiling for motion-dashboard publication.

The profiler is observational only: failures disable profiling and never alter
dashboard generation.  Phase values use ``time.perf_counter_ns()`` and are
reported in milliseconds.  ``plot_total_ms`` overlaps the detailed plotting
phases and is therefore excluded from ``component_sum_ms``.
"""

import csv
import logging
import os
from datetime import datetime


PROFILE_FIELDS = (
    "dashboard_index",
    "trigger_transform_index",
    "volume",
    "group",
    "motion_sample_count",
    "tsnr_count",
    "outcome",
    "trigger_setup_ms",
    "monitor_data_prep_ms",
    "tsnr_prepare_ms",
    "plot_data_extract_ms",
    "figure_create_ms",
    "plot_population_ms",
    "annotation_ms",
    "tsnr_panel_ms",
    "layout_ms",
    "savefig_ms",
    "cleanup_ms",
    "plot_total_ms",
    "output_validate_ms",
    "atomic_replace_ms",
    "published_read_ms",
    "stream_push_ms",
    "component_sum_ms",
    "unattributed_ms",
    "total_ms",
    "dashboard_path",
)

COMPONENT_FIELDS = (
    "trigger_setup_ms",
    "monitor_data_prep_ms",
    "tsnr_prepare_ms",
    "plot_data_extract_ms",
    "figure_create_ms",
    "plot_population_ms",
    "annotation_ms",
    "tsnr_panel_ms",
    "layout_ms",
    "savefig_ms",
    "cleanup_ms",
    "output_validate_ms",
    "atomic_replace_ms",
    "published_read_ms",
    "stream_push_ms",
)


def dashboard_profiling_enabled(value=None):
    """Use the project's existing on/off spelling for the diagnostic flag."""
    if value is None:
        value = os.environ.get("MOTION_DASHBOARD_PROFILE_FLAG", "off")
    return str(value).lower() == "on"


class DashboardProfiler:
    """Write one buffered CSV row for each attempted dashboard publication."""

    def __init__(self, directory):
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        self.path = os.path.join(
            directory,
            f"motion_dashboard_profile_{stamp}.csv",
        )
        self.stream = open(
            self.path,
            "w",
            newline="",
            buffering=1024 * 1024,
        )
        self.writer = csv.DictWriter(self.stream, fieldnames=PROFILE_FIELDS)
        self.writer.writeheader()
        self.count = 0
        self.disabled = False

    def record(self, values):
        if self.disabled:
            return
        row = dict(values)
        component_sum = sum(
            float(row.get(field, 0.0) or 0.0)
            for field in COMPONENT_FIELDS
        )
        total = float(row.get("total_ms", 0.0) or 0.0)
        row["component_sum_ms"] = component_sum
        row["unattributed_ms"] = total - component_sum
        try:
            self.writer.writerow(row)
            self.count += 1
            if self.count % 64 == 0:
                self.stream.flush()
        except (OSError, ValueError) as error:
            self.disabled = True
            logging.warning(
                "Motion-dashboard profiling disabled after write failure: %s",
                error,
            )

    def close(self):
        try:
            self.stream.flush()
            self.stream.close()
        except (OSError, ValueError) as error:
            logging.warning("Failed to close motion-dashboard profile: %s", error)

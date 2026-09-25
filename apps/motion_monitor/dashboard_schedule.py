"""Acquisition-local monotonic scheduling for periodic motion dashboards."""

import logging
import math
import os
import time


DEFAULT_DASHBOARD_INTERVAL_SEC = 5.0
DASHBOARD_INTERVAL_ENV = "MOTION_DASHBOARD_INTERVAL_SEC"


def load_dashboard_interval(value=None):
    """Return a finite positive dashboard interval, falling back safely."""
    raw_value = os.environ.get(DASHBOARD_INTERVAL_ENV) if value is None else value
    if raw_value is None:
        return DEFAULT_DASHBOARD_INTERVAL_SEC

    try:
        interval = float(raw_value)
    except (TypeError, ValueError):
        interval = None

    if interval is None or not math.isfinite(interval) or interval <= 0.0:
        logging.warning(
            "Invalid %s=%r; using default %.3f s.",
            DASHBOARD_INTERVAL_ENV,
            raw_value,
            DEFAULT_DASHBOARD_INTERVAL_SEC,
        )
        return DEFAULT_DASHBOARD_INTERVAL_SEC
    return interval


class PeriodicDashboardSchedule:
    """Run at most one dashboard after each completed monotonic interval."""

    def __init__(self, interval_sec, clock=time.monotonic):
        self.interval_sec = interval_sec
        self.clock = clock
        self.last_completed = None

    def start(self):
        """Start a fresh acquisition interval without generating immediately."""
        self.last_completed = self.clock()

    def reset(self):
        """Discard timing state at the acquisition boundary."""
        self.last_completed = None

    def due(self):
        if self.last_completed is None:
            return False
        return self.clock() - self.last_completed >= self.interval_sec

    def run_if_due(self, generate_dashboard):
        """Generate once when due and restart timing after successful completion."""
        if not self.due():
            return False
        if not generate_dashboard():
            return False
        self.last_completed = self.clock()
        return True


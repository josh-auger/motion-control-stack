"""Focused tests for monotonic periodic-dashboard scheduling."""

import os
import unittest
from unittest import mock

from apps.motion_monitor.dashboard_schedule import (
    DEFAULT_DASHBOARD_INTERVAL_SEC,
    PeriodicDashboardSchedule,
    load_dashboard_interval,
)


class FakeClock:
    def __init__(self, value=0.0):
        self.value = value

    def __call__(self):
        return self.value

    def advance(self, seconds):
        self.value += seconds


class DashboardScheduleTests(unittest.TestCase):
    def setUp(self):
        self.clock = FakeClock(100.0)
        self.schedule = PeriodicDashboardSchedule(5.0, clock=self.clock)
        self.schedule.start()

    def test_interval_not_reached_does_not_generate(self):
        generate = mock.Mock(return_value=True)
        self.clock.advance(4.999)

        self.assertFalse(self.schedule.run_if_due(generate))
        generate.assert_not_called()

    def test_interval_reached_generates_exactly_once(self):
        generate = mock.Mock(return_value=True)
        self.clock.advance(5.0)

        self.assertTrue(self.schedule.run_if_due(generate))
        generate.assert_called_once_with()

    def test_timer_restarts_after_completed_generation(self):
        def generate():
            self.clock.advance(0.5)
            return True

        self.clock.advance(5.0)
        self.assertTrue(self.schedule.run_if_due(generate))
        self.clock.advance(4.9)
        self.assertFalse(self.schedule.run_if_due(mock.Mock(return_value=True)))
        self.clock.advance(0.1)
        self.assertTrue(self.schedule.due())

    def test_large_elapsed_interval_does_not_catch_up(self):
        generate = mock.Mock(return_value=True)
        self.clock.advance(20.0)

        self.assertTrue(self.schedule.run_if_due(generate))
        self.assertFalse(self.schedule.run_if_due(generate))
        generate.assert_called_once_with()

    def test_failed_generation_does_not_advance_timer(self):
        self.clock.advance(5.0)
        self.assertFalse(self.schedule.run_if_due(lambda: None))
        self.assertTrue(self.schedule.due())

    def test_reset_gives_next_acquisition_a_fresh_interval(self):
        self.clock.advance(20.0)
        self.assertTrue(self.schedule.due())
        self.schedule.reset()
        self.clock.advance(100.0)
        self.assertFalse(self.schedule.due())

        self.schedule.start()
        self.clock.advance(4.9)
        self.assertFalse(self.schedule.due())
        self.clock.advance(0.1)
        self.assertTrue(self.schedule.due())


class DashboardIntervalConfigurationTests(unittest.TestCase):
    def test_default_interval(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(
                load_dashboard_interval(),
                DEFAULT_DASHBOARD_INTERVAL_SEC,
            )

    def test_valid_custom_interval(self):
        with mock.patch.dict(
            os.environ,
            {"MOTION_DASHBOARD_INTERVAL_SEC": "3.0"},
            clear=True,
        ):
            self.assertEqual(load_dashboard_interval(), 3.0)

    def test_invalid_intervals_fall_back_with_warning(self):
        for value in ("abc", "0", "-1", "nan", "inf"):
            with self.subTest(value=value), self.assertLogs(level="WARNING"):
                self.assertEqual(
                    load_dashboard_interval(value),
                    DEFAULT_DASHBOARD_INTERVAL_SEC,
                )


if __name__ == "__main__":
    unittest.main()


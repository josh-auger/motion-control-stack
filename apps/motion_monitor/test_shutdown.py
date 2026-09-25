"""Focused motion-monitor DONE acknowledgement tests."""

from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import monitor_directory as monitor


class StopAfterAcknowledgement(Exception):
    pass


class MotionMonitorShutdownTests(unittest.TestCase):
    def test_done_marker_is_deleted_after_final_teardown_and_monitor_freezes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            old_metadata = root / "protocol_metadata.json"
            old_metadata.write_text("{}")
            close_path = root / "protocol_CLOSE_20260925_120000.closeM"
            close_path.write_text("close")

            processor = MagicMock()
            processor.accumulator = SimpleNamespace(count=0)
            processor.pending_volumes = set()
            processor.processed_volumes = set()
            retirer = MagicMock()
            retirer.mode = "on"
            remove_observations = []
            real_unlink = Path.unlink

            def remove_after_teardown(path):
                remove_observations.append(
                    (processor.reset.called, retirer.reset.called)
                )
                real_unlink(Path(path))

            def stop_in_closed_state(_seconds):
                self.assertFalse(close_path.exists())
                self.assertTrue(processor.reset.called)
                self.assertTrue(retirer.reset.called)
                raise StopAfterAcknowledgement

            with (
                patch.object(monitor, "TSNRVolumeProcessor", return_value=processor),
                patch.object(
                    monitor,
                    "TransformRetirementCoordinator",
                    return_value=retirer,
                ),
                patch.object(monitor, "reset_logging", return_value=None),
                patch.object(monitor.os, "remove", side_effect=remove_after_teardown),
                patch.object(monitor.time, "sleep", side_effect=stop_in_closed_state),
            ):
                with self.assertRaises(StopAfterAcknowledgement):
                    monitor.monitor_directory(directory, 50.0, 0.3, 8080, "off")

            self.assertEqual(remove_observations, [(True, True)])
            self.assertTrue(old_metadata.is_file())
            processor.retry_pending.assert_called_once()
            processor.reset.assert_called_once()
            retirer.retire_pending_released.assert_called_once()
            retirer.reset.assert_called_once()


if __name__ == "__main__":
    unittest.main()

"""Focused checks for immediate consumption of atomically published pointers."""

from __future__ import annotations

import builtins
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))

from apps.python_fire_server.pointer_file import write_pointer_file_atomically
from apps.queue_processor import process_queue_directory as queue


class StopAfterClose(Exception):
    pass


class AtomicPointerConsumptionTests(unittest.TestCase):
    def _run_until_close(self, root: Path) -> list[float]:
        close_path = root / "done.closeQ"
        close_path.write_text("")
        sleep_calls: list[float] = []

        def make_identity(_input_dir, _reference, index, volume, group):
            identity = root / (
                f"alignTransform_{index:04d}_{volume:04d}-{group:04d}_identity.tfm"
            )
            identity.write_text("identity")
            return str(identity)

        def fake_sleep(seconds):
            sleep_calls.append(seconds)
            if close_path.exists():
                self.fail("queue slept before processing the acquisition close")
            raise StopAfterClose

        with patch.object(queue, "resample_nrrd_volume", side_effect=lambda path, **_: path), \
             patch.object(queue, "create_identityTransformFile", side_effect=make_identity), \
             patch.object(queue, "read_slice_timings_from_json", return_value=[0.0]), \
             patch.object(queue, "reset_logging", return_value=None), \
             patch.object(queue.time, "sleep", side_effect=fake_sleep):
            with self.assertRaises(StopAfterClose):
                queue.monitor_directory(root, "on", "cuda")

        return sleep_calls

    def test_final_pointer_is_read_without_wait_and_temporary_name_is_ignored(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            (root / "reference.nhdr").write_text("reference")
            pointer = root / "protocol_volume_0000_group_0000.txt"
            write_pointer_file_atomically(pointer, ["reference.nhdr"])
            temporary_pointer = root / f"{pointer.name}.deadbeef.tmp"
            temporary_pointer.write_text("not-yet-published.nhdr\n")

            sleep_calls = self._run_until_close(root)

            status = json.loads((root / "registration_status.json").read_text())
            self.assertEqual(
                status["volumes"]["0"]["groups"]["0"]["status"],
                "registered",
            )
            self.assertTrue(temporary_pointer.is_file())
            self.assertEqual(sleep_calls, [0.005])

    def test_pointer_disappearing_before_open_is_nonfatal(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            pointer = root / "protocol_volume_0000_group_0000.txt"
            write_pointer_file_atomically(pointer, ["reference.nhdr"])
            real_open = builtins.open

            def disappearing_open(filepath, *args, **kwargs):
                if os.path.abspath(os.fspath(filepath)) == str(pointer):
                    pointer.unlink()
                    raise FileNotFoundError("injected disappearance")
                return real_open(filepath, *args, **kwargs)

            with patch("builtins.open", side_effect=disappearing_open), \
                 self.assertLogs(level="WARNING") as captured:
                sleep_calls = self._run_until_close(root)

            self.assertEqual(sleep_calls, [0.005])
            self.assertFalse(pointer.exists())
            self.assertTrue(
                any(
                    "Pointer disappeared before it could be read" in line
                    for line in captured.output
                )
            )

    def test_empty_pointer_is_rejected_without_crashing_queue(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            pointer = root / "protocol_volume_0000_group_0000.txt"
            write_pointer_file_atomically(pointer, [])

            with self.assertLogs(level="ERROR") as captured:
                sleep_calls = self._run_until_close(root)

            self.assertEqual(sleep_calls, [0.005])
            self.assertTrue(pointer.is_file())
            self.assertTrue(
                any("pointer contains no target filenames" in line for line in captured.output)
            )

    def test_pointer_read_oserror_is_nonfatal(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            pointer = root / "protocol_volume_0000_group_0000.txt"
            write_pointer_file_atomically(pointer, ["reference.nhdr"])
            real_open = builtins.open

            def failing_open(filepath, *args, **kwargs):
                if os.path.abspath(os.fspath(filepath)) == str(pointer):
                    raise OSError("injected read failure")
                return real_open(filepath, *args, **kwargs)

            with patch("builtins.open", side_effect=failing_open), \
                 self.assertLogs(level="ERROR") as captured:
                sleep_calls = self._run_until_close(root)

            self.assertEqual(sleep_calls, [0.005])
            self.assertTrue(pointer.is_file())
            self.assertTrue(
                any("injected read failure" in line for line in captured.output)
            )


if __name__ == "__main__":
    unittest.main()

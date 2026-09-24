"""Focused observational checks for FIFO-off queue profiling."""

import csv
import json
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from apps.queue_processor import process_queue_directory as queue


class StopAfterClose(Exception):
    pass


class QueueProfileTests(unittest.TestCase):
    def run_small_acquisition(self, profile_flag):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            (root / "ref.nhdr").write_text("reference")
            (root / "volume_0000_group_0002.txt").write_text("ref.nhdr\n")
            calls = []
            sleep_calls = []

            def make_reference(_input_dir, _reference, _index, _volume, _group):
                identity = root / "alignTransform_0000_0000-0002_identity.tfm"
                identity.write_text("identity")
                (root / "target0.nhdr").write_text("target 0")
                (root / "target1.nhdr").write_text("target 1")
                old = root / "volume_0001_group_0000.txt"
                newest = root / "volume_0001_group_0001.txt"
                old.write_text("target0.nhdr\n")
                newest.write_text("target1.nhdr\n")
                # Match FIRE's unpublished pointer convention. The final
                # extension is .tmp, so queue discovery must ignore it.
                (root / "volume_0001_group_0002.txt.deadbeef.tmp").write_text(
                    "target1.nhdr\n"
                )
                now = time.time()
                os.utime(old, (now - 2, now - 2))
                os.utime(newest, (now - 1, now - 1))
                return str(identity)

            def fake_registration(_reference, _targets, _initial, label,
                                  _engine, persistent_cuda):
                calls.append((label, persistent_cuda))
                (root / f"alignTransform_{label}.tfm").write_text("registered")
                (root / "done.closeQ").write_text("")
                return True

            def fake_sleep(_seconds):
                sleep_calls.append(_seconds)
                # The close trigger has been handled and the acquisition reset.
                if (root / "registration_status.json").exists() and not (root / "done.closeQ").exists():
                    document = json.loads((root / "registration_status.json").read_text())
                    if "1" in document["volumes"]:
                        raise StopAfterClose

            persistent = object() if profile_flag == "on" else None
            with patch.object(queue, "resample_nrrd_volume", side_effect=lambda path, **_: path), \
                 patch.object(queue, "create_identityTransformFile", side_effect=make_reference), \
                 patch.object(queue, "read_slice_timings_from_json", return_value=[0.0, 1.0]), \
                 patch.object(queue, "run_MIregistration", side_effect=fake_registration), \
                 patch.object(queue, "reset_logging", return_value=None), \
                 patch.object(queue.time, "sleep", side_effect=fake_sleep):
                with self.assertRaises(StopAfterClose):
                    queue.monitor_directory(
                        directory, "off", "cuda", persistent,
                        profile_flag=profile_flag,
                        cuda_execution_mode="persistent" if persistent else "standalone",
                    )

            self.assertEqual(len(calls), 1)
            # The only sleep is the idle wait after acquisition close. Final
            # pointers and LIFO-skipped pointers incur no stability sleep.
            self.assertEqual(sleep_calls, [0.005])
            self.assertEqual(calls[0][0], "0001_0001-0001")
            self.assertIs(calls[0][1], persistent)
            status = json.loads((root / "registration_status.json").read_text())
            self.assertEqual(status["volumes"]["1"]["groups"]["0"]["status"], "skipped")
            self.assertEqual(status["volumes"]["1"]["groups"]["1"]["status"], "registered")
            self.assertTrue((root / "volume_0001_group_0002.txt.deadbeef.tmp").is_file())
            item_files = list(root.glob("queue_profile_items_*.csv"))
            scan_files = list(root.glob("queue_profile_scans_*.csv"))
            if profile_flag == "off":
                self.assertEqual(item_files, [])
                self.assertEqual(scan_files, [])
                return
            self.assertEqual(len(item_files), 1)
            self.assertEqual(len(scan_files), 1)
            with item_files[0].open(newline="") as source:
                rows = list(csv.DictReader(source))
            with scan_files[0].open(newline="") as source:
                scans = list(csv.DictReader(source))
            self.assertEqual(len(rows), 3)
            by_group = {(int(row["volume"]), int(row["group"])): row for row in rows}
            skipped = by_group[(1, 0)]
            registered = by_group[(1, 1)]
            self.assertEqual(skipped["outcome"], "skipped_by_lifo")
            self.assertEqual(skipped["skip_reason"], "Old file found. Skipping")
            self.assertEqual(skipped["selected_perf_ns"], "")
            self.assertEqual(skipped["pointer_candidate_count"], "2")
            self.assertTrue(skipped["file_mtime_ns"])
            self.assertTrue(skipped["status_end_perf_ns"])
            self.assertEqual(registered["outcome"], "registered")
            self.assertEqual(registered["output_label"], "0001_0001-0001")
            self.assertLess(int(registered["selected_perf_ns"]),
                            int(registered["registration_start_perf_ns"]))
            self.assertLess(int(registered["registration_start_perf_ns"]),
                            int(registered["registration_end_perf_ns"]))
            self.assertLess(int(registered["registration_end_perf_ns"]),
                            int(registered["loop_ready_perf_ns"]))
            self.assertLess(int(registered["reference_update_start_perf_ns"]),
                            int(registered["reference_update_end_perf_ns"]))
            self.assertTrue(any(scan["pointer_candidate_count"] == "2" and
                                scan["lifo_old_pointer_count"] == "1" for scan in scans))

    def test_disabled_keeps_skip_and_registration_without_files(self):
        self.run_small_acquisition("off")

    def test_enabled_records_skip_registration_and_scan(self):
        self.run_small_acquisition("on")

    def test_invalid_flag_fails_before_startup(self):
        with tempfile.TemporaryDirectory() as directory, \
             patch.dict(os.environ, {"QUEUE_PROFILE_FLAG": "unexpected"}), \
             patch.object(sys, "argv", ["process_queue_directory.py", directory]), \
             self.assertRaisesRegex(SystemExit, "2"):
            queue.main()


if __name__ == "__main__":
    unittest.main()

"""Regression coverage for interrupted-acquisition queue finalization."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared"))

from apps.queue_processor import process_queue_directory as queue
from apps.queue_processor.registration_status import RegistrationStatusTracker


class StopAfterClose(Exception):
    pass


class InterruptedAcquisitionShutdownTests(unittest.TestCase):
    protocol = "ep2d_bold_fire_noMoCo_tr1210_REF"

    def _pointer(self, root: Path, volume: int, group: int, group_size: int = 2) -> None:
        path = root / (
            f"{self.protocol}_volume_{volume:04d}_group_{group:04d}.txt"
        )
        first_slice = group * group_size
        slice_names = [
            f"{self.protocol}_volume_{volume:04d}_slice_{index:04d}.nhdr"
            for index in range(first_slice, first_slice + group_size)
        ]
        path.write_text("".join(f"{name}\n" for name in slice_names))
        for name in slice_names:
            (root / name).write_text("synthetic slice")

    def _run_queue(self, root: Path, slice_count: int = 20):
        sleep_calls: list[float] = []
        close_removals: list[str] = []
        real_remove = os.remove

        def make_identity(_input_dir, _reference, index, volume, group):
            identity = root / (
                f"alignTransform_{index:04d}_{volume:04d}-{group:04d}_identity.tfm"
            )
            identity.write_text("identity")
            return str(identity)

        def register(_reference, _targets, _initial, label, _engine, _persistent):
            (root / f"alignTransform_{label}.tfm").write_text("registered")
            return True

        def observe_remove(path):
            if str(path).endswith(".closeQ"):
                close_removals.append(os.path.basename(path))
            real_remove(path)

        def stop_after_acknowledgement(seconds):
            sleep_calls.append(seconds)
            if not list(root.glob("*.closeQ")):
                raise StopAfterClose

        with (
            patch.object(queue, "resample_nrrd_volume", side_effect=lambda path, **_: path),
            patch.object(queue, "create_identityTransformFile", side_effect=make_identity),
            patch.object(
                queue,
                "read_slice_timings_from_json",
                return_value=list(range(slice_count)),
            ),
            patch.object(queue, "run_MIregistration", side_effect=register),
            patch.object(queue, "reset_logging", return_value=None),
            patch.object(queue.os, "remove", side_effect=observe_remove),
            patch.object(queue.time, "sleep", side_effect=stop_after_acknowledgement),
        ):
            with self.assertRaises(StopAfterClose):
                queue.monitor_directory(root, "off", "cuda")

        status = json.loads((root / queue.REGISTRATION_STATUS_FILENAME).read_text())
        return status, sleep_calls, close_removals

    def test_october_8_partial_tail_is_terminalized_and_acknowledged_once(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            self._pointer(root, 0, 0)
            for slice_index in range(2, 20):
                (root / f"{self.protocol}_volume_0000_slice_{slice_index:04d}.nhdr").write_text(
                    "synthetic reference slice"
                )
            for volume in range(1, 11):
                for group in range(10):
                    self._pointer(root, volume, group)
            for group in range(7):
                self._pointer(root, 11, group)
            close_path = root / f"{self.protocol}_CLOSE_20261008_172909.closeQ"
            close_path.write_text("close")

            status, sleep_calls, close_removals = self._run_queue(root)

            self.assertEqual(
                len(list(root.glob(f"{self.protocol}_volume_0011_slice_*.nhdr"))),
                14,
            )
            groups = status["volumes"]["11"]["groups"]
            for group in range(7):
                self.assertEqual(groups[str(group)]["status"], "registered")
                self.assertNotIn("reason", groups[str(group)])
            for group in range(7, 10):
                self.assertEqual(
                    groups[str(group)],
                    {
                        "status": "failed",
                        "reason": "aborted_missing_pointer",
                    },
                )
            self.assertEqual(status["revision"], 12)
            self.assertFalse(close_path.exists())
            self.assertEqual(close_removals, [close_path.name])
            self.assertEqual(sleep_calls, [0.005])

    def test_complete_tail_adds_no_aborted_groups(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            self._pointer(root, 0, 0)
            for group in range(10):
                self._pointer(root, 1, group)
            (root / "done.closeQ").write_text("close")

            status, _, close_removals = self._run_queue(root)

            groups = status["volumes"]["1"]["groups"]
            self.assertEqual(len(groups), 10)
            self.assertTrue(all(record["status"] == "registered" for record in groups.values()))
            self.assertTrue(all("reason" not in record for record in groups.values()))
            self.assertEqual(close_removals, ["done.closeQ"])

    def test_close_before_first_image_publishes_empty_status_and_acknowledges(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            close_path = root / "metadata_only.closeQ"
            close_path.write_text("close")

            status, _, close_removals = self._run_queue(root)

            self.assertEqual(status["revision"], 0)
            self.assertEqual(status["volumes"], {})
            self.assertEqual(close_removals, [close_path.name])

    def test_pending_close_retries_at_bounded_rate_without_republication_loop(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            close_path = root / "pending.closeQ"
            close_path.write_text("close")
            real_publish = RegistrationStatusTracker.publish
            publish_attempts = 0
            sleep_calls = 0
            monotonic_now = 0

            def flaky_publish(tracker):
                nonlocal publish_attempts
                publish_attempts += 1
                if publish_attempts == 1:
                    raise OSError("injected publication failure")
                return real_publish(tracker)

            def advance_clock():
                nonlocal monotonic_now
                monotonic_now += 10_000_000
                return monotonic_now

            def stop_after_acknowledgement(_seconds):
                nonlocal sleep_calls
                sleep_calls += 1
                if not close_path.exists():
                    raise StopAfterClose
                if sleep_calls > 100:
                    self.fail("pending close did not resolve")

            with (
                patch.object(queue, "reset_logging", return_value=None),
                patch.object(
                    RegistrationStatusTracker,
                    "publish",
                    autospec=True,
                    side_effect=flaky_publish,
                ),
                patch.object(queue.time, "monotonic_ns", side_effect=advance_clock),
                patch.object(queue.time, "sleep", side_effect=stop_after_acknowledgement),
                self.assertLogs(level="ERROR") as captured,
            ):
                with self.assertRaises(StopAfterClose):
                    queue.monitor_directory(root, "on", "cuda")

            self.assertEqual(publish_attempts, 2)
            self.assertGreaterEqual(sleep_calls, 25)
            self.assertLess(sleep_calls, 35)
            self.assertEqual(
                sum("injected publication failure" in line for line in captured.output),
                1,
            )
            status = json.loads((root / queue.REGISTRATION_STATUS_FILENAME).read_text())
            self.assertEqual(status, {"schema_version": 1, "revision": 0, "volumes": {}})

    def test_marker_delete_retry_does_not_repeat_teardown_or_status_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            close_path = root / "delete_retry.closeQ"
            close_path.write_text("close")
            real_remove = os.remove
            real_replace = os.replace
            delete_attempts = 0
            status_replacements = 0
            monotonic_now = 0

            def flaky_remove(path):
                nonlocal delete_attempts
                if str(path).endswith(".closeQ"):
                    delete_attempts += 1
                    if delete_attempts == 1:
                        raise PermissionError("injected marker delete failure")
                real_remove(path)

            def observe_replace(source, destination):
                nonlocal status_replacements
                if os.path.basename(destination) == queue.REGISTRATION_STATUS_FILENAME:
                    status_replacements += 1
                real_replace(source, destination)

            def advance_clock():
                nonlocal monotonic_now
                monotonic_now += queue.CLOSE_RETRY_INTERVAL_NS
                return monotonic_now

            def stop_after_acknowledgement(_seconds):
                if not close_path.exists():
                    raise StopAfterClose

            with (
                patch.object(queue, "reset_logging", return_value=None),
                patch.object(queue.os, "remove", side_effect=flaky_remove),
                patch.object(queue.os, "replace", side_effect=observe_replace),
                patch.object(queue.time, "monotonic_ns", side_effect=advance_clock),
                patch.object(queue.time, "sleep", side_effect=stop_after_acknowledgement),
                self.assertLogs(level="INFO") as captured,
            ):
                with self.assertRaises(StopAfterClose):
                    queue.monitor_directory(root, "off", "cuda")

            self.assertEqual(delete_attempts, 2)
            self.assertEqual(status_replacements, 1)
            self.assertEqual(
                sum("Queue DONE" in line for line in captured.output),
                1,
            )

    def test_lifo_never_discards_close_control_file_on_mtime_tie(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")
            self._pointer(root, 0, 0, group_size=1)
            close_path = root / "same_time.closeQ"

            def make_identity(_input_dir, _reference, index, volume, group):
                identity = root / (
                    f"alignTransform_{index:04d}_{volume:04d}-{group:04d}_identity.tfm"
                )
                identity.write_text("identity")
                self._pointer(root, 1, 0, group_size=1)
                close_path.write_text("close")
                shared_mtime_ns = 1_800_000_000_000_000_000
                os.utime(
                    root / f"{self.protocol}_volume_0001_group_0000.txt",
                    ns=(shared_mtime_ns, shared_mtime_ns),
                )
                os.utime(close_path, ns=(shared_mtime_ns, shared_mtime_ns))
                return str(identity)

            def register(_reference, _targets, _initial, label, _engine, _persistent):
                (root / f"alignTransform_{label}.tfm").write_text("registered")
                return True

            def stop_after_acknowledgement(_seconds):
                if not close_path.exists():
                    raise StopAfterClose

            with (
                patch.object(queue, "resample_nrrd_volume", side_effect=lambda path, **_: path),
                patch.object(queue, "create_identityTransformFile", side_effect=make_identity),
                patch.object(queue, "read_slice_timings_from_json", return_value=[0]),
                patch.object(queue, "run_MIregistration", side_effect=register),
                patch.object(queue, "reset_logging", return_value=None),
                patch.object(queue.time, "sleep", side_effect=stop_after_acknowledgement),
            ):
                with self.assertRaises(StopAfterClose):
                    queue.monitor_directory(root, "off", "cuda")

            status = json.loads((root / queue.REGISTRATION_STATUS_FILENAME).read_text())
            self.assertEqual(
                status["volumes"]["1"]["groups"]["0"]["status"],
                "registered",
            )


if __name__ == "__main__":
    unittest.main()

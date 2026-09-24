"""Focused tests for monotonic FIRE MOCO transform advancement."""

import csv
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import handle_data
import moco_state
from handle_data import handleData
from moco_profile import FireMocoProfiler
from moco_state import FIRE_MOCO_STATE_FILENAME, FireMocoStatePublisher
from moco_transform_discovery import (
    discover_newest_registration_transform,
    parse_registration_transform_filename,
)


def transform_name(index, volume=1, group=0, suffix=""):
    return f"alignTransform_{index:04d}_{volume:04d}-{group:04d}{suffix}.tfm"


class RecordingProfiler:
    def __init__(self):
        self.rows = []

    def record(self, row):
        self.rows.append(dict(row))


class FakeConnection:
    def __init__(self, failures=0):
        self.failures = failures
        self.sent = []

    def send_feedback(self, name, data):
        if self.failures:
            self.failures -= 1
            raise OSError("simulated send failure")
        self.sent.append((name, data))


def build_handler(directory, connection=None):
    instance = handleData.__new__(handleData)
    instance.datafolder = str(directory)
    instance.imageNo = 45
    instance.sliceNo = 1
    instance.volcount = 1
    instance.groupsize = 1
    instance.framenumber = 0
    instance.frameNumberLookup = []
    instance.last_consumed_registration_index = None
    instance.last_consumed_transform_filename = None
    instance.moco_profiler = RecordingProfiler()
    instance.moco_state_publisher = FireMocoStatePublisher(directory)
    instance.refImgCoordFrame = object()
    instance.connection = connection if connection is not None else FakeConnection()
    return instance


def write_transform(directory, index, volume=1, group=0):
    path = Path(directory, transform_name(index, volume, group))
    path.write_bytes(f"transform-{index}".encode())
    return path


class TransformFilenameTests(unittest.TestCase):
    def test_production_filename_parsing(self):
        parsed = parse_registration_transform_filename(
            "alignTransform_0105_0012-0043.tfm"
        )
        self.assertEqual(parsed.registration_index, 105)
        self.assertEqual(parsed.volume, 12)
        self.assertEqual(parsed.group, 43)

        # Formatting is a minimum width, so long acquisitions may exceed four digits.
        parsed = parse_registration_transform_filename(
            "alignTransform_10005_0234-0001.tfm"
        )
        self.assertEqual(parsed.registration_index, 10005)

    def test_identity_malformed_and_unrelated_files_are_ignored(self):
        invalid = (
            "alignTransform_0000_0000-0000_identity.tfm",
            "alignTransform_0001.tfm",
            "alignTransform_x001_0001-0000.tfm",
            "foo.tfm",
            "alignTransform_0001_0001-0000.tfm.tmp",
        )
        for filename in invalid:
            with self.subTest(filename=filename):
                self.assertIsNone(parse_registration_transform_filename(filename))


class TransformDiscoveryTests(unittest.TestCase):
    def test_single_newer_transform_is_selected(self):
        with tempfile.TemporaryDirectory() as directory:
            write_transform(directory, 1)
            discovery = discover_newest_registration_transform(directory, None)
            self.assertEqual(discovery.root_entry_count, 1)
            self.assertEqual(discovery.valid_transform_count, 1)
            self.assertEqual(discovery.newer_candidate_count, 1)
            self.assertEqual(discovery.selected.registration_index, 1)

    def test_highest_index_wins_across_gaps_without_mtime_ordering(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = [write_transform(directory, index) for index in (101, 105, 109)]
            os.utime(paths[0], (300, 300))
            os.utime(paths[1], (200, 200))
            os.utime(paths[2], (100, 100))

            with patch("os.path.getmtime") as getmtime:
                discovery = discover_newest_registration_transform(directory, 100)

            getmtime.assert_not_called()
            self.assertEqual(discovery.newer_candidate_count, 3)
            self.assertEqual(discovery.selected.registration_index, 109)

    def test_consumed_history_cannot_resurface(self):
        with tempfile.TemporaryDirectory() as directory:
            for index in range(100, 106):
                write_transform(directory, index)
            discovery = discover_newest_registration_transform(directory, 105)
            self.assertEqual(discovery.valid_transform_count, 6)
            self.assertEqual(discovery.newer_candidate_count, 0)
            self.assertIsNone(discovery.selected)

            write_transform(directory, 106)
            discovery = discover_newest_registration_transform(directory, 105)
            self.assertEqual(discovery.selected.registration_index, 106)


class MocoAdvancementTests(unittest.TestCase):
    def run_feedback(self, instance, *, read=None, convert=None, package=None):
        read = object() if read is None else read
        convert = object() if convert is None else convert
        package = object() if package is None else package
        with (
            patch.object(handle_data.sitk, "ReadTransform", return_value=read),
            patch.object(handle_data, "convert_transform_for_moco", return_value=convert),
            patch.object(handle_data, "package_transform_as_datastruct", return_value=package),
            patch.object(handle_data, "format_moco_struct", return_value="feedback"),
        ):
            instance.package_and_send_moco_feedback()
        return read, convert, package

    def test_success_commits_state_and_frame_mapping_after_send(self):
        with tempfile.TemporaryDirectory() as directory:
            path = write_transform(directory, 5, volume=2, group=3)
            original = path.read_bytes()
            instance = build_handler(directory)

            registration, converted, packaged = self.run_feedback(instance)

            self.assertEqual(instance.last_consumed_registration_index, 5)
            self.assertEqual(instance.last_consumed_transform_filename, path.name)
            self.assertEqual(instance.framenumber, 1)
            self.assertEqual(len(instance.frameNumberLookup), 1)
            entry = instance.frameNumberLookup[0]
            self.assertEqual(entry["frameNumber"], 1)
            self.assertIs(entry["regTransformObject"], registration)
            self.assertIs(entry["mocoTransformObject"], converted)
            self.assertIs(entry["mocoStruct"], packaged)
            self.assertEqual(instance.moco_profiler.rows[-1]["outcome"], "feedback_sent")
            self.assertEqual(path.read_bytes(), original)
            self.assertTrue(path.exists())
            with open(Path(directory, FIRE_MOCO_STATE_FILENAME)) as stream:
                published = json.load(stream)
            self.assertEqual(published["committed_registration_index"], 5)
            self.assertEqual(published["release_before_registration_index"], 5)
            self.assertEqual(published["committed_transform_filename"], path.name)

    def test_no_new_transform_does_not_reopen_consumed_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = write_transform(directory, 5)
            instance = build_handler(directory)
            instance.last_consumed_registration_index = 5
            instance.last_consumed_transform_filename = path.name

            with patch.object(handle_data.sitk, "ReadTransform") as read:
                instance.package_and_send_moco_feedback()

            read.assert_not_called()
            self.assertEqual(instance.connection.sent, [])
            self.assertEqual(instance.moco_profiler.rows[-1]["outcome"], "no_new_transform")

    def test_read_conversion_and_packaging_failures_do_not_advance(self):
        cases = ("read", "conversion", "packaging")
        for failure in cases:
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as directory:
                write_transform(directory, 5)
                instance = build_handler(directory)
                read_effect = OSError("read failed") if failure == "read" else object()
                convert_effect = RuntimeError("conversion failed") if failure == "conversion" else object()
                package_effect = RuntimeError("packaging failed") if failure == "packaging" else object()

                with (
                    patch.object(
                        handle_data.sitk,
                        "ReadTransform",
                        side_effect=read_effect if isinstance(read_effect, Exception) else None,
                        return_value=None if isinstance(read_effect, Exception) else read_effect,
                    ),
                    patch.object(
                        handle_data,
                        "convert_transform_for_moco",
                        side_effect=convert_effect if isinstance(convert_effect, Exception) else None,
                        return_value=None if isinstance(convert_effect, Exception) else convert_effect,
                    ),
                    patch.object(
                        handle_data,
                        "package_transform_as_datastruct",
                        side_effect=package_effect if isinstance(package_effect, Exception) else None,
                        return_value=None if isinstance(package_effect, Exception) else package_effect,
                    ),
                    patch.object(handle_data, "format_moco_struct", return_value="feedback"),
                ):
                    instance.package_and_send_moco_feedback()

                self.assertIsNone(instance.last_consumed_registration_index)
                self.assertEqual(instance.frameNumberLookup, [])
                self.assertEqual(instance.framenumber, 0)
                self.assertEqual(instance.connection.sent, [])
                self.assertEqual(
                    instance.moco_profiler.rows[-1]["outcome"],
                    f"{failure if failure != 'read' else 'transform_read'}_failed",
                )
                self.assertFalse(Path(directory, FIRE_MOCO_STATE_FILENAME).exists())

    def test_send_failure_is_retryable_and_does_not_commit_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            write_transform(directory, 105)
            connection = FakeConnection(failures=1)
            instance = build_handler(directory, connection)

            self.run_feedback(instance)
            self.assertIsNone(instance.last_consumed_registration_index)
            self.assertEqual(instance.frameNumberLookup, [])
            self.assertEqual(instance.framenumber, 0)
            self.assertEqual(instance.moco_profiler.rows[-1]["outcome"], "send_failed")
            self.assertFalse(Path(directory, FIRE_MOCO_STATE_FILENAME).exists())

            self.run_feedback(instance)
            self.assertEqual(instance.last_consumed_registration_index, 105)
            self.assertEqual(instance.framenumber, 1)
            self.assertEqual(len(instance.frameNumberLookup), 1)

    def test_read_failure_leaves_same_candidate_eligible_for_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            write_transform(directory, 105)
            instance = build_handler(directory)

            with patch.object(
                handle_data.sitk,
                "ReadTransform",
                side_effect=OSError("read failed"),
            ):
                instance.package_and_send_moco_feedback()
            self.assertIsNone(instance.last_consumed_registration_index)

            self.run_feedback(instance)
            self.assertEqual(instance.last_consumed_registration_index, 105)

    def test_newer_candidate_supersedes_failed_candidate(self):
        with tempfile.TemporaryDirectory() as directory:
            write_transform(directory, 105)
            instance = build_handler(directory, FakeConnection(failures=1))
            self.run_feedback(instance)
            write_transform(directory, 106)

            self.run_feedback(instance)
            self.assertEqual(instance.last_consumed_registration_index, 106)
            self.assertEqual(
                instance.frameNumberLookup[-1]["regTransformFilename"],
                transform_name(106),
            )
            published = json.loads(
                Path(directory, FIRE_MOCO_STATE_FILENAME).read_text()
            )
            self.assertEqual(published["committed_registration_index"], 106)

    def test_publication_failure_does_not_roll_back_commit_and_later_catches_up(self):
        with tempfile.TemporaryDirectory() as directory:
            write_transform(directory, 105)
            instance = build_handler(directory)

            with patch.object(
                moco_state.os,
                "replace",
                side_effect=OSError("state publication failed"),
            ):
                self.run_feedback(instance)

            self.assertEqual(instance.last_consumed_registration_index, 105)
            self.assertEqual(len(instance.connection.sent), 1)
            self.assertFalse(Path(directory, FIRE_MOCO_STATE_FILENAME).exists())

            write_transform(directory, 109)
            self.run_feedback(instance)
            published = json.loads(
                Path(directory, FIRE_MOCO_STATE_FILENAME).read_text()
            )
            self.assertEqual(instance.last_consumed_registration_index, 109)
            self.assertEqual(published["committed_registration_index"], 109)
            self.assertEqual(published["release_before_registration_index"], 109)

    def test_feedback_log_failure_prevents_advancement(self):
        with tempfile.TemporaryDirectory() as directory:
            write_transform(directory, 5)
            instance = build_handler(directory)
            real_open = open

            def fail_feedback_log(path, *args, **kwargs):
                if str(path).endswith("log_moco_feedback_sent.log"):
                    raise OSError("log failed")
                return real_open(path, *args, **kwargs)

            with (
                patch.object(handle_data.sitk, "ReadTransform", return_value=object()),
                patch.object(handle_data, "convert_transform_for_moco", return_value=object()),
                patch.object(handle_data, "package_transform_as_datastruct", return_value=object()),
                patch.object(handle_data, "format_moco_struct", return_value="feedback"),
                patch("builtins.open", side_effect=fail_feedback_log),
            ):
                instance.package_and_send_moco_feedback()

            self.assertIsNone(instance.last_consumed_registration_index)
            self.assertEqual(instance.connection.sent, [])
            self.assertEqual(
                instance.moco_profiler.rows[-1]["outcome"],
                "feedback_logging_failed",
            )
            self.assertFalse(Path(directory, FIRE_MOCO_STATE_FILENAME).exists())

    def test_transform_files_are_not_moved_deleted_or_rewritten(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = [write_transform(directory, index) for index in (1, 5, 9)]
            before = {path.name: path.read_bytes() for path in paths}
            instance = build_handler(directory)

            self.run_feedback(instance)

            after = {path.name: path.read_bytes() for path in paths}
            self.assertEqual(after, before)
            self.assertFalse(Path(directory, "processed_protocol").exists())

    def test_profiler_failure_does_not_change_successful_feedback(self):
        with tempfile.TemporaryDirectory() as directory:
            write_transform(directory, 5)
            instance = build_handler(directory)
            instance.moco_profiler = MagicMock()
            instance.moco_profiler.record.side_effect = OSError("profile failed")

            self.run_feedback(instance)

            self.assertEqual(instance.last_consumed_registration_index, 5)
            self.assertEqual(len(instance.connection.sent), 1)
            self.assertEqual(len(instance.frameNumberLookup), 1)

    def test_moco_disabled_constructor_does_not_create_profiler(self):
        connection = MagicMock()
        connection.is_exhausted = False
        with tempfile.TemporaryDirectory() as directory:
            with (
                patch.object(handle_data.h5py, "File", return_value=MagicMock()),
                patch.object(handle_data, "FireMocoProfiler") as profiler,
                patch.object(handle_data, "FireMocoStatePublisher") as publisher,
            ):
                instance = handleData(connection, directory, moco_enabled=False)

        profiler.assert_not_called()
        publisher.assert_not_called()
        self.assertIsNone(instance.moco_profiler)
        self.assertIsNone(instance.last_consumed_registration_index)


class MocoProfilerTests(unittest.TestCase):
    def test_csv_contains_timing_and_outcome_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            profiler = FireMocoProfiler(directory)
            profiler.record({
                "volume": 3,
                "selected_registration_index": 9,
                "discovery_selection_ms": 1.25,
                "outcome": "feedback_sent",
            })
            path = profiler.path
            profiler.close()

            with open(path, newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["volume"], "3")
            self.assertEqual(rows[0]["selected_registration_index"], "9")
            self.assertEqual(rows[0]["outcome"], "feedback_sent")

    def test_writer_failure_is_nonfatal_and_disables_profiler(self):
        with tempfile.TemporaryDirectory() as directory:
            profiler = FireMocoProfiler(directory)
            profiler._writer = MagicMock()
            profiler._writer.writerow.side_effect = OSError("disk full")

            profiler.record({"outcome": "no_new_transform"})

            self.assertTrue(profiler.disabled)
            profiler.close()


if __name__ == "__main__":
    unittest.main()

"""Focused tests for downstream barriers and complete acquisition archival."""

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import handle_data
from handle_data import AcquisitionFinalizationResult, handleData


class AcquisitionFinalizationTests(unittest.TestCase):
    def _handler(self, root: Path, protocol: str = "protocol_A") -> handleData:
        instance = handleData.__new__(handleData)
        instance.datafolder = str(root)
        instance.protocol_name = protocol
        return instance

    def test_close_request_is_atomically_published(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            handler = self._handler(root)
            real_replace = handle_data.os.replace
            with patch.object(
                handle_data.os,
                "replace",
                wraps=real_replace,
            ) as replace:
                close_path = Path(handler.generate_close_file("closeQ"))

            self.assertTrue(close_path.is_file())
            self.assertEqual(close_path.suffix, ".closeQ")
            replace.assert_called_once()
            temporary, final = replace.call_args.args
            self.assertTrue(temporary.endswith(".tmp"))
            self.assertEqual(final, str(close_path))
            self.assertFalse(Path(temporary).exists())

    def test_complete_archive_preserves_nested_processed_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            handler = self._handler(root)
            root_files = {
                "metadata.json": b"metadata",
                "alignTransform_0105_0002-0003.tfm": b"active transform",
                "fire_moco_state.moco": b"state",
                "measurement-20260925T120000.hdf5": b"closed hdf5",
            }
            for name, contents in root_files.items():
                (root / name).write_bytes(contents)

            processed = root / "processed_protocol_A"
            transforms = processed / "transforms"
            transforms.mkdir(parents=True)
            nested_files = {
                "volume_0001_group_0001.txt": b"pointer",
                "slice.nhdr": b"header",
                "slice.raw": b"raw",
                "transforms/alignTransform_0001_0001-0000.tfm": b"transform 1",
                "transforms/alignTransform_0002_0001-0001.tfm": b"transform 2",
            }
            for relative, contents in nested_files.items():
                (processed / relative).write_bytes(contents)

            stale = root / "processed_protocol_B"
            stale.mkdir()
            (stale / "unrelated.raw").write_bytes(b"unrelated")

            result = handler.consolidate_outputs_in_directory("savedData_test_protocol_A")

            destination = root / "savedData_test_protocol_A"
            self.assertTrue(result.success)
            self.assertEqual(result.root_files_moved, len(root_files))
            self.assertTrue(result.processed_directory_moved)
            self.assertFalse(processed.exists())
            self.assertTrue(stale.is_dir())
            for name, contents in root_files.items():
                self.assertEqual((destination / name).read_bytes(), contents)
            for relative, contents in nested_files.items():
                self.assertEqual(
                    (destination / "processed_protocol_A" / relative).read_bytes(),
                    contents,
                )
            self.assertFalse((destination / "volume_0001_group_0001.txt").exists())

    def test_missing_processed_directory_is_successful(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.json").write_text("{}")

            result = self._handler(root).consolidate_outputs_in_directory(
                "savedData_without_processed"
            )

            self.assertTrue(result.success)
            self.assertEqual(result.processed_directory_status, "missing")
            self.assertTrue(
                (root / "savedData_without_processed" / "metadata.json").is_file()
            )

    def test_processed_destination_collision_preserves_both_trees(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            handler = self._handler(root)
            source = root / "processed_protocol_A"
            source.mkdir()
            (source / "source.raw").write_bytes(b"source")
            saved = root / "savedData_collision"
            destination = saved / "processed_protocol_A"
            destination.mkdir(parents=True)
            (destination / "existing.raw").write_bytes(b"destination")

            status, error = handler.move_processed_acquisition_directory(saved)

            self.assertEqual(status, "collision")
            self.assertIsNotNone(error)
            self.assertEqual((source / "source.raw").read_bytes(), b"source")
            self.assertEqual(
                (destination / "existing.raw").read_bytes(), b"destination"
            )

    def test_processed_move_failure_is_nonfatal_and_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            handler = self._handler(root)
            source = root / "processed_protocol_A"
            source.mkdir()
            (source / "source.raw").write_bytes(b"source")
            saved = root / "savedData_failure"
            saved.mkdir()

            with patch.object(
                handle_data.os,
                "rename",
                side_effect=OSError("injected directory move failure"),
            ):
                status, error = handler.move_processed_acquisition_directory(saved)

            self.assertEqual(status, "failed")
            self.assertIn("injected directory move failure", error)
            self.assertTrue((source / "source.raw").is_file())
            self.assertFalse((saved / "processed_protocol_A").exists())

    def test_existing_saved_data_destination_is_not_merged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "metadata.json"
            source.write_text("source")
            destination = root / "savedData_existing"
            destination.mkdir()
            (destination / "metadata.json").write_text("destination")

            result = self._handler(root).consolidate_outputs_in_directory(
                destination.name
            )

            self.assertFalse(result.success)
            self.assertEqual(source.read_text(), "source")
            self.assertEqual((destination / "metadata.json").read_text(), "destination")

    def test_root_move_failure_preserves_source_and_other_moves_continue(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            handler = self._handler(root)
            good = root / "good.json"
            failed = root / "failed.tfm"
            good.write_text("good")
            failed.write_text("failed")
            processed = root / "processed_protocol_A"
            processed.mkdir()
            (processed / "retired.raw").write_text("retired")
            real_rename = handle_data.os.rename

            def selective_failure(source, destination):
                if Path(source).name == failed.name:
                    raise OSError("injected root move failure")
                return real_rename(source, destination)

            with patch.object(
                handle_data.os,
                "rename",
                side_effect=selective_failure,
            ):
                result = handler.consolidate_outputs_in_directory(
                    "savedData_partial"
                )

            destination = root / "savedData_partial"
            self.assertFalse(result.success)
            self.assertTrue(failed.is_file())
            self.assertTrue((destination / good.name).is_file())
            self.assertTrue(
                (destination / "processed_protocol_A" / "retired.raw").is_file()
            )
            self.assertIn("injected root move failure", result.errors[0])

    def test_finalize_waits_for_queue_then_monitor_before_consolidation(self):
        handler = self._handler(Path("/unused"))
        events = []

        def generate(extension):
            events.append(f"request_{extension}")
            return f"/{extension}"

        def wait(path):
            events.append(f"done_{Path(path).name}")

        expected = AcquisitionFinalizationResult(
            "savedData_test", "/unused/savedData_test", 2, "moved", ()
        )
        handler.generate_close_file = MagicMock(side_effect=generate)
        handler.wait_for_close_acknowledgement = MagicMock(side_effect=wait)
        handler.close_acquisition_file_for_finalization = MagicMock(
            side_effect=lambda: events.append("close_hdf5")
        )
        handler.consolidate_outputs_in_directory = MagicMock(
            side_effect=lambda: events.append("consolidate") or expected
        )

        result = handler.finalize_acquisition_outputs()

        self.assertIs(result, expected)
        self.assertEqual(
            events,
            [
                "request_closeQ",
                "done_closeQ",
                "request_closeM",
                "done_closeM",
                "close_hdf5",
                "consolidate",
            ],
        )

    def test_acquisition_hdf5_is_flushed_and_closed_before_movement(self):
        handler = self._handler(Path("/unused"))
        handler.hf = MagicMock()

        handler.close_acquisition_file_for_finalization()

        handler.hf.flush.assert_called_once_with()
        handler.hf.close.assert_called_once_with()

    def test_queue_or_monitor_timeout_never_starts_consolidation(self):
        for failing_wait in (1, 2):
            with self.subTest(failing_wait=failing_wait):
                handler = self._handler(Path("/unused"))
                handler.generate_close_file = MagicMock(
                    side_effect=lambda extension: f"/{extension}"
                )
                waits = 0

                def wait(_path):
                    nonlocal waits
                    waits += 1
                    if waits == failing_wait:
                        raise TimeoutError("downstream not done")

                handler.wait_for_close_acknowledgement = MagicMock(side_effect=wait)
                handler.close_acquisition_file_for_finalization = MagicMock()
                handler.consolidate_outputs_in_directory = MagicMock()

                with self.assertRaisesRegex(TimeoutError, "downstream not done"):
                    handler.finalize_acquisition_outputs()

                handler.consolidate_outputs_in_directory.assert_not_called()
                self.assertEqual(
                    handler.generate_close_file.call_count,
                    1 if failing_wait == 1 else 2,
                )


if __name__ == "__main__":
    unittest.main()

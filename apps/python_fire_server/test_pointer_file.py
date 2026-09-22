"""Focused tests for atomic queue-pointer publication."""

import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from apps.python_fire_server.pointer_file import write_pointer_file_atomically


class AtomicPointerFileTests(unittest.TestCase):
    def test_sequential_writes_are_complete_before_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            real_replace = os.replace
            publications = []
            expected_publications = iter([
                "slice_0000.nhdr\nslice_0001.nhdr\n",
                "slice_0002.nhdr\n",
            ])

            def observe_replace(source, destination):
                expected = next(expected_publications)
                self.assertEqual(os.path.dirname(source), directory)
                self.assertEqual(os.path.splitext(source)[1], ".tmp")
                self.assertFalse(os.path.exists(destination))
                self.assertEqual(Path(source).read_text(), expected)
                real_replace(source, destination)
                publications.append(expected)

            pointer_specs = [
                ("protocol_volume_0001_group_0000.txt", ["slice_0000.nhdr", "slice_0001.nhdr"]),
                ("protocol_volume_0001_group_0001.txt", ["slice_0002.nhdr"]),
            ]
            expected_contents = [
                "slice_0000.nhdr\nslice_0001.nhdr\n",
                "slice_0002.nhdr\n",
            ]

            with patch(
                "apps.python_fire_server.pointer_file.os.replace",
                side_effect=observe_replace,
            ):
                for name, entries in pointer_specs:
                    write_pointer_file_atomically(os.path.join(directory, name), entries)

            self.assertEqual(len(publications), 2)
            for (name, _), expected in zip(pointer_specs, expected_contents):
                self.assertEqual(Path(directory, name).read_text(), expected)
            self.assertEqual(
                sorted(os.listdir(directory)),
                sorted(name for name, _ in pointer_specs),
            )

    def test_failed_write_removes_temporary_file_without_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            pointer_path = os.path.join(directory, "protocol_volume_0001_group_0000.txt")

            def failing_entries():
                yield "slice_0000.nhdr"
                raise OSError("write failed")

            with patch("apps.python_fire_server.pointer_file.os.replace") as replace:
                with self.assertRaisesRegex(OSError, "write failed"):
                    write_pointer_file_atomically(pointer_path, failing_entries())

            replace.assert_not_called()
            self.assertFalse(os.path.exists(pointer_path))
            self.assertEqual(os.listdir(directory), [])

    def test_failed_publication_removes_temporary_file(self):
        with tempfile.TemporaryDirectory() as directory:
            pointer_path = os.path.join(directory, "protocol_volume_0001_group_0000.txt")

            with patch(
                "apps.python_fire_server.pointer_file.os.replace",
                side_effect=OSError("publish failed"),
            ):
                with self.assertRaisesRegex(OSError, "publish failed"):
                    write_pointer_file_atomically(pointer_path, ["slice_0000.nhdr"])

            self.assertFalse(os.path.exists(pointer_path))
            self.assertEqual(os.listdir(directory), [])


if __name__ == "__main__":
    unittest.main()

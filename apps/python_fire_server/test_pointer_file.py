"""Focused tests for atomic queue-pointer publication."""

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
from handle_data import handleData
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

    def test_slice_payload_is_ready_before_pointer_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            handler = handleData.__new__(handleData)
            handler.connection = [
                (1022, MagicMock(), b"header", b"attributes", b"pixels")
            ]
            image = handler.connection[0][1]
            image.getHead.return_value.user_float = [0.0] * 8
            image.position = [0.0, 0.0, 0.0]
            image.acquisition_time_stamp = 1
            handler.datafolder = str(root)
            handler.protocol_name = "protocol"
            handler.imageNo = 0
            handler.sliceNo = 0
            handler.volcount = 1
            handler.groupcount = 0
            handler.groupsize = 1
            handler.nslices_per_volume = 2
            handler.moco_enabled = False
            events = []

            def save_raw(data):
                self.assertEqual(data, b"pixels")
                raw_name = "protocol_volume_0001_slice_0000.raw"
                (root / raw_name).write_bytes(data)
                events.append("raw_closed")
                return raw_name

            def save_header(_image, raw_name):
                self.assertEqual(
                    (root / raw_name).read_bytes(),
                    b"pixels",
                )
                (root / "protocol_volume_0001_slice_0000.nhdr").write_text(
                    f"data file: {raw_name}\n"
                )
                events.append("header_closed")

            real_publish = handle_data.write_pointer_file_atomically

            def publish(pointer_path, header_names):
                self.assertEqual(events, ["raw_closed", "header_closed"])
                self.assertEqual(
                    (root / header_names[0]).read_text(),
                    "data file: protocol_volume_0001_slice_0000.raw\n",
                )
                events.append("pointer_published")
                real_publish(pointer_path, header_names)

            with patch.object(handler, "save_image", return_value=None), \
                 patch.object(
                     handler,
                     "save_raw_image_from_ismrmrd_data",
                     side_effect=save_raw,
                 ), \
                 patch.object(
                     handler,
                     "save_nrrd_slice_from_ismrmrd_data",
                     side_effect=save_header,
                 ), \
                 patch.object(
                     handle_data,
                     "write_pointer_file_atomically",
                     side_effect=publish,
                 ):
                handler.next()

            pointer = root / "protocol_volume_0001_group_0000.txt"
            self.assertEqual(events, ["raw_closed", "header_closed", "pointer_published"])
            self.assertEqual(
                pointer.read_text(),
                "protocol_volume_0001_slice_0000.nhdr\n",
            )

    def test_reference_header_is_ready_before_pointer_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw_names = [
                "protocol_volume_0000_slice_0000.raw",
                "protocol_volume_0000_slice_0001.raw",
            ]
            for index, raw_name in enumerate(raw_names):
                (root / raw_name).write_bytes(f"slice {index}".encode())

            handler = handleData.__new__(handleData)
            handler.datafolder = str(root)
            handler.protocol_name = "protocol"
            handler.volcount = 0
            handler.sliceNo = 2
            handler.groupsize = 2
            items = [
                {"raw_filename": raw_names[0], "image": MagicMock()},
                {"raw_filename": raw_names[1], "image": MagicMock()},
            ]
            events = []

            def generate_header(_image, header_path, listed_raws, _next_image):
                self.assertEqual(listed_raws, raw_names)
                self.assertEqual(
                    [(root / name).read_bytes() for name in listed_raws],
                    [b"slice 0", b"slice 1"],
                )
                Path(header_path).write_text("\n".join(listed_raws) + "\n")
                events.append("reference_header_closed")

            real_publish = handle_data.write_pointer_file_atomically

            def publish(pointer_path, header_names):
                self.assertEqual(events, ["reference_header_closed"])
                self.assertEqual(len(header_names), 1)
                self.assertEqual(
                    (root / header_names[0]).read_text(),
                    "\n".join(raw_names) + "\n",
                )
                events.append("pointer_published")
                real_publish(pointer_path, header_names)

            with patch.object(
                     handler,
                     "generate_nrrd_header",
                     side_effect=generate_header,
                 ), \
                 patch.object(
                     handle_data,
                     "write_pointer_file_atomically",
                     side_effect=publish,
                 ):
                handler.save_nrrd_volume_from_ismrmrd_data(items)

            self.assertEqual(events, ["reference_header_closed", "pointer_published"])
            pointers = list(root.glob("protocol_volume_0000_group_*.txt"))
            self.assertEqual(len(pointers), 1)
            published_header = pointers[0].read_text().strip()
            self.assertTrue((root / published_header).is_file())

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

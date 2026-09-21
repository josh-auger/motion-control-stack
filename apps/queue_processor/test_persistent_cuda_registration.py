"""Focused checks for the SLIMM persistent process protocol."""

import os
import stat
import tempfile
import unittest

try:
    from .persistent_cuda_registration import PersistentCudaRegistrationProcess
except ImportError:
    from persistent_cuda_registration import PersistentCudaRegistrationProcess


FAKE_SLIMM = '''#!/usr/bin/env python3
import sys

print("ordinary startup diagnostic", flush=True)
print("SLIMM_PROTOCOL\\tREADY\\tversion=1", flush=True)
for line in sys.stdin:
    fields = line.strip().split()
    if not fields:
        continue
    command = fields[0]
    if command == "EXIT":
        print("SLIMM_PROTOCOL\\tSUCCESS\\tcommand=EXIT", flush=True)
        break
    if command == "LOAD_REFERENCE":
        if "bad-reference" in line:
            print("SLIMM_PROTOCOL\\tERROR\\tcommand=LOAD_REFERENCE", flush=True)
        else:
            print("SLIMM_PROTOCOL\\tSUCCESS\\tcommand=LOAD_REFERENCE", flush=True)
    elif command == "REGISTER":
        print("ordinary registration diagnostic", flush=True)
        if "bad-register" in line:
            print("SLIMM_PROTOCOL\\tERROR\\tcommand=REGISTER", flush=True)
        elif "wrong-command" in line:
            print("SLIMM_PROTOCOL\\tSUCCESS\\tcommand=LOAD_REFERENCE", flush=True)
        elif "unexpected-eof" in line:
            sys.exit(2)
        else:
            print("SLIMM_PROTOCOL\\tSUCCESS\\tcommand=REGISTER", flush=True)
'''


class PersistentCudaRegistrationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.binary = os.path.join(self.directory.name, "fake-slimm")
        with open(self.binary, "w", encoding="utf-8") as output:
            output.write(FAKE_SLIMM)
        os.chmod(self.binary, stat.S_IRWXU)

    def tearDown(self):
        self.directory.cleanup()

    def test_reuses_process_and_reference_across_requests(self):
        manager = PersistentCudaRegistrationProcess(self.binary)
        process = manager.process
        try:
            pid = manager.process.pid
            with self.assertLogs(level="INFO") as captured:
                manager.load_reference("/tmp/reference.nhdr")
                manager.load_reference("/tmp/reference.nhdr")
                manager.register("first", ["0"] * 6, ["/tmp/one.nhdr"])
                manager.register("second", ["0"] * 6, ["/tmp/two.nhdr"])
            self.assertEqual(sum("LOAD_REFERENCE completed" in line for line in captured.output), 1)
            self.assertEqual(sum("ordinary registration diagnostic" in line for line in captured.output), 2)
            self.assertEqual(manager.process.pid, pid)
            self.assertIsNone(manager.process.poll())
        finally:
            manager.close()
        self.assertEqual(process.returncode, 0)

    def test_protocol_error_and_wrong_response_fail(self):
        manager = PersistentCudaRegistrationProcess(self.binary)
        try:
            with self.assertRaisesRegex(RuntimeError, "SLIMM protocol error"):
                manager.load_reference("bad-reference")
            self.assertIsNone(manager.reference_path)
            manager.load_reference("good-reference")
            with self.assertRaisesRegex(RuntimeError, "SLIMM protocol error"):
                manager.register("bad-register", ["0"] * 6, ["target"])
            with self.assertRaisesRegex(RuntimeError, "Unexpected SLIMM protocol response"):
                manager.register("wrong-command", ["0"] * 6, ["target"])
        finally:
            manager.close()

    def test_unexpected_eof_fails(self):
        manager = PersistentCudaRegistrationProcess(self.binary)
        try:
            manager.load_reference("reference")
            with self.assertRaisesRegex(RuntimeError, "exited before protocol response"):
                manager.register("unexpected-eof", ["0"] * 6, ["target"])
        finally:
            manager.close()

    def test_start_failure_fails(self):
        with self.assertRaises(FileNotFoundError):
            PersistentCudaRegistrationProcess("/nonexistent/cuda-registration")


if __name__ == "__main__":
    unittest.main()

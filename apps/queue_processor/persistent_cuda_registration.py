"""Synchronous protocol client for the long-lived SLIMM CUDA registration binary."""

import json
import logging
import os
import select
import subprocess
import time


CUDA_BINARY = "/opt/moco/bin/cuda-standalone-registration"


class PersistentCudaRegistrationProcess:
    def __init__(self, binary=CUDA_BINARY):
        self.process = None
        self.reference_path = None
        self._pending_output = b""
        start = time.monotonic()
        try:
            self.process = subprocess.Popen(
                [binary, "--persistent"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=0,
            )
            logging.info("Persistent CUDA registration process started (PID=%d)", self.process.pid)
            self._await_response("READY", timeout=30)
            logging.info("SLIMM persistent registration ready (%.3f s)", time.monotonic() - start)
        except (OSError, RuntimeError):
            self.close()
            raise

    def _read_line(self, deadline):
        """Read complete lines without losing bytes buffered after a protocol response."""
        while b"\n" not in self._pending_output:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise RuntimeError("Timed out waiting for SLIMM protocol response")
            ready, _, _ = select.select([self.process.stdout], [], [], remaining)
            if not ready:
                raise RuntimeError("Timed out waiting for SLIMM protocol response")
            try:
                chunk = os.read(self.process.stdout.fileno(), 4096)
            except OSError as error:
                raise RuntimeError(f"SLIMM stdout pipe failed: {error}") from error
            if not chunk:
                raise RuntimeError(
                    f"SLIMM process exited before protocol response (exit={self.process.poll()})"
                )
            self._pending_output += chunk
        line, self._pending_output = self._pending_output.split(b"\n", 1)
        return line.decode("utf-8", errors="replace").rstrip("\r")

    def _await_response(self, expected, timeout=120):
        deadline = time.monotonic() + timeout
        while True:
            line = self._read_line(deadline)
            if not line.startswith("SLIMM_PROTOCOL"):
                if line:
                    logging.info("SLIMM: %s", line)
                continue
            fields = line.split()
            if len(fields) < 2 or fields[0] != "SLIMM_PROTOCOL":
                raise RuntimeError(f"Malformed SLIMM protocol response: {line}")
            if fields[1] == "ERROR":
                raise RuntimeError(f"SLIMM protocol error: {line}")
            if expected == "READY" and fields[1] == "READY":
                return line
            if expected != "READY" and fields[1] == "SUCCESS":
                commands = [field[8:] for field in fields[2:] if field.startswith("command=")]
                if commands == [expected]:
                    return line
            raise RuntimeError(f"Unexpected SLIMM protocol response (expected {expected}): {line}")

    def _command(self, command, expected):
        if self.process.poll() is not None:
            raise RuntimeError(f"SLIMM process is not running (exit={self.process.returncode})")
        try:
            self.process.stdin.write((command + "\n").encode("utf-8"))
            self.process.stdin.flush()
        except (OSError, ValueError) as error:
            raise RuntimeError(f"SLIMM stdin pipe failed: {error}") from error
        return self._await_response(expected)

    def load_reference(self, path):
        if path == self.reference_path:
            return
        start = time.monotonic()
        response = self._command(f"LOAD_REFERENCE {json.dumps(path, ensure_ascii=False)}", "LOAD_REFERENCE")
        self.reference_path = path
        logging.info("SLIMM LOAD_REFERENCE completed in %.3f s: %s", time.monotonic() - start, response)

    def register(self, label, init_parameters, targets):
        if self.reference_path is None:
            raise RuntimeError("SLIMM reference has not been loaded")
        if len(init_parameters) != 6:
            raise ValueError("CUDA registration requires six initialization parameters")
        command = " ".join(
            ["REGISTER", label, *init_parameters, str(len(targets))]
            + [json.dumps(path, ensure_ascii=False) for path in targets]
        )
        start = time.monotonic()
        response = self._command(command, "REGISTER")
        logging.info("SLIMM REGISTER %s completed in %.3f s: %s", label, time.monotonic() - start, response)

    def close(self):
        process = self.process
        if process is None:
            return
        self.process = None
        if process.poll() is None:
            try:
                process.stdin.write(b"EXIT\n")
                process.stdin.flush()
                process.wait(timeout=5)
            except (OSError, ValueError, subprocess.TimeoutExpired):
                if process.poll() is None:
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
        if process.stdin:
            process.stdin.close()
        if process.stdout:
            process.stdout.close()
        logging.info("Persistent CUDA registration process stopped (PID=%d)", process.pid)

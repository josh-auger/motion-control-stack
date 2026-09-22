"""Opt-in, acquisition-local CSV profiling for the queue processor.

All *_perf_ns fields use time.perf_counter_ns() and may be subtracted from
each other within this process. *_wall_ns fields use time.time_ns();
file_mtime_ns is filesystem wall-clock metadata. Wall and file times may be
compared as an arrival proxy, subject to filesystem timestamp accuracy, but
neither may be subtracted from a perf-counter timestamp.

An item row describes one pointer. first_observed is the first successful
directory scan containing it; first_eligible is its first candidate snapshot.
selection_snapshot is when newest-file selection finished. selected is set
only for the file actually processed; LIFO discards have decision and skip
timestamps instead. loop_ready is after the whole candidate batch finishes,
immediately before the next directory scan. A scan row describes one nonempty
snapshot; pointer_candidate_count is that snapshot's count, not a durable
queue depth. Empty scans are aggregated into the next nonempty scan row.
"""

import csv
import os
import time
from datetime import datetime


ITEM_FIELDS = (
    "scan_id", "volume", "group", "pointer_filename", "pointer_path",
    "outcome", "skip_reason", "status_result", "output_label",
    "fifo_flag", "reg_engine", "cuda_execution_mode",
    "candidate_file_count", "pointer_candidate_count", "lifo_old_pointer_count",
    "file_mtime_ns", "first_observed_wall_ns", "first_observed_perf_ns",
    "first_eligible_perf_ns", "selection_snapshot_perf_ns", "selected_perf_ns",
    "decision_wall_ns", "decision_perf_ns", "preparation_start_perf_ns",
    "preparation_end_perf_ns", "reference_setup_start_perf_ns",
    "reference_setup_end_perf_ns", "registration_start_perf_ns",
    "registration_end_perf_ns", "status_start_perf_ns", "status_end_perf_ns",
    "reference_update_start_perf_ns", "reference_update_end_perf_ns",
    "item_complete_perf_ns", "loop_ready_perf_ns",
)

SCAN_FIELDS = (
    "scan_id", "scan_start_perf_ns", "scan_end_perf_ns", "loop_ready_perf_ns",
    "directory_entry_count", "candidate_file_count", "pointer_candidate_count", "selected_file_count",
    "lifo_old_pointer_count", "scan_and_build_ms", "newest_selection_ms",
    "skip_handling_ms", "total_cycle_ms", "empty_polls_since_previous",
    "empty_scan_ms_since_previous", "observed_pointer_count_total",
    "registered_total", "skipped_total",
)


class QueueProfiler:
    """Write two buffered CSV streams without participating in queue decisions."""

    def __init__(self, directory, fifo_flag, reg_engine, cuda_execution_mode):
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        self.item_path = os.path.join(directory, f"queue_profile_items_{stamp}.csv")
        self.scan_path = os.path.join(directory, f"queue_profile_scans_{stamp}.csv")
        self.item_file = open(self.item_path, "w", newline="", buffering=1024 * 1024)
        try:
            self.scan_file = open(self.scan_path, "w", newline="", buffering=1024 * 1024)
        except OSError:
            self.item_file.close()
            raise
        self.item_writer = csv.DictWriter(self.item_file, fieldnames=ITEM_FIELDS)
        self.scan_writer = csv.DictWriter(self.scan_file, fieldnames=SCAN_FIELDS)
        self.item_writer.writeheader()
        self.scan_writer.writeheader()
        self.common = {"fifo_flag": fifo_flag, "reg_engine": reg_engine,
                       "cuda_execution_mode": cuda_execution_mode}
        self.items = {}
        self.completed = []
        self.scan = None
        self.scan_number = 0
        self.empty_polls = 0
        self.empty_scan_ns = 0
        self.observed_total = 0
        self.last_directory_entry_count = 0
        self.registered_total = 0
        self.skipped_total = 0
        self.written = 0
        self.disabled = False

    def observe(self, filename, path):
        if self.disabled or filename in self.items:
            return
        perf_ns = time.perf_counter_ns()
        wall_ns = time.time_ns()
        try:
            mtime_ns = os.stat(path).st_mtime_ns
        except OSError:
            mtime_ns = None
        self.items[filename] = {
            **self.common, "pointer_filename": filename, "pointer_path": path,
            "file_mtime_ns": mtime_ns, "first_observed_wall_ns": wall_ns,
            "first_observed_perf_ns": perf_ns,
        }
        self.observed_total += 1

    def empty_scan(self, elapsed_ns):
        if not self.disabled:
            self.empty_polls += 1
            self.empty_scan_ns += elapsed_ns

    def directory_entries(self, count):
        if not self.disabled:
            self.last_directory_entry_count = count

    def begin_scan(self, filenames, scan_start_ns, scan_end_ns):
        if self.disabled:
            return
        self.scan_number += 1
        pointer_names = [f for f in filenames if f.endswith(".txt")]
        self.scan = {
            "scan_id": self.scan_number,
            "scan_start_perf_ns": scan_start_ns,
            "scan_end_perf_ns": scan_end_ns,
            "directory_entry_count": self.last_directory_entry_count,
            "candidate_file_count": len(filenames),
            "pointer_candidate_count": len(pointer_names),
            "selected_file_count": len(filenames),
            "lifo_old_pointer_count": 0,
            "scan_and_build_ms": (scan_end_ns - scan_start_ns) / 1e6,
            "empty_polls_since_previous": self.empty_polls,
            "empty_scan_ms_since_previous": self.empty_scan_ns / 1e6,
        }
        self.empty_polls = 0
        self.empty_scan_ns = 0
        for name in pointer_names:
            row = self.items.get(name)
            if row is not None:
                row.setdefault("first_eligible_perf_ns", scan_end_ns)
                row.update({"scan_id": self.scan_number,
                            "candidate_file_count": len(filenames),
                            "pointer_candidate_count": len(pointer_names)})

    def selection(self, start_ns, end_ns, filenames, newest):
        if self.disabled:
            return
        old_names = [f for f in filenames if f != newest and f.endswith(".txt")]
        self.scan["selected_file_count"] = 1
        self.scan["lifo_old_pointer_count"] = len(old_names)
        self.scan["newest_selection_ms"] = (end_ns - start_ns) / 1e6
        for name in filenames:
            row = self.items.get(name)
            if row is not None:
                row["selection_snapshot_perf_ns"] = end_ns
                row["lifo_old_pointer_count"] = len(old_names)

    def skip_handling(self, start_ns, end_ns):
        if not self.disabled:
            self.scan["skip_handling_ms"] = (end_ns - start_ns) / 1e6

    def mark(self, filename, field, value=None):
        if self.disabled:
            return
        row = self.items.get(filename)
        if row is not None:
            row[field] = time.perf_counter_ns() if value is None else value

    def identify(self, filename, volume, group, output_label=None):
        if self.disabled:
            return
        row = self.items.get(filename)
        if row is not None:
            row.update({"volume": volume, "group": group})
            if output_label is not None:
                row["output_label"] = output_label

    def decision(self, filename, outcome, reason=""):
        if self.disabled:
            return
        row = self.items.get(filename)
        if row is not None:
            row.update({"outcome": outcome, "skip_reason": reason,
                        "decision_perf_ns": time.perf_counter_ns(),
                        "decision_wall_ns": time.time_ns()})
            if outcome == "skipped_by_lifo":
                self.skipped_total += 1

    def complete(self, filename, outcome=None):
        if self.disabled:
            return
        row = self.items.get(filename)
        if row is not None:
            if outcome is not None:
                row["outcome"] = outcome
            row["item_complete_perf_ns"] = time.perf_counter_ns()
            self.completed.append(filename)
            if row.get("outcome") == "registered":
                self.registered_total += 1

    def end_scan(self, loop_ready_ns):
        if self.disabled or self.scan is None:
            return
        self.scan["loop_ready_perf_ns"] = loop_ready_ns
        self.scan["total_cycle_ms"] = (
            loop_ready_ns - self.scan["scan_start_perf_ns"]) / 1e6
        self.scan["observed_pointer_count_total"] = self.observed_total
        self.scan["registered_total"] = self.registered_total
        self.scan["skipped_total"] = self.skipped_total
        try:
            for name in self.completed:
                row = self.items.pop(name, None)
                if row is not None:
                    row["loop_ready_perf_ns"] = loop_ready_ns
                    self.item_writer.writerow(row)
                    self.written += 1
            self.scan_writer.writerow(self.scan)
            if self.written and self.written % 256 < len(self.completed):
                self.item_file.flush()
                self.scan_file.flush()
        except (OSError, ValueError) as error:
            self.disabled = True
            # Profiling failures must never change registration scheduling.
            import logging
            logging.warning("Queue profiling disabled after write failure: %s", error)
        self.completed = []
        self.scan = None

    def close(self):
        for stream in (self.item_file, self.scan_file):
            try:
                stream.close()
            except OSError:
                pass

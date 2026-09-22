# Queue processor profiling

Set `QUEUE_PROFILE_FLAG=on` for a controlled run. The default is `off`.
The launcher passes this variable into the container. Each acquisition writes
`queue_profile_items_<timestamp>.csv` and `queue_profile_scans_<timestamp>.csv`
beside its queue logfile in `/data`; preserve both files after consolidation.
The CSV files are separate from `registration_status.json`. They are buffered,
flushed periodically, and closed before a successful close-trigger acknowledgement.

## Clocks and interpretation

- Every `*_perf_ns` field uses `time.perf_counter_ns()`. Subtract these only
  within the same queue-processor process to calculate elapsed time.
- `first_observed_wall_ns` and `decision_wall_ns` use `time.time_ns()`.
  `file_mtime_ns` comes from `os.stat(...).st_mtime_ns` on the first scan that
  observed the pointer. Subtract wall values from each other only. File mtime
  is an availability proxy, subject to timestamp precision and clock behavior
  of the source filesystem; it is not necessarily image arrival time.
- Do not subtract a wall or file timestamp from a perf-counter timestamp.
- Profiling calls and one additional `stat` per newly observed pointer add
  some work when enabled. `scan_and_build_ms` and `total_cycle_ms` include this
  profiling overhead; use an off/on control run to assess impact.

## Item CSV

One row is emitted for each observed pointer that reaches an existing terminal
processing decision. `scan_id` links it to a scan row. `volume`, `group`, and
`output_label` come from the queue processor's existing parsed identifiers.
`pointer_filename` and `pointer_path` identify the source pointer.
`fifo_flag`, `reg_engine`, and `cuda_execution_mode` record runtime settings.

`candidate_file_count` is the count of unseen valid `.json`, `.txt`, and
`.closeQ` files in the scan snapshot. `pointer_candidate_count` counts only
`.txt` files in that snapshot. `lifo_old_pointer_count` counts pointer files
discarded by the existing newest-file selection. These are snapshot counts,
**not** the depth of a persistent queue.

| Field | Event or meaning |
|---|---|
| `file_mtime_ns` | Pointer mtime on first observation; wall clock domain. |
| `first_observed_wall_ns`, `first_observed_perf_ns` | First successful directory scan containing the unseen pointer. |
| `first_eligible_perf_ns` | First candidate snapshot containing the pointer. |
| `selection_snapshot_perf_ns` | Existing LIFO newest-file `max` call finished; blank for FIFO. |
| `selected_perf_ns` | File entered the existing processing loop; blank for LIFO-discarded pointers. |
| `decision_wall_ns`, `decision_perf_ns` | Existing register/reference/skip/refusal decision. |
| `preparation_start_perf_ns`, `preparation_end_perf_ns` | Pointer write-stability wait, filename parsing, and target-list read. For a skipped pointer, end also includes expected-group validation. |
| `reference_setup_start_perf_ns`, `reference_setup_end_perf_ns` | Initial upsampling/identity setup; only the reference pointer. |
| `registration_start_perf_ns`, `registration_end_perf_ns` | Immediately around the existing `run_MIregistration` call; blank for skips/reference. |
| `status_start_perf_ns`, `status_end_perf_ns` | Existing status record and possible volume finalization/publication. |
| `reference_update_start_perf_ns`, `reference_update_end_perf_ns` | Existing post-registration reference check; blank for skips/reference. |
| `item_complete_perf_ns` | Selected pointer's status/reference work and seen-file marking completed, or a LIFO skip finished; before normal progress logs. |
| `loop_ready_perf_ns` | Whole candidate batch finished, after the selected item and progress logs, before the next scan. A skipped pointer therefore shares this timestamp with the selected file from its batch. |
| `outcome` | `registered`, `failed`, `reference`, `failed_reference`, `skipped_by_lifo`, or an existing refusal/disappearance path. |
| `skip_reason` | Literal existing skip message for LIFO skips. |
| `status_result` | `record_attempted` if the status call returned, or `error` if skipped-pointer status handling raised; this is not a guarantee of durable publication. |

Useful differences: `first_observed_wall_ns - file_mtime_ns` estimates
availability-to-discovery; `selected_perf_ns - first_eligible_perf_ns` is
discovery-to-selection; `registration_start_perf_ns - selected_perf_ns` is
selection-to-registration; `registration_end_perf_ns - registration_start_perf_ns`
is registration duration; `loop_ready_perf_ns - registration_end_perf_ns`
includes post-registration handling and progress logging. For skips, compare
`decision_wall_ns - file_mtime_ns` and `decision_perf_ns - first_observed_perf_ns`.
Blank fields mean that stage did not apply or was not observable.

## Scan CSV

One row describes a nonempty call to `list_new_files` and the work done before
the next scan. `scan_start_perf_ns`, `scan_end_perf_ns`, and
`loop_ready_perf_ns` bound the cycle. `scan_and_build_ms` includes `os.listdir`,
extension/seen filtering, per-candidate `getmtime`, acquisition-order sorting,
and enabled profiling overhead. `newest_selection_ms` times only the LIFO
`max(...getmtime...)` operation. `skip_handling_ms` covers the existing
old-file skip loop, including pointer reading and status work.
`total_cycle_ms` runs from scan start to loop readiness.

`directory_entry_count` is the size of the same `os.listdir` result already
used by the processor; it includes unrelated files in the directory.
`candidate_file_count`, `pointer_candidate_count`, `selected_file_count`, and
`lifo_old_pointer_count` describe that one snapshot. `empty_polls_since_previous`
and `empty_scan_ms_since_previous` aggregate empty scans since the preceding
nonempty row; the 5 ms sleeps are **not** included in `empty_scan_ms`.
`observed_pointer_count_total`, `registered_total`, and `skipped_total` are
profiling counters, not replacements for `registration_status.json`.

Normal acquisition completion closes the CSV files before the close trigger is
removed. Abrupt termination before acquisition close may leave the latest
buffered rows unwritten; partial rows are flushed approximately every 256
completed pointers.

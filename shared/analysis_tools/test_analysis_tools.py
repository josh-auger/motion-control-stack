"""Focused synthetic tests for the reusable post-run analyzers."""

from __future__ import annotations

import csv
import json
import os
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from analysis_common import (  # noqa: E402
    AnalysisError, adaptive_break_limits, describe, discover_queue_profile_pair,
    discover_unique, optional_files, ordinary_least_squares, read_csv, read_json,
    resolve_output_directory, running_median, temporal_third, validate_input_directory,
)
from analyze_fire_moco_profile import (  # noqa: E402
    _plots as plot_fire_moco, discovery_medians_by_volume,
)
from analyze_motion_dashboard_profile import dashboard_x_axis  # noqa: E402
from analyze_queue_profile import (  # noqa: E402
    _plot_scans as plot_queue_scans, attach_scan_volumes,
    match_receipts_to_queue_items, parse_items, scan_discovery_medians_by_volume,
)
from analyze_registration_completeness import (  # noqa: E402
    OUTCOME_CODES, OUTCOME_COLORS, _plot as plot_completeness,
    build_completeness_figure, outcome_matrix, parse_status,
)
from analyze_registration_calls import (  # noqa: E402
    _plots as plot_registration_calls, log_only_rows, parse_log_components,
    split_reference_initialization,
)

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")


def write_csv_fixture(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class CommonTests(unittest.TestCase):
    def test_valid_and_nonexistent_input_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(validate_input_directory(Path(directory)), Path(directory).resolve())
            with self.assertRaisesRegex(AnalysisError, "does not exist"):
                validate_input_directory(Path(directory) / "missing")

    def test_output_directory_contract(self):
        with tempfile.TemporaryDirectory() as directory, tempfile.TemporaryDirectory() as absolute:
            source = Path(directory)
            root, leaf = resolve_output_directory(source, Path("analysis_test1"), "queue_profile")
            self.assertEqual(root, source / "analysis_test1")
            self.assertEqual(leaf, source / "analysis_test1" / "queue_profile")
            root, leaf = resolve_output_directory(source, Path(absolute), "queue_profile")
            self.assertEqual(root, Path(absolute))
            self.assertEqual(leaf, Path(absolute) / "queue_profile")
            root, _ = resolve_output_directory(source, None, "queue_profile", now=datetime(2026, 9, 30, 12, 15, 3))
            self.assertEqual(root.name, "analysis_20260930_121503")

    def test_temporal_thirds_for_nondivisible_and_short_datasets(self):
        self.assertEqual([temporal_third(i, 5) for i in range(5)], ["early", "early", "middle", "middle", "late"])
        self.assertEqual([temporal_third(i, 2) for i in range(2)], ["early", "middle"])
        self.assertEqual(describe([])["count"], 0)

    def test_running_median_and_linear_fit(self):
        medians, window = running_median([1, 9, 3, 5, 7], window=3)
        self.assertEqual(window, 3)
        self.assertEqual(medians, [5.0, 3.0, 5.0, 5.0, 6.0])
        fit = ordinary_least_squares([1, 2, 3], [3, 5, 7])
        self.assertAlmostEqual(fit["slope"], 2.0)
        self.assertAlmostEqual(fit["intercept"], 1.0)
        self.assertAlmostEqual(fit["r_squared"], 1.0)

    def test_adaptive_broken_axis_rule(self):
        self.assertIsNone(adaptive_break_limits(list(range(1, 31))))
        broken = adaptive_break_limits([1, 2, 2, 3, 3, 4, 4, 5, 100])
        self.assertIsNotNone(broken)
        self.assertEqual(broken["upper_count"], 1)

    def test_missing_and_ambiguous_artifact(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(AnalysisError, "was not found"):
                discover_unique(root, "*.csv", "profile")

    def test_profiler_disabled_hint_and_missing_optional_artifact(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(AnalysisError, "MOTION_DASHBOARD_PROFILE_FLAG"):
                discover_unique(
                    root, "motion_dashboard_profile_*.csv", "motion-dashboard profile",
                    missing_hint="MOTION_DASHBOARD_PROFILE_FLAG may have been disabled.",
                )
            self.assertEqual(optional_files(root, "fire_moco_state.moco"), [])
            (root / "a.csv").touch()
            (root / "b.csv").touch()
            with self.assertRaisesRegex(AnalysisError, "Ambiguous"):
                discover_unique(root, "*.csv", "profile")

    def test_queue_pair_discovery_and_ambiguity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "queue_profile_items_a.csv").touch()
            (root / "queue_profile_scans_a.csv").touch()
            self.assertEqual(discover_queue_profile_pair(root)[0].name, "queue_profile_items_a.csv")
            (root / "queue_profile_items_b.csv").touch()
            (root / "queue_profile_scans_b.csv").touch()
            with self.assertRaisesRegex(AnalysisError, "Ambiguous"):
                discover_queue_profile_pair(root)

    def test_required_column_and_malformed_csv(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "test.csv"
            path.write_text("a\n1\n", encoding="utf-8")
            with self.assertRaisesRegex(AnalysisError, "column"):
                read_csv(path, {"a", "b"})
            path.write_text('a,b\n"unterminated,2\n', encoding="utf-8")
            with self.assertRaisesRegex(AnalysisError, "Could not read CSV"):
                read_csv(path, {"a", "b"})

    def test_malformed_json(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text("{", encoding="utf-8")
            with self.assertRaisesRegex(AnalysisError, "Could not read JSON"):
                read_json(path)


class QueueTests(unittest.TestCase):
    def test_variable_acquisition_and_group_identifiers(self):
        fields = sorted({
            "scan_id", "volume", "group", "pointer_filename", "outcome",
            "first_observed_perf_ns", "first_eligible_perf_ns", "selected_perf_ns",
            "decision_perf_ns", "item_complete_perf_ns", "registration_start_perf_ns",
            "registration_end_perf_ns", "candidate_file_count", "pointer_candidate_count",
        })
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "items.csv"
            rows = []
            for index, (volume, group) in enumerate([(3, 0), (3, 7), (9, 2), (41, 99)], 1):
                rows.append({
                    "scan_id": index, "volume": volume, "group": group,
                    "pointer_filename": f"p_{volume}_{group}.txt", "outcome": "registered",
                    "first_observed_perf_ns": index * 1000, "first_eligible_perf_ns": index * 1000,
                    "selected_perf_ns": index * 1000 + 1, "decision_perf_ns": index * 1000 + 2,
                    "registration_start_perf_ns": index * 1000 + 3,
                    "registration_end_perf_ns": index * 1000 + 8,
                    "item_complete_perf_ns": index * 1000 + 10,
                    "candidate_file_count": 1, "pointer_candidate_count": 1,
                })
            write_csv_fixture(path, fields, rows)
            parsed = parse_items(path)
            self.assertEqual([(r["volume"], r["group"]) for r in parsed], [(3, 0), (3, 7), (9, 2), (41, 99)])

    def test_exact_fire_queue_matching_and_unmatched_records(self):
        receipts = [{"image_filename": "a.nhdr"}, {"image_filename": "orphan.nhdr"}]
        queue = [{"pointer_filename": "p1.txt"}, {"pointer_filename": "p2.txt"}]
        result = match_receipts_to_queue_items(receipts, queue, {"p1.txt": ["a.nhdr"], "p2.txt": ["missing.nhdr"]})
        self.assertEqual(len(result["matched"]), 1)
        self.assertEqual(result["unmatched_fire"][0]["image_filename"], "orphan.nhdr")
        self.assertEqual(result["unmatched_queue"][0]["pointer_filename"], "p2.txt")

    def test_duplicate_identifier_and_timestamp_only_refusal(self):
        with self.assertRaisesRegex(AnalysisError, "Duplicate/ambiguous FIRE"):
            match_receipts_to_queue_items([{"image_filename": "a"}, {"image_filename": "a"}], [{"pointer_filename": "p"}], {"p": ["a"]})
        with self.assertRaisesRegex(AnalysisError, "timestamp-only"):
            match_receipts_to_queue_items([{"timestamp": 1}], [{"pointer_filename": "p"}], {"p": ["a"]})
        with self.assertRaisesRegex(AnalysisError, "pointer membership"):
            match_receipts_to_queue_items([{"image_filename": "a"}], [{"pointer_filename": "p1"}, {"pointer_filename": "p2"}], {"p1": ["a"], "p2": ["a"]})

    def test_scan_volume_mapping_and_separate_plots(self):
        scans = [
            {"scan_id": 1, "directory_entry_count": 10, "scan_and_build_ms": 21.0},
            {"scan_id": 2, "directory_entry_count": 20, "scan_and_build_ms": 41.0},
            {"scan_id": 3, "directory_entry_count": 30, "scan_and_build_ms": 61.0},
        ]
        items = [
            {"scan_id": 1, "volume": 2},
            {"scan_id": 2, "volume": 3},
            {"scan_id": 2, "volume": 4},
            {"scan_id": 3, "volume": 5},
        ]
        self.assertEqual(attach_scan_volumes(scans, items), 0)
        self.assertEqual([row["volume"] for row in scans], [2, 4, 5])
        self.assertEqual(
            scan_discovery_medians_by_volume([
                *scans,
                {"volume": 4, "scan_and_build_ms": 45.0},
            ]),
            [(2, 21.0), (4, 43.0), (5, 61.0)],
        )
        with tempfile.TemporaryDirectory() as directory:
            info = plot_queue_scans(scans, Path(directory))
            self.assertTrue((Path(directory) / "queue_discovery_timing.png").is_file())
            self.assertTrue((Path(directory) / "queue_directory_entries.png").is_file())
            self.assertTrue((Path(directory) / "queue_discovery_vs_directory_entries.png").is_file())
            self.assertAlmostEqual(info["fit"]["slope"], 2.0)
            self.assertAlmostEqual(info["fit"]["intercept"], 1.0)
            self.assertEqual(info["volume_medians"], [(2, 21.0), (4, 41.0), (5, 61.0)])


class FireMocoPlotTests(unittest.TestCase):
    def test_volume_axis_medians_and_removed_consumed_plot(self):
        rows = [
            {"volume": 2, "discovery_selection_ms": 4.0},
            {"volume": 1, "discovery_selection_ms": 3.0},
            {"volume": 1, "discovery_selection_ms": 1.0},
            {"volume": 2, "discovery_selection_ms": 8.0},
        ]
        self.assertEqual(discovery_medians_by_volume(rows), [(1, 2.0), (2, 6.0)])
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            stale = output / "fire_moco_lookup_timing.png"
            stale.touch()
            info = plot_fire_moco(rows, output)
            self.assertEqual(info["x_field"], "volume")
            self.assertTrue((output / "fire_moco_discovery_timing.png").is_file())
            self.assertFalse(stale.exists())
            self.assertFalse((output / "consumed_registration_index.png").exists())


class DashboardPlotTests(unittest.TestCase):
    def test_volume_axis_is_preferred(self):
        values, label, mapping = dashboard_x_axis([
            {"volume": 4, "profile_start_perf_ns": 100},
            {"volume": 8, "profile_start_perf_ns": 200},
        ])
        self.assertEqual(values, [4.0, 8.0])
        self.assertEqual(label, "Volume number")
        self.assertEqual(mapping, "volume")

    def test_runtime_fallback_requires_explicit_timestamps(self):
        values, label, mapping = dashboard_x_axis([
            {"volume": None, "profile_start_perf_ns": 1_000_000_000},
            {"volume": None, "profile_start_perf_ns": 3_500_000_000},
        ])
        self.assertEqual(values, [0.0, 2.5])
        self.assertEqual(label, "Run time (s)")
        self.assertIn("profile_start_perf_ns", mapping)
        with self.assertRaisesRegex(AnalysisError, "cannot be established"):
            dashboard_x_axis([{"volume": None}])


class RegistrationCompletenessTests(unittest.TestCase):
    @staticmethod
    def _row(volume, statuses):
        group_outcomes = dict(enumerate(statuses))
        return {
            "volume": volume,
            "expected": len(statuses),
            "registered": statuses.count("registered"),
            "skipped": statuses.count("skipped"),
            "failed": statuses.count("failed"),
            "registration_finalized": True,
            "group_outcomes": group_outcomes,
        }

    def test_registered_skipped_failed_variable_shapes(self):
        payload = {
            "schema_version": 1, "revision": 2,
            "volumes": {
                "4": {"registration_finalized": True, "groups": {
                    "0": {"status": "registered", "transform_filename": "x.tfm"},
                    "5": {"status": "skipped"}, "11": {"status": "failed"},
                }},
                "20": {"registration_finalized": True, "groups": {
                    "2": {"status": "registered", "transform_filename": "y.tfm"},
                }},
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "registration_status.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            volumes, exceptions = parse_status(path)
            self.assertEqual([r["volume"] for r in volumes], [4, 20])
            self.assertEqual(volumes[0]["expected"], 3)
            self.assertEqual(
                volumes[0]["group_outcomes"],
                {0: "registered", 5: "skipped", 11: "failed"},
            )
            self.assertEqual({r["status"] for r in exceptions}, {"skipped", "failed"})

    def test_empty_status_and_invalid_schema_are_handled(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "registration_status.json"
            path.write_text(json.dumps({"schema_version": 1, "volumes": {}}), encoding="utf-8")
            self.assertEqual(parse_status(path), ([], []))
            path.write_text(json.dumps({"schema_version": 99, "volumes": {}}), encoding="utf-8")
            with self.assertRaisesRegex(AnalysisError, "schema_version"):
                parse_status(path)

    def test_positional_outcome_matrix_and_discrete_rendering(self):
        rows = [
            self._row(0, ["registered", "registered", "registered", "registered"]),
            self._row(1, ["skipped", "registered", "registered", "registered"]),
            self._row(2, ["registered", "registered", "registered", "skipped"]),
            self._row(3, ["registered", "registered", "skipped", "registered"]),
            self._row(4, ["registered", "skipped", "registered", "skipped"]),
            self._row(5, ["registered", "registered", "failed", "registered"]),
        ]
        matrix, volumes, groups = outcome_matrix(rows)
        r, s, f = (OUTCOME_CODES[status] for status in ("registered", "skipped", "failed"))
        self.assertEqual(volumes, [0, 1, 2, 3, 4, 5])
        self.assertEqual(groups, [0, 1, 2, 3])
        self.assertEqual(matrix.tolist(), [
            [r, s, r, r, r, r],
            [r, r, r, r, s, r],
            [r, r, r, s, r, f],
            [r, r, s, r, s, r],
        ])

        fig, ax, secondary, expected, legend, mesh = build_completeness_figure(rows)
        fig.canvas.draw()
        self.assertEqual(expected, 4)
        self.assertIsNotNone(secondary)
        ticks = [tick for tick in ax.get_yticks() if ax.get_ylim()[0] <= tick <= ax.get_ylim()[1]]
        self.assertTrue(all(float(tick).is_integer() for tick in ticks))
        self.assertAlmostEqual(secondary._functions[0](4), 100.0)
        self.assertEqual(ax.get_ylim(), (0.0, 4.0))
        self.assertEqual(mesh._shading, "flat")
        from matplotlib.collections import QuadMesh
        from matplotlib.colors import to_rgba
        self.assertIsInstance(mesh, QuadMesh)
        for status, code in OUTCOME_CODES.items():
            self.assertEqual(mesh.cmap(mesh.norm(code)), to_rgba(OUTCOME_COLORS[status]))
        self.assertEqual(legend._ncols, 3)
        self.assertLess(legend.get_bbox_to_anchor()._bbox.y0, 0)
        self.assertEqual([text.get_text() for text in legend.get_texts()], [
            "registered (18)", "skipped (5)", "failed (1)",
        ])
        fig.clf()
        with tempfile.TemporaryDirectory() as directory:
            info = plot_completeness(rows, Path(directory) / "plot.png")
            self.assertTrue(info["percent_axis"])
            self.assertTrue((Path(directory) / "plot.png").is_file())

    def test_variable_expected_counts_omit_percent_axis(self):
        rows = [
            {
                "volume": 0, "expected": 1, "registered": 1,
                "skipped": 0, "failed": 0, "registration_finalized": True,
                "group_outcomes": {22: "registered"},
            },
            self._row(1, ["registered", "registered", "registered", "registered"]),
        ]
        matrix, volumes, groups = outcome_matrix(rows)
        self.assertEqual((volumes, groups), ([0, 1], list(range(23))))
        self.assertTrue(np.isnan(matrix[0, 0]))
        self.assertEqual(matrix[22, 0], OUTCOME_CODES["registered"])
        fig, _ax, secondary, expected, _legend, _mesh = build_completeness_figure(rows)
        self.assertIsNone(expected)
        self.assertIsNone(secondary)
        fig.clf()


class RegistrationCallTests(unittest.TestCase):
    def test_persistent_and_standalone_components_use_explicit_labels(self):
        text = "\n".join([
            "Running registration call 1",
            "Alignment transform saved as /data/alignTransform_0001_0002-0003.tfm",
            "[TIMING] target_read_ms: 1.25",
            "[TIMING] optimizer_ms: 4.50",
            "Running registration call 2",
            "SLIMM LOAD_REFERENCE completed in 0.210 s: SLIMM_PROTOCOL SUCCESS "
            "command=LOAD_REFERENCE total_ms=209.5 read_ms=4.0",
            "SLIMM REGISTER 0002_0002-0004 completed in 0.006 s: "
            "SLIMM_PROTOCOL SUCCESS command=REGISTER label=0002_0002-0004 "
            "request_ms=5.9 optimizer_ms=4.1 target_read_ms=0.8 target_add_ms=0.02",
        ])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "log_local_queue_processor_test.log"
            path.write_text(text, encoding="utf-8")
            records, warnings = parse_log_components([path])
            self.assertFalse(warnings)
            self.assertEqual(records["0001_0002-0003"]["standalone_optimizer_ms"], 4.5)
            self.assertEqual(records["0002_0002-0004"]["persistent_request_ms"], 5.9)
            self.assertEqual(records["0002_0002-0004"]["persistent_reference_initialization_ms"], 209.5)

    def test_log_components_do_not_match_by_timestamp_or_call_order(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "log.log"
            path.write_text("2026-01-01 - [TIMING] optimizer_ms: 3.0\n", encoding="utf-8")
            records, _ = parse_log_components([path])
            self.assertEqual(records, {})

    def test_log_only_fallback_preserves_explicit_volume_group_and_backend(self):
        components = {
            "0012_0007-0009": {
                "log_full_registration_elapsed_ms": 6.2,
                "persistent_request_ms": 5.8,
            }
        }
        rows = log_only_rows(components)
        self.assertEqual((rows[0]["registration_index"], rows[0]["volume"], rows[0]["group"]), (12, 7, 9))
        self.assertEqual((rows[0]["reg_engine"], rows[0]["cuda_execution_mode"]), ("cuda", "persistent"))

    def test_reference_initialization_excluded_only_when_explicit(self):
        ordinary = {"output_label": "0002_0001-0001", "volume": 1, "full_registration_elapsed_ms": 6.0, "persistent_request_ms": 5.0}
        initialization = {"output_label": "0001_0001-0000", "volume": 1, "full_registration_elapsed_ms": 220.0, "persistent_request_ms": 6.0, "persistent_reference_initialization_ms": 210.0}
        included, excluded = split_reference_initialization([initialization, ordinary])
        self.assertEqual(included, [ordinary])
        self.assertEqual(excluded, [initialization])
        large_but_unmarked = {"output_label": "0003_0001-0002", "volume": 1, "full_registration_elapsed_ms": 500.0}
        included, excluded = split_reference_initialization([ordinary, large_but_unmarked])
        self.assertEqual(len(included), 2)
        self.assertFalse(excluded)

    def test_backend_plot_counts_every_individual_component_observation(self):
        rows = [
            {"output_label": "1_1-0", "volume": 1, "full_registration_elapsed_ms": 6.0, "persistent_request_ms": 5.0, "persistent_optimizer_ms": 3.0},
            {"output_label": "2_1-1", "volume": 1, "full_registration_elapsed_ms": 7.0, "persistent_request_ms": 6.0, "persistent_optimizer_ms": 4.0},
        ]
        with tempfile.TemporaryDirectory() as directory:
            info = plot_registration_calls(
                rows, ["persistent_request_ms", "persistent_optimizer_ms"], Path(directory)
            )
            self.assertEqual(info["component_point_count"], 4)
            self.assertTrue((Path(directory) / "registration_timing_by_volume.png").is_file())


if __name__ == "__main__":
    unittest.main()

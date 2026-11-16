"""Numerical, acquisition, missing-data, and artifact checks for retrospective motion."""

from __future__ import annotations

import ast
import base64
import csv
import hashlib
import html
import json
from pathlib import Path
import re

import h5py
import ismrmrd
import numpy as np
import pytest
import SimpleITK as sitk

from shared.analysis_tools import compare_subject_motion as motion


def euler(angles=(0, 0, 0), translation=(0, 0, 0), center=(10, -20, 30)):
    transform = sitk.Euler3DTransform()
    transform.SetCenter(center)
    transform.SetRotation(*angles)
    transform.SetTranslation(translation)
    return transform


@pytest.fixture(scope="module")
def monitor_functions():
    # Compile only pure arithmetic functions from the existing monitor, avoiding
    # its service imports, logging, global state, and plotting setup.
    path = motion.HERE.parents[1] / "apps/motion_monitor/monitor_directory.py"
    tree = ast.parse(path.read_text())
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                 and node.name in {"compose_transform_pair", "calculate_displacement"}]
    assert len(functions) == 2
    namespace = {"np": np, "sitk": sitk}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)
    return namespace


@pytest.mark.parametrize("axis", range(3))
def test_known_single_axis_rotation_and_translation(axis):
    angles = [0., 0., 0.]
    angles[axis] = np.deg2rad(12)
    tfm = euler(angles, (3, 4, 0))
    assert motion.displacement(tfm) == pytest.approx(5 + 100 * np.sin(np.deg2rad(6)))
    params = motion.motion_parameters(tfm)
    assert params["translation_magnitude_mm"] == 5
    assert params["rotation_euler_norm_deg"] == pytest.approx(12)
    assert params["rotation_angle_deg"] == pytest.approx(12)


def test_zero_tiny_and_invalid_radius():
    assert motion.displacement(euler()) == 0
    assert motion.displacement(euler((1e-12, 0, 0))) == pytest.approx(0, abs=1e-9)
    assert motion.displacement(euler(translation=(3, 4, 12))) == 13
    for radius in (0, -1, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            motion.displacement(euler(), radius)


def test_composition_matches_physical_point_mapping_and_centers():
    previous = euler((.12, -.21, .17), (2, -3, 4), (5, 7, -9))
    current = euler((-.16, .23, -.18), (-4, 1, 3), (-2, 10, 11))
    relative = motion.relative_transform(previous, current)
    assert relative.GetCenter() == previous.GetCenter()
    assert relative.GetComputeZYX() is False
    for point in [(0., 0., 0.), (12., 5., -8.), (-7., 22., 14.)]:
        expected = current.TransformPoint(previous.GetInverse().TransformPoint(point))
        assert relative.TransformPoint(point) == pytest.approx(expected, abs=1e-12)
    # Same rigid mapping expressed about another center has zero relative motion.
    recentered = euler(center=(20, -10, 12))
    recentered.SetMatrix(previous.GetMatrix())
    r = np.array(previous.GetMatrix()).reshape(3, 3)
    old, new = np.array(previous.GetCenter()), np.array(recentered.GetCenter())
    recentered.SetTranslation((np.array(previous.GetTranslation()) + old - r @ old - new + r @ new).tolist())
    assert motion.displacement(motion.relative_transform(previous, recentered)) == pytest.approx(0, abs=1e-6)


def test_multiaxis_calculations_match_established_monitor(monitor_functions):
    rng = np.random.default_rng(19)
    for _ in range(25):
        a, b = [euler(rng.uniform(-.3, .3, 3).tolist(), rng.normal(size=3).tolist(), rng.normal(size=3).tolist()) for _ in range(2)]
        actual = motion.relative_transform(a, b)
        reference = monitor_functions["compose_transform_pair"](a, b)
        assert actual.GetParameters() == pytest.approx(reference.GetParameters(), abs=1e-12)
        assert motion.displacement(actual) == pytest.approx(monitor_functions["calculate_displacement"](reference, 50), abs=1e-11)


def test_saved_compute_zyx_and_versor_physical_conventions(tmp_path):
    tfm = euler((.1, -.2, .3), (1, 2, -3))
    tfm.SetComputeZYX(True)
    sitk.WriteTransform(tfm, str(tmp_path / "saved.tfm"))
    loaded, _ = motion.load_transform(tmp_path, "saved.tfm", [])
    assert loaded.GetComputeZYX() is True
    params = motion.motion_parameters(loaded)
    assert [params[f"{a}_rotation_deg"] for a in "XYZ"] == pytest.approx(np.degrees(tfm.GetParameters()[:3]))
    versor = sitk.VersorRigid3DTransform()
    versor.SetMatrix(tfm.GetMatrix())
    versor.SetTranslation(tfm.GetTranslation())
    versor.SetCenter(tfm.GetCenter())
    assert motion.motion_parameters(versor)["rotation_angle_deg"] == pytest.approx(params["rotation_angle_deg"])
    assert motion.displacement(motion.relative_transform(versor, tfm)) == pytest.approx(0, abs=1e-6)


def test_missing_observations_never_bridge_gaps():
    tfms = [euler(), euler(translation=(1, 0, 0)), None, euler(translation=(3, 0, 0)), euler(translation=(4, 0, 0))]
    assert motion.consecutive_displacements(tfms) == [None, 1, None, None, 1]
    assert motion.consecutive_displacements([]) == []
    assert all(v is None for v in motion.motion_parameters(None).values())


def test_statistics_missing_values_events_and_percentiles():
    st = motion.metric_statistics([None, .1, .4, .5, None, .6, .3], .3)
    assert st["count"] == 5
    assert st["mean"] == pytest.approx(.38)
    assert st["median"] == .4
    assert st["p95"] == pytest.approx(.58)
    assert st["max"] == .6
    assert st["exceedance_count"] == 3
    assert st["exceedance_pct"] == 60
    assert st["event_runs"] == 2
    assert st["missing_count"] == 2
    assert motion.metric_statistics([None], .3)["exceedance_pct"] is None


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_subject(root: Path, subject="GBS001", nvol=4, final_offset=.91):
    directory = root / f"savedData_20261110T120000_protocol_{subject}"
    directory.mkdir()
    processed = directory / "processed_protocol"
    archive = processed / "transforms"
    archive.mkdir(parents=True)
    offsets = [i * .07 for i in range(13)] + [final_offset]
    metadata = {
        "encoding": [{"encodingLimits": {"repetition": {"minimum": 0, "maximum": nvol - 1}, "slice": {"minimum": 0, "maximum": 55}},
                      "reconSpace": {"matrixSize": {"x": 90, "y": 90, "z": 1}, "fieldOfView_mm": {"x": 216., "y": 216., "z": 2.4}}}],
        "sequenceParameters": {"TR": [1000.], "TE": [29.]},
        "userParameters": {"userParameterLong": [{"name": "MultiBandFactor", "value": 4}], "userParameterDouble": [{"name": "SpacingBetweenSlices", "value": 2.4}]},
        # Deliberately unordered spatial keys; arrival groups are sorted in time.
        "SliceTiming": {str(i): offsets[(i * 5) % 14] for i in range(56)},
    }
    (directory / "protocol_metadata.json").write_text(json.dumps(metadata))
    (directory / "protocol_volume_0000_reference.nhdr").write_text("NRRD0004\nspace: left-posterior-superior\n")
    status = {"schema_version": 1, "volumes": {}}
    queue = []
    for v, g in [(0, 14), *((v, g) for v in range(1, nvol) for g in range(14))]:
        name = f"alignTransform_{0 if v == 0 else (v-1)*14+g+1:04d}_{v:04d}-{g:04d}" + ("_identity" if v == 0 else "") + ".tfm"
        tfm = sitk.VersorRigid3DTransform() if v == 0 else euler((g * .001, v * .002, -g * .003), (v + g * .02, 0, 0))
        if v == 0:
            tfm.SetCenter((10, -20, 30))
        sitk.WriteTransform(tfm, str((directory if v == 0 else archive) / name))
        status["volumes"].setdefault(str(v), {"registration_finalized": True, "groups": {}})["groups"][str(g)] = {"status": "registered", "transform_filename": name}
        if v:
            (processed / f"protocol_volume_{v:04d}_group_{g:04d}.txt").write_text("\n".join(f"protocol_volume_{v:04d}_slice_{i:04d}.nhdr" for i in range(g * 4, (g + 1) * 4)) + "\n")
        queue.append({"volume": v, "group": g, "outcome": "reference" if v == 0 else "registered", "output_label": name.removeprefix("alignTransform_").removesuffix(".tfm"),
                      "reg_engine": "cuda", "cuda_execution_mode": "persistent", "fifo_flag": "off",
                      "registration_start_perf_ns": str((v * 14 + g) * 100_000_000) if v else "",
                      "registration_end_perf_ns": str((v * 14 + g) * 100_000_000 + 8_000_000 + g * 1000) if v else ""})
    (directory / "registration_status.json").write_text(json.dumps(status))
    write_csv(directory / "queue_profile_items_test.csv", queue)
    with h5py.File(directory / "measurement-test.hdf5", "w") as stream:
        for i in range(nvol * 56):
            v, ordinal = divmod(i, 56)
            header = ismrmrd.ImageHeader()
            header.repetition = v
            header.slice = ordinal
            header.acquisition_time_stamp = 1000 + round(400 * (v + offsets[ordinal // 4]))
            group = stream.create_group(f"image_{i}")
            group.create_dataset("header", data=np.frombuffer(bytes(header), dtype=np.uint8))
    return directory


def snapshot(directory):
    return {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.rglob("*") if p.is_file()}


def assert_no_obsolete_metrics(value):
    text = (value if isinstance(value, str) else json.dumps(value)).lower()
    for obsolete in ["volume_fd", "volume_statistics", "volume_rows", "representative_group", "volume_motion", "end-of-volume", "volume-level fd"]:
        assert obsolete not in text


def displacement_summary_table(report):
    section = report.split("<h2>Consecutive slice-group displacement: all expected registration estimates</h2>", 1)[1]
    table = section.split("</table>", 1)[0]
    headers = [html.unescape(value) for value in re.findall(r"<th>(.*?)</th>", table)]
    rows = [[html.unescape(value) for value in re.findall(r"<td>(.*?)</td>", row)]
            for row in re.findall(r"<tr>(.*?)</tr>", table) if "<td>" in row]
    return headers, rows


def test_subject_geometry_counts_group_sequence_and_read_only(tmp_path):
    directory = make_subject(tmp_path)
    before = snapshot(directory)
    s = motion.analyze_subject("GBS001", directory)
    assert snapshot(directory) == before
    assert len(s["group_rows"]) == 43
    assert "volume_rows" not in s and "volume_statistics" not in s
    assert s["processing"]["expected_registrations"] == 42
    assert s["processing"]["registered"] == 42
    assert s["processing"]["reference_records"] == 1
    assert s["processing"]["status_record_count"] == 43
    assert s["processing"]["verified_pointers"] == 42
    assert s["processing"]["complete_registered_volumes_excluding_reference"] == 3
    assert s["metadata"]["observed_mrd_slice_images"] == 224
    assert s["runtime_statistics_ms"]["count"] == 42
    assert s["group_statistics"]["group_displacement_mm"]["count"] == 42
    assert [(r["volume"], r["group"]) for r in s["group_rows"]] == [(0, 14), *((v, g) for v in range(1, 4) for g in range(14))]
    assert s["group_rows"][0]["acquisition_time_s"] == .91
    assert s["group_rows"][0]["group_displacement_mm"] is None
    assert s["group_rows"][1]["acquisition_time_s"] == 1
    # Both volume-boundary transitions survive; none are reset to an identity.
    for v in [2, 3]:
        row = next(r for r in s["group_rows"] if r["volume"] == v and r["group"] == 0)
        previous = euler((.013, (v - 1) * .002, -.039), (v - 1 + .26, 0, 0))
        current = euler((0, v * .002, 0), (v, 0, 0))
        expected = motion.displacement(motion.relative_transform(previous, current))
        assert row["group_displacement_mm"] == pytest.approx(expected)
    # Reference is excluded, with exactly one undefined baseline observation.
    st = s["group_statistics"]["group_displacement_mm"]
    assert st["missing_count"] == 1
    assert st["mean"] == pytest.approx(np.mean([r["group_displacement_mm"] for r in s["group_rows"][1:]]))
    assert s["warnings"] == ["Motion CSV missing: recorded group-displacement cross-check unavailable."]


def test_failed_skipped_absent_and_unreadable_observations(tmp_path):
    directory = make_subject(tmp_path)
    path = directory / "registration_status.json"
    status = json.loads(path.read_text())
    groups = status["volumes"]["2"]["groups"]
    groups["10"]["status"] = "failed"
    groups["11"]["status"] = "skipped"
    del groups["12"]
    filename = groups["13"]["transform_filename"]
    (directory / "processed_protocol/transforms" / filename).write_text("corrupt transform")
    path.write_text(json.dumps(status))
    s = motion.analyze_subject("GBS001", directory)
    assert {key: s["processing"][key] for key in ["registered", "failed", "skipped", "missing"]} == {"registered": 39, "failed": 1, "skipped": 1, "missing": 1}
    assert s["processing"]["valid_call_transforms"] == 38
    rows = {(r["volume"], r["group"]): r for r in s["group_rows"]}
    for group in range(10, 14):
        assert rows[2, group]["translation_magnitude_mm"] is None
        assert rows[2, group]["group_displacement_mm"] is None
    assert rows[3, 0]["group_displacement_mm"] is None
    assert rows[3, 1]["group_displacement_mm"] is not None
    assert s["runtime_statistics_ms"]["count"] == 42  # measured runtimes survive transform loss


def test_absent_status_file_still_represents_expected_slots(tmp_path):
    directory = make_subject(tmp_path)
    (directory / "registration_status.json").unlink()
    s = motion.analyze_subject("GBS001", directory)
    assert s["processing"]["missing"] == 42
    assert s["processing"]["finalized_volumes"] == 0
    assert s["processing"]["registered"] == 0
    assert s["group_statistics"]["group_displacement_mm"]["count"] == 0
    assert len(s["group_rows"]) == 43


def test_preserve_subject_timing_and_missing_mrd_images(tmp_path):
    directory = make_subject(tmp_path, final_offset=.9075)
    meta = motion.acquisition_metadata(directory / "protocol_metadata.json")
    with h5py.File(directory / "measurement-test.hdf5", "a") as stream:
        del stream["image_57"]
    warnings = []
    times, count = motion.acquisition_times(directory, meta, warnings)
    assert times[0, 13][0] == .9075
    assert times[1, 13][0] == 1.9075
    assert times[1, 0][0] is None
    assert count == 223
    assert warnings
    (directory / "measurement-test.hdf5").unlink()
    warnings = []
    fallback, count = motion.acquisition_times(directory, meta, warnings)
    assert fallback[2, 13][0] == pytest.approx(2.9075)
    assert "derived_metadata" in fallback[2, 13][1]
    assert count is None and warnings


def test_discovery_duplicate_subjects_and_output_protection(tmp_path):
    directory = make_subject(tmp_path)
    (tmp_path / "GBS_six_subject_comparison").mkdir()
    (tmp_path / "savedData_fake_GBSabc").mkdir()
    assert motion.discover_subjects(tmp_path) == [("GBS001", directory.resolve())]
    with pytest.raises(motion.AnalysisError, match="Expected 6"):
        motion.run(tmp_path, tmp_path / "out")
    with pytest.raises(motion.AnalysisError, match="separate"):
        motion.run(tmp_path, directory / "analysis", expected_subjects=1)
    (tmp_path / "savedData_another_GBS001").mkdir()
    with pytest.raises(motion.AnalysisError, match="repeated subject"):
        motion.discover_subjects(tmp_path)


def test_aggregation_equal_subject_weight_and_unavailable_subject():
    subjects = [{"subject": name, "processing": {"expected_registrations": 14, "registered": 12},
                 "group_statistics": {"group_displacement_mm": {"count": count, "mean": mean}}}
                for name, count, mean in [("GBS001", 1000, 1), ("GBS002", 10, 3)]]
    subjects.append({"subject": "GBS003"})
    aggregate = motion.aggregate_subjects(subjects)
    assert aggregate["subject_count"] == 3
    assert aggregate["total_expected_registrations"] == 28
    assert set(aggregate["comparisons"]) == {"group_displacement_mm"}
    displacement = aggregate["comparisons"]["group_displacement_mm"]
    assert displacement["subjects_with_data"] == 2
    assert displacement["subject_mean_distribution"]["mean"] == 2
    assert displacement["ascending_subject_means"] == ["GBS001", "GBS002"]


def test_exports_figures_html_consistency_and_unavailable_subject(tmp_path):
    directory = make_subject(tmp_path)
    s = motion.analyze_subject("GBS001", directory)
    unavailable = {"subject": "GBS002", "input_directory": "unavailable", "warnings": ["Missing metadata"], "group_rows": []}
    output = tmp_path / "comparison"
    output.mkdir()
    motion.make_figures([s, unavailable], output, .3)
    summary = motion.write_reports([s, unavailable], output, 50, .3, ["Subject unavailable"])
    assert json.loads((output / "summary.json").read_text()) == summary
    groups = json.loads((output / "group_motion.json").read_text())
    with (output / "group_motion.csv").open() as stream:
        csv_groups = list(csv.DictReader(stream))
    assert len(csv_groups) == len(groups) == 43
    assert csv_groups[0]["group_displacement_mm"] == ""
    assert groups[0]["group_displacement_mm"] is None
    with (output / "motion_statistics.csv").open() as stream:
        stats = list(csv.DictReader(stream))
    displacement = next(r for r in stats if r["subject"] == "GBS001" and r["metric"] == "group_displacement_mm")
    assert float(displacement["mean"]) == summary["subjects"][0]["group_statistics"]["group_displacement_mm"]["mean"]
    assert {r["scope"] for r in stats} == {"group_statistics", "runtime"}
    assert_no_obsolete_metrics(summary)
    assert len(summary["subjects"]) == 2
    report = (output / "report.html").read_text()
    assert "GBS001" in report and "GBS002" in report
    assert "Consecutive slice-group displacement" in report
    assert_no_obsolete_metrics(report)
    assert "Missing metadata" in report and "radius_mm" in report and "t_delta" in report
    images = re.findall(r'src="data:image/png;base64,([^"]+)"', report)
    assert len(images) == 5
    for image, (name, _) in zip(images, motion.FIGURES):
        data = (output / name).read_bytes()
        assert data.startswith(b"\x89PNG")
        assert base64.b64decode(image) == data
    assert not re.search(r'(src|href)="https?://', report)


def test_html_affected_volumes_reuses_existing_six_subject_counts(tmp_path):
    expected = {"GBS024": 76, "GBS029": 4, "GBS061": 45, "GBS067": 93, "GBS069": 56, "GBS071": 227}
    subjects = [{"subject": subject, "group_event_affected_volumes": count, "warnings": [],
                 "group_statistics": {"group_displacement_mm": {
                     "count": 6706, "mean": .12, "median": .11, "p95": .23, "max": .43,
                     "exceedance_count": 79, "exceedance_pct": 1.178, "event_runs": 999}}}
                for subject, count in expected.items()]
    summary = {"subjects": subjects, "generated_at_utc": "2026-11-16T12:00:00+00:00",
               "methods": motion.methods(50, .3), "warnings": [], "cross_subject": motion.aggregate_subjects(subjects)}
    (tmp_path / "summary.json").write_text(json.dumps(summary))
    for name, _ in motion.FIGURES:
        (tmp_path / name).write_bytes(b"existing figure bytes")
    before = snapshot(tmp_path)
    motion.write_html_report(summary, tmp_path)
    assert {name: digest for name, digest in snapshot(tmp_path).items() if name != "report.html"} == before
    assert json.loads((tmp_path / "summary.json").read_text()) == summary
    report = (tmp_path / "report.html").read_text()
    headers, rows = displacement_summary_table(report)
    assert headers == ["Subject", "Valid transitions", "Mean (mm)", "Median (mm)", "P95 (mm)", "Max (mm)",
                       ">0.3 mm count", "Exceedance (%)", "Affected volumes (>0.3 mm)"]
    assert "Event runs" not in headers
    assert {row[0]: int(row[-1]) for row in rows} == expected
    assert [row[0] for row in rows] == list(expected)
    assert all(row[1:-1] == ["6706", "0.12", "0.11", "0.23", "0.43", "79", "1.178"] for row in rows)
    assert "exceeding 0.30 mm. Each volume is counted once" in report
    assert "Cross-volume transitions are included and attributed to the volume containing the current slice group." in report


def test_duplicate_records_and_conflicting_transforms_refused(tmp_path):
    with pytest.raises(motion.AnalysisError, match="Duplicate"):
        motion.unique_by_key([{"v": 1, "g": 2}] * 2, ("v", "g"), "test")
    directory = make_subject(tmp_path)
    source = next((directory / "processed_protocol/transforms").glob("*.tfm"))
    sitk.WriteTransform(euler(translation=(99, 0, 0)), str(directory / source.name))
    with pytest.raises(motion.AnalysisError, match="Conflicting duplicate"):
        motion.load_transform(directory, source.name, [])


def test_obsolete_output_cleanup_is_specific_and_repeatable(tmp_path):
    obsolete = ["volume_motion.csv", "volume_motion.json", "figure_2_framewise_displacement.png"]
    retained = ["group_motion.csv", "validation.json", "validation_tests.xml", "notes.txt"]
    for name in obsolete + retained:
        (tmp_path / name).write_text(name)
    motion.remove_obsolete_outputs(tmp_path)
    motion.remove_obsolete_outputs(tmp_path)
    assert all(not (tmp_path / name).exists() for name in obsolete)
    assert all((tmp_path / name).read_text() == name for name in retained)


def test_displacement_figures_use_every_group_and_include_requested_statistics(tmp_path, monkeypatch):
    directory = make_subject(tmp_path)
    s = motion.analyze_subject("GBS001", directory)
    subjects = [{**s, "subject": f"GBS{i:03d}"} for i in range(1, 7)]
    # Include a large valid group value to exercise the retained symlog behavior.
    for subject in subjects:
        subject["group_rows"] = [dict(row) for row in s["group_rows"]]
        subject["group_rows"][5]["group_displacement_mm"] = 2.5
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes

    figures, box_data = {}, []
    def capture_figure(fig, path, **kwargs):
        figures[Path(path).name] = fig
    original_boxplot = Axes.boxplot
    def capture_boxplot(ax, data, *args, **kwargs):
        box_data.append(data)
        return original_boxplot(ax, data, *args, **kwargs)
    monkeypatch.setattr(Figure, "savefig", capture_figure)
    monkeypatch.setattr(Axes, "boxplot", capture_boxplot)
    motion.make_figures(subjects, tmp_path, .3)
    fig2 = figures[motion.FIGURES[1][0]]
    assert len(fig2.axes) == 6
    for ax, subject in zip(fig2.axes, subjects):
        rows = subject["group_rows"]
        np.testing.assert_equal(ax.lines[0].get_xdata(), [r["acquisition_time_s"] for r in rows])
        np.testing.assert_equal(ax.lines[0].get_ydata(), [np.nan if r["group_displacement_mm"] is None else r["group_displacement_mm"] for r in rows])
        assert ax.get_yscale() == "symlog"
        assert ax.get_ylim() == fig2.axes[0].get_ylim()
        assert "Consecutive slice-group" in ax.get_ylabel()
        assert list(ax.lines[1].get_ydata()) == [.3, .3]
    assert "Acquisition time (s)" == fig2.axes[-1].get_xlabel()
    fig3 = figures[motion.FIGURES[2][0]]
    assert len(box_data) == 2  # group distributions, then registration runtimes
    assert len(box_data[0]) == 6
    for values, subject in zip(box_data[0], subjects):
        np.testing.assert_equal(values, [r["group_displacement_mm"] for r in subject["group_rows"][1:]])
    cells = list(fig3.axes[1].tables[0].get_celld().values())
    text = " ".join(cell.get_text().get_text() for cell in cells)
    for label in ["Mean (mm)", "Median (mm)", "P95 (mm)", "Maximum (mm)", "> 0.3 mm (%)"]:
        assert label in text
    for subject in subjects:
        st = motion.metric_statistics([r["group_displacement_mm"] for r in subject["group_rows"]], .3)
        for metric in ["mean", "median", "p95", "max", "exceedance_pct"]:
            assert f"{st[metric]:.3f}" in text
    assert "Consecutive slice-group displacement" in fig3.axes[0].get_ylabel()
    for fig in [fig2, fig3]:
        texts = [fig._suptitle.get_text(), *[ax.get_title() for ax in fig.axes]]
        assert_no_obsolete_metrics(" ".join(texts))


def test_generated_six_subject_artifacts_if_available():
    """Audit actual exports independently; skip when replay data aren't installed."""
    data_root = motion.HERE.parents[1] / "data"
    output = data_root / "GBS_six_subject_comparison"
    if not (output / "summary.json").exists():
        pytest.skip("Generate the six-subject report to validate its real-data exports")
    summary = json.loads((output / "summary.json").read_text())
    groups = json.loads((output / "group_motion.json").read_text())
    expected_subjects = {"GBS024", "GBS029", "GBS061", "GBS067", "GBS069", "GBS071"}
    assert {s["subject"] for s in summary["subjects"]} == expected_subjects
    assert {r["subject"] for r in groups} == expected_subjects
    assert_no_obsolete_metrics(summary)
    assert summary["cross_subject"]["subject_count"] == 6
    assert summary["cross_subject"]["total_expected_registrations"] == 6 * 6706
    for s in summary["subjects"]:
        subject = s["subject"]
        gr = [r for r in groups if r["subject"] == subject]
        assert len(gr) == 6707
        assert gr[0]["is_reference"] is True and gr[0]["group_displacement_mm"] is None
        assert all({r["group"] for r in gr if r["volume"] == v} == set(range(14)) for v in range(1, 480))
        assert s["warnings"] == []
        status = json.loads((Path(s["input_directory"]) / "registration_status.json").read_text())
        outcomes = [g["status"] for v, record in status["volumes"].items() if int(v) != 0 for g in record["groups"].values()]
        for key in ("registered", "failed", "skipped"):
            assert s["processing"][key] == outcomes.count(key)
        assert s["processing"]["missing"] == 0
        assert s["processing"]["valid_call_transforms"] == 6706
        assert s["recorded_group_displacement_crosscheck"]["paired_count"] == 6706
        assert s["recorded_group_displacement_crosscheck"]["max_abs_error_mm"] < 1e-6
        for rows, metric, stats in [(gr, "group_displacement_mm", s["group_statistics"]["group_displacement_mm"]),
                                   (gr, "runtime_ms", s["runtime_statistics_ms"])]:
            values = [r[metric] for r in rows if r[metric] is not None]
            assert stats["count"] == len(values)
            assert stats["mean"] == pytest.approx(np.mean(values))
            assert stats["median"] == pytest.approx(np.median(values))
            assert stats["p95"] == pytest.approx(np.percentile(values, 95))
            assert stats["max"] == max(values)
            if "exceedance_count" in stats:
                assert stats["exceedance_count"] == sum(v > .3 for v in values)
                assert stats["exceedance_pct"] == pytest.approx(100 * sum(v > .3 for v in values) / len(values))
        assert all(r["group_displacement_mm"] is not None for r in gr[1:])
        assert sum(r["group"] == 0 and r["volume"] > 1 for r in gr) == 478
        assert gr[-1]["acquisition_time_s"] == (479.9075 if subject in {"GBS024", "GBS071"} else 479.91)
    with (output / "group_motion.csv").open() as stream:
        assert len(list(csv.DictReader(stream))) == len(groups)
    with (output / "subject_summary.csv").open() as stream:
        assert {r["subject"] for r in csv.DictReader(stream)} == expected_subjects
    report = (output / "report.html").read_text()
    headers, rows = displacement_summary_table(report)
    assert headers[-1] == "Affected volumes (>0.3 mm)"
    assert "Event runs" not in headers
    expected_affected = {"GBS024": 76, "GBS029": 4, "GBS061": 45, "GBS067": 93, "GBS069": 56, "GBS071": 227}
    assert {row[0]: int(row[-1]) for row in rows} == expected_affected
    assert {s["subject"]: s["group_event_affected_volumes"] for s in summary["subjects"]} == expected_affected
    assert_no_obsolete_metrics(report)
    assert_no_obsolete_metrics((output / "motion_statistics.csv").read_text())
    for name in ["volume_motion.csv", "volume_motion.json", "figure_2_framewise_displacement.png"]:
        assert not (output / name).exists()
    assert all(subject in report for subject in expected_subjects)
    images = re.findall(r'src="data:image/png;base64,([^"]+)"', report)
    assert len(images) == 5
    for encoded, (filename, _) in zip(images, motion.FIGURES):
        assert base64.b64decode(encoded) == (output / filename).read_bytes()

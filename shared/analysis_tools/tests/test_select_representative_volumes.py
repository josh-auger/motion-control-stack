"""Focused checks for read-only, reproducible representative-volume selection."""
from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from shared.analysis_tools import select_representative_volumes as selection


def subject(volumes=3):
    rows = [{"is_reference": True, "volume": 0, "group": 14}]
    for volume in range(1, volumes + 1):
        for group in range(14):
            rows.append({
                "subject": "GBS001", "volume": volume, "group": group, "is_reference": False,
                "group_displacement_mm": .05 if group else 2., "status": "registered",
                "transform_valid": True, "registration_finalized": True, "pointer_valid": True,
                "translation_magnitude_mm": float(group), "rotation_angle_deg": group / 10,
            })
    return {"subject": "GBS001", "metadata": {"volumes": volumes + 1}, "input_directory": "/source/GBS001", "group_rows": rows}


@pytest.mark.parametrize("value,expected", [(.1, "minimal"), (.10001, "unclassified"), (.14999, "unclassified"),
                                          (.15, "low_moderate"), (.3, "low_moderate"), (.30001, "high"), (None, "incomplete")])
def test_exact_category_boundaries(value, expected):
    assert selection.category(value) == expected


def test_group_zero_is_separate_and_each_volume_has_thirteen_intra_transitions():
    data = subject()
    data["group_rows"][1 + 14 + 3]["group_displacement_mm"] = .4  # v2/group3
    before = copy.deepcopy(data)
    rows = selection.volume_characteristics(data)
    assert data == before
    assert len(rows) == 3
    assert rows[0]["category"] == "minimal"  # 2 mm at group 0 does not classify the volume high.
    assert rows[0]["cross_volume_displacement_mm"] == 2
    assert rows[0]["exceedance_count"] == 0
    assert rows[0]["maximum_reference_relative_translation_mm"] == 13
    assert rows[0]["maximum_reference_relative_rotation_deg"] == 1.3
    assert rows[1]["category"] == "high"
    assert rows[1]["exceedance_group_indices"] == [3]
    assert rows[1]["valid_intra_volume_transitions"] == 13
    assert rows[1]["mean_intra_volume_displacement_mm"] == pytest.approx((12 * .05 + .4) / 13)
    assert rows[1]["p95_intra_volume_displacement_mm"] == pytest.approx(np.percentile([.05] * 12 + [.4], 95))
    assert 0 not in rows[1]["exceedance_group_indices"]


def test_missing_observations_prevent_selection_without_zero_filling():
    data = subject(1)
    data["group_rows"][4]["group_displacement_mm"] = None
    data["group_rows"][4]["transform_valid"] = False
    row = selection.volume_characteristics(data)[0]
    assert row["valid_intra_volume_transitions"] == 12
    assert row["category"] == "incomplete"
    assert not row["complete_registration_records"]
    assert row["mean_intra_volume_displacement_mm"] == pytest.approx(.05)


def test_explicit_metadata_order_and_duplicates():
    data = subject(1)
    baseline = selection.volume_characteristics(data)
    data["group_rows"][1:] = list(reversed(data["group_rows"][1:]))
    assert selection.volume_characteristics(data) == baseline
    data["group_rows"].append(data["group_rows"][-1])
    with pytest.raises(selection.AnalysisError, match="Duplicate"):
        selection.volume_characteristics(data)


def candidate(subject_id, volume, maximum):
    cat = selection.category(maximum)
    return {"subject": subject_id, "volume_index": volume, "category": cat, "complete_registration_records": True,
            "cross_volume_displacement_mm": .05, "adjacent_large_event_volumes": [],
            "mean_intra_volume_displacement_mm": maximum / 2, "maximum_intra_volume_displacement_mm": maximum,
            "p95_intra_volume_displacement_mm": maximum * .9, "exceedance_count": 2 if cat == "high" else 0,
            "moderately_elevated_transition_count": 5 if cat == "low_moderate" else 0}


def test_deterministic_primary_diversity_and_category_counts():
    candidates = [candidate(s, v, maximum) for s in selection.SUBJECTS for v, maximum in [(10, .08), (20, .2), (30, .5)]]
    candidates.append(candidate("GBS071", 31, 1.8))
    original = copy.deepcopy(candidates)
    first = selection.select_primary(candidates)
    second = selection.select_primary(list(reversed(candidates)))
    assert first == second
    assert candidates == original
    assert len({r["subject"] for r in first}) == 6
    assert {cat: sum(r["category"] == cat for r in first) for cat in selection.CATEGORIES} == {
        "minimal": 2, "low_moderate": 2, "high": 2}
    assert any(r["subject"] == "GBS071" and r["volume_index"] == 31 and r["category"] == "high" for r in first)


def test_image_source_missing_or_ambiguous_is_unresolved(tmp_path):
    replay = tmp_path / "savedData_test_GBS061"
    replay.mkdir()
    (replay / "scan_metadata.json").write_text(json.dumps({"userParameters": {"userParameterString": [{"name": "SourceNIfTI", "value": "source.nii.gz"}]}}))
    subjects = [{"subject": "GBS061", "input_directory": str(replay)}]
    assert selection.find_nifti_sources(tmp_path, subjects) == {"GBS061": None}
    a = tmp_path / "scan_GBS061_first"
    a.mkdir()
    path = a / "source.nii.gz"
    path.write_text("placeholder")
    assert selection.find_nifti_sources(tmp_path, subjects) == {"GBS061": path.resolve()}
    b = tmp_path / "scan_GBS061_second"
    b.mkdir()
    (b / "source.nii.gz").write_text("placeholder")
    assert selection.find_nifti_sources(tmp_path, subjects) == {"GBS061": None}


def test_existing_comparison_directory_cannot_be_an_output(tmp_path):
    with pytest.raises(selection.AnalysisError, match="new subdirectory"):
        selection.run(tmp_path, tmp_path, None)


def test_global_maximum_includes_boundaries_and_strict_intra_excludes_them():
    records = selection.volume_characteristics(subject(2))
    # Group 0 is 2 mm; all intra-volume values are 0.05 mm.
    records[1]["_trajectory"][3]["group_displacement_mm"] = .8
    first, row = selection.largest_transition(records)
    assert first["volume_index"] == 1 and row["group"] == 0
    assert row["group_displacement_mm"] == 2
    intra, row = selection.largest_transition(records, intra_only=True)
    assert intra["volume_index"] == 2 and row["group"] == 3
    assert row["group_displacement_mm"] == .8
    records[0]["_trajectory"][0]["group_displacement_mm"] = None
    records[1]["_trajectory"][0]["group_displacement_mm"] = float("nan")
    assert selection.largest_transition(records)[1]["group"] == 3
    with pytest.raises(selection.AnalysisError, match="No valid"):
        selection.largest_transition([])


def test_maximum_role_reuses_image_mapping_without_changing_original_selection():
    record = {"subject": "GBS071", "volume_index": 162}
    original = {**record, "category": "high", "nifti_path": "/verified/source.nii.gz", "nifti_volume_index": 162,
                "selection_rationale": "Original rationale", "reference_relative_trajectory": [{"group": 13}]}
    before = copy.deepcopy(original)
    extra = selection.maximum_selection(record, {"group": 13, "group_displacement_mm": 1.8796},
                                        selection.INTRA_MAXIMUM_CATEGORY, [original], None)
    assert original == before
    assert extra["category"] == selection.INTRA_MAXIMUM_CATEGORY
    assert extra["maximum_transition_type"] == "intra_volume"
    assert extra["maximum_transition_slice_group_index"] == 13
    assert extra["existing_selection_reference"] == {**record, "category": "high"}
    assert extra["nifti_path"] == original["nifti_path"]
    extra["reference_relative_trajectory"][0]["group"] = 0
    assert original == before  # The additional role must not alias mutable original data.


def test_maximum_count_includes_boundary_uses_strict_threshold_and_reports_ties():
    records = selection.volume_characteristics(subject(2))
    # Fourteen valid measurements, one boundary and two intra-volume exceedances.
    for record in records:
        for row in record["_trajectory"]:
            row["group_displacement_mm"] = .3
        for group in (0, 2, 13):
            record["_trajectory"][group]["group_displacement_mm"] = .4
    before = copy.deepcopy(records)
    winner, stats, ties = selection.most_exceedances(list(reversed(records)))
    assert winner["volume_index"] == 1  # Stable subject/volume tie-break, independent of input order.
    assert stats["all_group_exceedance_count"] == 3
    assert stats["all_group_exceedance_indices"] == [0, 2, 13]
    assert stats["all_group_exceedance_pct"] == pytest.approx(100 * 3 / 14)
    assert stats["group_0_exceedance_count"] == 1
    assert stats["intra_volume_exceedance_count"] == 2
    assert stats["mean_all_group_displacement_mm"] == pytest.approx(np.mean([.4] * 3 + [.3] * 11))
    assert stats["p95_all_group_displacement_mm"] == pytest.approx(.4)
    assert stats["maximum_all_group_displacement_mm"] == .4
    assert [(r["subject"], r["volume_index"]) for r in ties] == [("GBS001", 2)]
    assert records == before
    with pytest.raises(selection.AnalysisError, match="No non-reference"):
        selection.most_exceedances([])


def test_maximum_count_does_not_turn_missing_transitions_into_zero():
    record = selection.volume_characteristics(subject(1))[0]
    record["_trajectory"][1]["group_displacement_mm"] = None
    stats = selection.all_group_exceedance_statistics(record)
    assert stats["valid_all_group_transition_count"] == 13
    assert stats["all_group_exceedance_pct"] == pytest.approx(100 / 13)
    record["_trajectory"][2]["group"] = 1
    with pytest.raises(selection.AnalysisError, match="unique groups"):
        selection.all_group_exceedance_statistics(record)


def test_actual_selection_manifest_if_available():
    root = selection.HERE.parents[1] / "data/GBS_six_subject_comparison"
    path = root / "representative_volumes/representative_volume_selection.json"
    if not path.exists():
        pytest.skip("Generate the representative selection manifest first")
    document = json.loads(path.read_text())
    selected = document["selections"]
    assert len(selected) == 12
    primary = [r for r in selected if r["category"] in selection.CATEGORIES]
    assert len(primary) == 6 and {r["subject"] for r in primary} == set(selection.SUBJECTS)
    assert len({(r["subject"], r["volume_index"]) for r in selected}) == 11
    assert len(document["all_non_reference_volume_characteristics"]) == 2874
    assert document["eligibility_totals"] == {"minimal": 30, "low_moderate": 1891, "high": 432, "unclassified": 521}
    saved_rows = json.loads((root / "group_motion.json").read_text())
    global_maximum = next(r for r in selected if r["category"] == selection.GLOBAL_MAXIMUM_CATEGORY)
    assert (global_maximum["subject"], global_maximum["volume_index"], global_maximum["maximum_transition_slice_group_index"]) == ("GBS071", 163, 0)
    assert global_maximum["maximum_transition_type"] == "cross_volume"
    assert global_maximum["maximum_transition_displacement_mm"] == max(r["group_displacement_mm"] for r in saved_rows if r["group_displacement_mm"] is not None)
    assert global_maximum["maximum_transition_displacement_mm"] == pytest.approx(2.619435689789227)
    assert global_maximum["exceedance_count"] == 3  # Its large boundary does not enter the intra-volume count.
    intra_maximum = next(r for r in selected if r["category"] == selection.INTRA_MAXIMUM_CATEGORY)
    assert (intra_maximum["subject"], intra_maximum["volume_index"], intra_maximum["maximum_transition_slice_group_index"]) == ("GBS071", 162, 10)
    assert intra_maximum["maximum_transition_type"] == "intra_volume"
    assert intra_maximum["maximum_transition_displacement_mm"] == max(r["group_displacement_mm"] for r in saved_rows if not r["is_reference"] and r["group"] > 0)
    assert intra_maximum["maximum_transition_displacement_mm"] == pytest.approx(1.8796206142685548)
    assert document["validation"]["global_maximum_exactly_matches_comparison"]
    assert document["validation"]["original_selections_unchanged"]
    count_peak = next(r for r in selected if r["category"] == selection.MAXIMUM_COUNT_CATEGORY)
    assert (count_peak["subject"], count_peak["volume_index"]) == ("GBS069", 218)
    assert count_peak["valid_all_group_transition_count"] == 14
    assert count_peak["all_group_exceedance_count"] == 7
    assert count_peak["all_group_exceedance_pct"] == 50
    assert count_peak["all_group_exceedance_indices"] == [0, 2, 5, 6, 10, 11, 13]
    assert count_peak["group_0_exceedance_count"] == 1
    assert count_peak["intra_volume_exceedance_count"] == count_peak["exceedance_count"] == 6
    assert count_peak["mean_all_group_displacement_mm"] == pytest.approx(.2968347809559768)
    assert count_peak["p95_all_group_displacement_mm"] == pytest.approx(.4133183288469402)
    assert count_peak["maximum_all_group_displacement_mm"] == pytest.approx(.4174655400416063)
    assert document["maximum_exceedance_count"]["other_tied_volumes"] == []
    assert document["maximum_exceedance_count"]["gbs071_volume_162_comparison"]["all_group_exceedance_count"] == 6
    markdown = (path.parent / "representative_volume_selection.md").read_text()
    count_section = markdown.split("## Maximum slice-group exceedance count\n", 1)[1].split("\n## Verified image indexing", 1)[0]
    table_rows = [line.split("|")[1:-1] for line in count_section.splitlines() if line.startswith("| GBS")]
    assert [[cell.strip() for cell in row[:4]] for row in table_rows] == [
        ["GBS069", "218", "7 / 14", "50.000000%"],
        ["GBS071", "162", "6 / 14", "42.857143%"],
    ]
    assert [[cell.strip() for cell in row[-2:]] for row in table_rows] == [["1", "6"], ["0", "6"]]
    assert document["validation"]["exceedance_count_search_volumes"] == 2874
    from collections import defaultdict
    observations = defaultdict(list)
    for row in saved_rows:
        if not row["is_reference"]:
            observations[row["subject"], row["volume"]].append(row["group_displacement_mm"])
    counts = {key: sum(value > .3 for value in values) for key, values in observations.items()}
    assert max(counts.values()) == 7
    assert [key for key, value in counts.items() if value == 7] == [("GBS069", 218)]
    assert [(r["subject"], r["volume_index"]) for r in selected[:9]] == [
        ("GBS061", 17), ("GBS069", 29), ("GBS024", 340), ("GBS029", 246), ("GBS071", 162), ("GBS067", 397),
        ("GBS061", 463), ("GBS067", 413), ("GBS071", 357)]
    for r in selected:
        group_rows = [g for g in saved_rows if g["subject"] == r["subject"] and g["volume"] == r["volume_index"]]
        assert [g["group"] for g in group_rows] == list(range(14))
        values = [g["group_displacement_mm"] for g in group_rows[1:]]
        assert r["mean_intra_volume_displacement_mm"] == pytest.approx(np.mean(values))
        assert r["maximum_intra_volume_displacement_mm"] == max(values)
        assert r["p95_intra_volume_displacement_mm"] == pytest.approx(np.percentile(values, 95))
        assert r["exceedance_count"] == sum(v > .3 for v in values)
        assert r["exceedance_group_indices"] == [g["group"] for g in group_rows[1:] if g["group_displacement_mm"] > .3]
        assert r["cross_volume_displacement_mm"] == group_rows[0]["group_displacement_mm"]
        assert r["valid_intra_volume_transitions"] == 13
        assert r["complete_registration_records"]
        assert r["image_mapping_verified"] and r["verified_slice_count"] == 56
        assert r["nifti_mrd_maximum_voxel_error"] == 0
        assert r["nifti_volume_index"] == r["volume_index"]
        assert {s["mrd_repetition"] for s in r["slice_mapping"]} == {r["volume_index"]}
        assert {s["spatial_slice_index"] for s in r["slice_mapping"]} == set(range(56))
        assert all(Path(s["nrrd_header"]).is_file() for s in r["slice_mapping"])
        assert len(r["reference_relative_trajectory"]) == 14
    for cat in selection.CATEGORIES:
        assert sum(r["category"] == cat for r in primary) == 2
    assert all(r["net_group_2_to_4_change_mm"] < .5 * (r["group_3_displacement_mm"] + r["group_4_displacement_mm"])
               for r in [s["diagnostic_window"] for s in selected if s["category"] == "diagnostic"])

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk

from shared.analysis_tools.compare_motion_control_images import (
    LoadedImage,
    _create_diagnostic_figure,
    _format_intensity_fit_annotation,
    _sitk_geometry,
    _symmetric_normalized_difference,
    compare_arrays,
    compare_directories,
    compare_physical_space,
    discover_images,
    match_images,
)


SCRIPT = Path(__file__).parents[1] / "compare_motion_control_images.py"


def loaded(image: sitk.Image) -> LoadedImage:
    array = sitk.GetArrayFromImage(image)
    return LoadedImage(array, _sitk_geometry(image, array.shape), image)


def sitk_image(array, *, spacing=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0), direction=None):
    image = sitk.GetImageFromArray(np.asarray(array))
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    if direction is not None:
        image.SetDirection(direction)
    return image


def write_nhdr(
    directory: Path,
    name: str,
    array: np.ndarray,
    *,
    origin=(0.0, 0.0, 0.0),
    directions=None,
) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    array = np.asarray(array)
    raw_path = directory / f"{name}.raw"
    array.tofile(raw_path)
    path = directory / f"{name}.nhdr"
    if directions is None:
        directions = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    sizes = tuple(reversed(array.shape))
    vector_text = []
    for item in directions:
        vector_text.append("none" if item is None else f"({','.join(map(str, item))})")
    dtype_names = {
        np.dtype("float32"): "float",
        np.dtype("float64"): "double",
        np.dtype("int16"): "short",
        np.dtype("uint16"): "ushort",
        np.dtype("int32"): "int",
    }
    path.write_text(
        "\n".join(
            [
                "NRRD0004",
                f"type: {dtype_names[array.dtype]}",
                f"dimension: {array.ndim}",
                f"sizes: {' '.join(map(str, sizes))}",
                "space: left-posterior-superior",
                f"space origin: ({','.join(map(str, origin))})",
                f"space directions: {' '.join(vector_text)}",
                "encoding: raw",
                "endian: little",
                f"data file: {raw_path.name}",
                "",
            ]
        ),
        encoding="ascii",
    )
    return path


def image_name(prefix: str, volume: int, slice_number: int = 0) -> str:
    return f"{prefix}_volume_{volume:04d}_slice_{slice_number:04d}"


def run_comparison(tmp_path: Path, arrays_a, arrays_b, **kwargs):
    dir_a, dir_b = tmp_path / "a", tmp_path / "b"
    for index, array in enumerate(arrays_a):
        write_nhdr(dir_a, image_name("native", index), array, **kwargs)
    for index, array in enumerate(arrays_b):
        write_nhdr(dir_b, image_name("converted", index), array, **kwargs)
    return compare_directories(dir_a, dir_b, tmp_path / "reports")


def test_identical_image_arrays(tmp_path):
    array = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    summary, rows, exit_code = run_comparison(tmp_path, [array], [array.copy()])
    assert exit_code == 0
    assert summary["overall_result"] == "PASS_EXACT"
    assert rows[0]["unequal_voxels"] == 0
    assert rows[0]["byte_equal"] is True


def test_single_voxel_difference(tmp_path):
    a = np.zeros((1, 3, 4), dtype=np.float32)
    b = a.copy()
    b[0, 1, 2] = 5
    summary, rows, exit_code = run_comparison(tmp_path, [a], [b])
    assert exit_code == 1
    assert rows[0]["unequal_voxels"] == 1
    assert rows[0]["max_absolute_difference"] == 5
    assert "FAIL_VOXELS" in summary["result_categories"]


def test_small_floating_point_difference_is_not_tolerated():
    a = np.array([1.0], dtype=np.float64)
    b = np.nextafter(a, np.inf)
    result = compare_arrays(a, b)
    assert result["exact_voxel_equal"] is False
    assert result["unequal_voxels"] == 1
    assert result["max_absolute_difference"] > 0


def test_different_array_shapes(tmp_path):
    a = np.zeros((1, 2, 3), dtype=np.float32)
    b = np.zeros((1, 2, 4), dtype=np.float32)
    _, rows, exit_code = run_comparison(tmp_path, [a], [b])
    assert exit_code == 1
    assert rows[0]["comparison_result"] == "FAIL_SHAPE"
    assert rows[0]["unequal_voxels"] is None


def test_different_data_types_fail_exact_equality():
    result = compare_arrays(np.arange(4, dtype=np.int16), np.arange(4, dtype=np.float32))
    assert result["dtype_equal"] is False
    assert result["unequal_voxels"] == 0
    assert result["exact_voxel_equal"] is False


def test_identical_voxels_with_different_geometry(tmp_path):
    array = np.arange(6, dtype=np.float32).reshape(1, 2, 3)
    a, b = tmp_path / "a", tmp_path / "b"
    write_nhdr(a, image_name("a", 0), array, origin=(0, 0, 0))
    write_nhdr(b, image_name("b", 0), array, origin=(10, 0, 0))
    summary, rows, exit_code = compare_directories(a, b, tmp_path / "reports")
    assert exit_code == 1
    assert summary["overall_result"] == "PASS_VOXELS_ONLY"
    assert set(summary["result_categories"]) == {"PASS_VOXELS_ONLY", "FAIL_GEOMETRY"}
    assert summary["spatial_overall_status"] in {"SPATIAL_FAIL_GRID", "SPATIAL_FAIL_COVERAGE"}
    assert summary["pixel_overall_status"] == "PIXEL_PASS_EXACT"
    assert rows[0]["geometry_status"] == "DIFFERENT"
    _, _, pixels_exit = compare_directories(
        a, b, tmp_path / "pixel-reports", validation_target="pixels"
    )
    assert pixels_exit == 0


def test_missing_and_extra_images_are_coverage_failure(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    array = np.zeros((1, 2, 2), dtype=np.float32)
    write_nhdr(a, image_name("a", 0), array)
    write_nhdr(a, image_name("a", 1), array)
    write_nhdr(b, image_name("b", 0), array)
    write_nhdr(b, image_name("b", 2), array)
    summary, _, exit_code = compare_directories(a, b, tmp_path / "reports")
    assert exit_code == 1
    assert summary["coverage"]["missing_from_a"] == ["volume:0002/slice:0000"]
    assert summary["coverage"]["missing_from_b"] == ["volume:0001/slice:0000"]
    assert "FAIL_COVERAGE" in summary["result_categories"]


def test_duplicate_image_identifiers_are_ambiguous(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    array = np.zeros((1, 2, 2), dtype=np.float32)
    write_nhdr(a, image_name("first", 3, 4), array)
    write_nhdr(a, image_name("second", 3, 4), array)
    write_nhdr(b, image_name("only", 3, 4), array)
    images_a, _ = discover_images(a)
    images_b, _ = discover_images(b)
    matches = match_images(images_a, images_b)
    assert matches["pairs"] == []
    assert matches["ambiguous_identifiers"] == ["volume:0003/slice:0004"]


def test_different_filenames_match_by_image_identity(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    array = np.ones((1, 2, 2), dtype=np.float32)
    write_nhdr(a, image_name("scanner_protocol", 12, 7), array)
    write_nhdr(b, image_name("renamed_converted", 12, 7), array)
    summary, rows, _ = compare_directories(a, b, tmp_path / "reports")
    assert summary["image_counts"]["matched_pairs"] == 1
    assert rows[0]["image_identifier"] == "volume:0012/slice:0007"


def test_nan_and_inf_handling():
    a = np.array([np.nan, np.inf, -np.inf, 1.0], dtype=np.float32)
    b = a.copy()
    equal = compare_arrays(a, b)
    assert equal["exact_voxel_equal"] is True
    assert equal["nonfinite_a"] == {"nan": 1, "positive_inf": 1, "negative_inf": 1}
    b[1] = -np.inf
    unequal = compare_arrays(a, b)
    assert unequal["unequal_voxels"] == 1
    assert unequal["finite_pairs_for_statistics"] == 1


def test_multiple_image_files_and_reports(tmp_path):
    arrays = [np.full((1, 2, 3), index, dtype=np.float32) for index in range(3)]
    summary, rows, exit_code = run_comparison(tmp_path, arrays, arrays)
    reports = tmp_path / "reports"
    assert exit_code == 0
    assert len(rows) == 3
    loaded = json.loads((reports / "comparison_summary.json").read_text())
    assert loaded["image_counts"]["matched_pairs"] == 3
    with (reports / "image_comparison.csv").open(newline="") as stream:
        assert len(list(csv.DictReader(stream))) == 3
    assert "PASS_EXACT" in (reports / "comparison_report.md").read_text()


def test_multidimensional_4d_image_arrays(tmp_path):
    array = np.arange(2 * 3 * 4 * 5, dtype=np.float32).reshape(2, 3, 4, 5)
    directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1), None]
    dir_a, dir_b = tmp_path / "a", tmp_path / "b"
    write_nhdr(dir_a, image_name("native", 0), array, directions=directions)
    write_nhdr(dir_b, image_name("converted", 0), array, directions=directions)
    summary, rows, exit_code = compare_directories(
        dir_a, dir_b, tmp_path / "reports", validation_target="pixels"
    )
    assert exit_code == 0
    assert rows[0]["shape_a"] == [2, 3, 4, 5]
    assert summary["aggregate"]["total_voxels"] == array.size


@pytest.mark.parametrize("different, expected", [(False, 0), (True, 1)])
def test_cli_process_exit_codes_and_outputs(tmp_path, different, expected):
    a, b = tmp_path / "a", tmp_path / "b"
    first = np.zeros((1, 2, 2), dtype=np.float32)
    second = first.copy()
    if different:
        second[0, 0, 0] = 1
    write_nhdr(a, image_name("a", 0), first)
    write_nhdr(b, image_name("b", 0), second)
    output = tmp_path / "cli-reports"
    process = subprocess.run(
        [
            sys.executable, str(SCRIPT), "--dir-a", str(a), "--dir-b", str(b),
            "--output-dir", str(output), "--label-a", "A", "--label-b", "B",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert process.returncode == expected, process.stderr
    assert (output / "comparison_summary.json").is_file()
    assert (output / "image_comparison.csv").is_file()
    assert (output / "comparison_report.md").is_file()


def test_image_pattern_and_derived_exclusion(tmp_path):
    directory = tmp_path / "images"
    array = np.zeros((1, 2, 2), dtype=np.float32)
    write_nhdr(directory, image_name("keep", 0), array)
    write_nhdr(directory, f"protocol_volume_0000_upsampled", array)
    records, excluded = discover_images(directory, "*.nhdr")
    assert list(records) == ["volume:0000/slice:0000"]
    assert any("derived image marker" in item["reason"] for item in excluded)


def test_simpleitk_identical_geometry_passes_physical_grid():
    image = sitk_image(np.arange(24, dtype=np.float32).reshape(2, 3, 4), spacing=(1.5, 2.0, 3.0))
    result = compare_physical_space(loaded(image), loaded(sitk.Image(image)))
    assert result["lossless_reorientation_feasible"] is True
    assert result["alignment_method"] == "lossless_discrete_axis_reorientation"
    assert result["spatial_status"] == "SPATIAL_PASS_EXACT"
    assert result["pixel_status"] == "PIXEL_PASS_EXACT"
    assert result["aligned_intensity_metrics"]["unequal_voxels"] == 0


def test_reversed_voxel_axis_with_compensating_geometry_is_lossless():
    array = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)
    native = sitk_image(array)
    converted = sitk_image(
        np.flip(array, axis=1),
        origin=(0.0, 2.0, 0.0),
        direction=(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0),
    )
    result = compare_physical_space(loaded(native), loaded(converted))
    assert result["axis_relation"]["reflected_source_axes"] == [1]
    assert result["lossless_reorientation_feasible"] is True
    assert result["axis_relation"]["corresponding_sample_max_discrepancy_mm"] == pytest.approx(0.0)
    assert result["aligned_intensity_metrics"]["unequal_voxels"] == 0
    # B index (x, y, z) corresponds to A index (x, 2-y, z).
    assert converted.TransformIndexToPhysicalPoint((1, 0, 1)) == native.TransformIndexToPhysicalPoint((1, 2, 1))


def test_axis_permutation_with_compensating_geometry_is_lossless():
    native = sitk_image(np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4))
    converted = sitk.PermuteAxes(native, [1, 0, 2])
    result = compare_physical_space(loaded(native), loaded(converted))
    assert result["axis_relation"]["source_axis_to_target_axis"] == [1, 0, 2]
    assert result["lossless_reorientation_feasible"] is True
    assert result["aligned_intensity_metrics"]["unequal_voxels"] == 0


def test_different_origin_is_genuine_physical_grid_failure():
    array = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    native = sitk_image(array)
    converted = sitk_image(array, origin=(0.25, 0.0, 0.0))
    result = compare_physical_space(loaded(native), loaded(converted))
    assert result["lossless_reorientation_feasible"] is False
    assert "SPATIAL_FAIL_GRID" in result["spatial_statuses"]
    assert result["pixel_status"] == "PIXEL_PASS_EXACT"
    assert result["center_difference_mm"] == pytest.approx(0.25)


def test_inconsistent_direction_matrix_is_not_axis_reorientation():
    array = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    native = sitk_image(array)
    angle = np.deg2rad(30.0)
    direction = (
        np.cos(angle), -np.sin(angle), 0.0,
        np.sin(angle), np.cos(angle), 0.0,
        0.0, 0.0, 1.0,
    )
    converted = sitk_image(array, direction=direction)
    result = compare_physical_space(loaded(native), loaded(converted))
    assert result["axis_relation"]["detected"] is False
    assert result["lossless_reorientation_feasible"] is False
    assert "SPATIAL_FAIL_GRID" in result["spatial_statuses"]


def test_partial_and_no_physical_overlap():
    array = np.arange(16, dtype=np.float32).reshape(1, 4, 4)
    native = sitk_image(array)
    partial = sitk_image(array, origin=(2.0, 0.0, 0.0))
    partial_result = compare_physical_space(loaded(native), loaded(partial))
    assert partial_result["physical_overlap"]["percentage_of_a"] == pytest.approx(50.0)
    assert partial_result["sampling_overlap_percentage_of_a"] == pytest.approx(50.0)

    separate = sitk_image(array, origin=(100.0, 0.0, 0.0))
    no_overlap = compare_physical_space(loaded(native), loaded(separate))
    assert no_overlap["physical_overlap"]["percentage_of_a"] == 0.0
    assert no_overlap["alignment_method"] == "none"
    assert "SPATIAL_FAIL_COVERAGE" in no_overlap["spatial_statuses"]


def test_intensity_scale_and_offset_fit_after_alignment():
    converted_array = np.arange(1, 25, dtype=np.float32).reshape(2, 3, 4)
    native_array = 2.5 * converted_array + 7.0
    result = compare_physical_space(
        loaded(sitk_image(native_array)), loaded(sitk_image(converted_array))
    )
    fit = result["diagnostic_affine_intensity_fit"]
    assert fit["scale"] == pytest.approx(2.5)
    assert fit["offset"] == pytest.approx(7.0)
    assert fit["residual_rmse"] == pytest.approx(0.0, abs=1e-10)
    assert "DIFFERENT_INTENSITY_SCALING" in result["statuses"]


def test_symmetric_normalized_difference_zero_handling_signs_and_finiteness():
    native = np.array([0.0, 5.0, 0.0, 8.0, -4.0, 1e308, np.nan, np.inf])
    converted = np.array([0.0, 0.0, 5.0, 4.0, -2.0, -1e308, 1.0, 2.0])
    difference = _symmetric_normalized_difference(native, converted)

    assert difference == pytest.approx(
        [0.0, 200.0, -200.0, 200.0 / 3.0, -200.0 / 3.0, 200.0, 0.0, 0.0]
    )
    assert np.all(np.isfinite(difference))


def test_symmetric_normalized_difference_is_antisymmetric():
    native = np.array([[1.0, -2.0], [3.0, -4.0]])
    converted = np.array([[2.0, -1.0], [-3.0, 4.0]])
    forward = _symmetric_normalized_difference(native, converted)
    reverse = _symmetric_normalized_difference(converted, native)
    assert forward == pytest.approx(-reverse)


def test_diagnostic_figure_panel_arrangement_and_fit_annotation():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    native = np.arange(1, 13, dtype=float).reshape(3, 4)
    aligned = (native - 1.25) / 0.875
    fit = {"available": True, "scale": 0.875123456, "offset": 1.251234567}
    correlation = 0.9980917654
    physical = {
        "alignment_method": "lossless_discrete_axis_reorientation",
        "interpolated": False,
        "sampling_overlap_percentage_of_a": 100.0,
        "aligned_intensity_metrics": {"pearson_correlation": correlation},
        "diagnostic_affine_intensity_fit": fit,
    }
    fig, axes = _create_diagnostic_figure(
        plt,
        {"image_identifier": "volume:0000"},
        physical,
        native,
        aligned,
        np.ones_like(native, dtype=bool),
        aligned.ravel(),
        native.ravel(),
    )
    try:
        assert axes.shape == (2, 3)
        assert [[axis.get_title() for axis in row] for row in axes] == [
            ["Native (A physical grid)", "Converted aligned (B)", "Numerical summary"],
            [
                "Absolute intensity difference: A − B",
                "Normalized difference (%)\n200(A−B)/(|A|+|B|), display saturated at ±100%",
                "Native vs aligned converted intensity",
            ],
        ]
        annotation = "\n".join(text.get_text() for text in axes[1, 2].texts)
        assert "Native = 0.875123 × Converted + 1.251235" in annotation
        assert "Pearson r = 0.998092" in annotation
        assert len(axes[1, 2].lines) == 1
        assert axes[0, 0].images[0].get_clim() == axes[0, 1].images[0].get_clim()
        assert axes[1, 0].images[0].get_clim()[0] == -axes[1, 0].images[0].get_clim()[1]
        assert axes[1, 1].images[0].get_clim() == (-100.0, 100.0)
    finally:
        plt.close(fig)


def test_intensity_fit_annotation_reports_unavailable_fit():
    annotation = _format_intensity_fit_annotation({"available": False}, 0.9123456789)
    assert annotation == "Least-squares fit unavailable\nPearson r = 0.912346"


def test_intensity_fit_annotation_formats_negative_intercept():
    annotation = _format_intensity_fit_annotation(
        {"available": True, "scale": 0.5, "offset": -0.5},
        1.0,
    )
    assert annotation == (
        "Native = 0.500000 × Converted − 0.500000\n"
        "Pearson r = 1.000000"
    )


def test_integer_valued_vs_fractional_float_diagnostics():
    converted = np.arange(24, dtype=np.int16).reshape(2, 3, 4)
    native = converted.astype(np.float32) + 0.25
    result = compare_physical_space(loaded(sitk_image(native)), loaded(sitk_image(converted)))
    assert result["intensity_distribution_a"]["integer_valued_percentage"] == 0.0
    assert result["intensity_distribution_b"]["integer_valued_percentage"] == 100.0
    assert result["diagnostic_affine_intensity_fit"]["offset"] == pytest.approx(0.25)


def test_simpleitk_numpy_axis_order_and_physical_formula():
    array = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)
    direction = (0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    image = sitk_image(array, spacing=(2.0, 3.0, 4.0), origin=(10.0, 20.0, 30.0), direction=direction)
    numpy_array = sitk.GetArrayFromImage(image)
    assert numpy_array.shape == (image.GetSize()[2], image.GetSize()[1], image.GetSize()[0])
    assert numpy_array[1, 2, 3] == image.GetPixel((3, 2, 1))
    index = np.asarray((3.0, 2.0, 1.0))
    matrix = np.asarray(direction).reshape(3, 3) @ np.diag(image.GetSpacing())
    expected = np.asarray(image.GetOrigin()) + matrix @ index
    assert image.TransformIndexToPhysicalPoint((3, 2, 1)) == pytest.approx(expected)


def test_resampling_uses_overlap_mask_and_identity_physical_transform():
    x = np.arange(4, dtype=np.float32)
    array = np.broadcast_to(x, (1, 4, 4)).copy()
    native = sitk_image(array)
    shifted = sitk_image(array, origin=(1.5, 0.0, 0.0))
    result = compare_physical_space(loaded(native), loaded(shifted))
    assert result["alignment_method"] == "identity_physical_resampling_linear"
    assert result["interpolated"] is True
    assert result["sampling_overlap_percentage_of_a"] == pytest.approx(75.0)
    assert result["aligned_intensity_metrics"]["voxel_count"] == 12


def test_spatial_equivalence_within_configured_tolerance_prefers_lossless_reindexing():
    array = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    native = sitk_image(array)
    converted = sitk_image(
        np.flip(array, axis=1),
        origin=(0.0005, 2.0, 0.0),
        direction=(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0),
    )
    accepted = compare_physical_space(
        loaded(native), loaded(converted), spatial_tolerance_mm=0.001
    )
    assert accepted["spatial_status"] == "SPATIAL_PASS_TOLERANCE"
    assert accepted["alignment_method"] == "lossless_discrete_axis_reorientation"
    assert accepted["interpolated"] is False
    assert accepted["axis_relation"]["corresponding_sample_max_discrepancy_mm"] == pytest.approx(0.0005)

    rejected = compare_physical_space(
        loaded(native), loaded(converted), spatial_tolerance_mm=0.0001
    )
    assert "SPATIAL_FAIL_GRID" in rejected["spatial_statuses"]


def test_high_intensity_correlation_does_not_override_spatial_failure():
    array = np.arange(1, 25, dtype=np.float32).reshape(2, 3, 4)
    native = sitk_image(array)
    converted = sitk_image(array * 2.0 + 1.0, origin=(10.0, 0.0, 0.0))
    result = compare_physical_space(
        loaded(native), loaded(converted), pixel_correlation_threshold=0.99
    )
    assert result["pixel_status"] == "PIXEL_SIMILAR"
    assert result["spatial_status"] in {"SPATIAL_FAIL_GRID", "SPATIAL_FAIL_COVERAGE"}


def test_different_voxel_spacing_fails_spatial_grid():
    array = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    native = sitk_image(array, spacing=(1.0, 1.0, 1.0))
    converted = sitk_image(array, spacing=(1.1, 1.0, 1.0))
    result = compare_physical_space(loaded(native), loaded(converted))
    assert "SPATIAL_FAIL_GRID" in result["spatial_statuses"]
    assert result["pixel_status"] == "PIXEL_PASS_EXACT"


def test_independent_summary_statuses_and_report_fields(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    native = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    converted = native * 2.0 + 1.0
    write_nhdr(a, image_name("a", 0), native)
    write_nhdr(b, image_name("b", 0), converted)
    summary, rows, exit_code = compare_directories(
        a,
        b,
        tmp_path / "reports",
        pixel_correlation_threshold=0.99,
    )
    assert exit_code == 0
    assert summary["coverage_status"] == "COVERAGE_PASS"
    assert summary["spatial_overall_status"] == "SPATIAL_PASS_EXACT"
    assert summary["pixel_overall_status"] == "PIXEL_SIMILAR"
    assert rows[0]["spatial_status"] == "SPATIAL_PASS_EXACT"
    assert rows[0]["pixel_status"] == "PIXEL_SIMILAR"
    loaded_summary = json.loads((tmp_path / "reports" / "comparison_summary.json").read_text())
    assert loaded_summary["validation_configuration"]["validation_target"] == "both"
    with (tmp_path / "reports" / "image_comparison.csv").open(newline="") as stream:
        csv_row = next(csv.DictReader(stream))
    assert csv_row["spatial_status"] == "SPATIAL_PASS_EXACT"
    assert csv_row["pixel_status"] == "PIXEL_SIMILAR"
    report = (tmp_path / "reports" / "comparison_report.md").read_text()
    assert "## B. Spatial correspondence" in report
    assert "## C. Pixel-value correspondence" in report


def test_validation_target_controls_exit_code(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    array = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    write_nhdr(a, image_name("a", 0), array, origin=(0, 0, 0))
    write_nhdr(b, image_name("b", 0), array, origin=(10, 0, 0))
    _, _, spatial_exit = compare_directories(
        a, b, tmp_path / "spatial", validation_target="spatial"
    )
    _, _, pixel_exit = compare_directories(
        a, b, tmp_path / "pixels", validation_target="pixels"
    )
    _, _, both_exit = compare_directories(
        a, b, tmp_path / "both", validation_target="both"
    )
    assert spatial_exit == 1
    assert pixel_exit == 0
    assert both_exit == 1


def test_cli_validation_target_pixels_ignores_spatial_failure_for_exit(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    array = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    write_nhdr(a, image_name("a", 0), array, origin=(0, 0, 0))
    write_nhdr(b, image_name("b", 0), array, origin=(10, 0, 0))
    process = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--dir-a", str(a),
            "--dir-b", str(b),
            "--output-dir", str(tmp_path / "reports"),
            "--validation-target", "pixels",
            "--spatial-tolerance-mm", "0.001",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert process.returncode == 0, process.stderr
    assert "Spatial: SPATIAL_FAIL" in process.stdout
    assert "Pixels: PIXEL_PASS_EXACT" in process.stdout

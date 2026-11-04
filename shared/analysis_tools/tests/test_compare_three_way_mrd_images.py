from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk

from shared.analysis_tools import compare_motion_control_images as two_way
from shared.analysis_tools import compare_three_way_mrd_images as three_way


def _loaded(image: sitk.Image) -> two_way.LoadedImage:
    array = sitk.GetArrayFromImage(image)
    return two_way.LoadedImage(array, two_way._sitk_geometry(image, array.shape), image)


def _write_image(
    path: Path,
    array: np.ndarray,
    *,
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> None:
    image = sitk.GetImageFromArray(array)
    image.SetSpacing((1.0, 1.0, 1.0))
    image.SetOrigin(origin)
    sitk.WriteImage(image, str(path), useCompression=False)


def _write_raw_4d_nrrd(path: Path, array_tzyx: np.ndarray) -> None:
    t, z, y, x = array_tzyx.shape
    header = (
        "NRRD0005\n"
        "# synthetic dcm2niix-like list-axis NRRD\n"
        "type: short\n"
        "dimension: 4\n"
        "space: right-anterior-superior\n"
        f"sizes: {x} {y} {z} {t}\n"
        "space directions: (-2,0,0) (0,2,0) (0,0,2) none\n"
        "kinds: space space space list\n"
        "endian: little\n"
        "encoding: raw\n"
        "space origin: (10,20,30)\n\n"
    ).encode("ascii")
    path.write_bytes(header + np.asarray(array_tzyx, dtype="<i2").tobytes(order="C"))


def _make_fire_tree(root: Path, *, z_size: int = 2, volumes: int = 2) -> None:
    _write_image(
        root / "protocol_volume_0000_20261104T000000.nrrd",
        np.arange(z_size * 6, dtype=np.float32).reshape(z_size, 2, 3),
    )
    for volume in range(1, volumes):
        for ordinal in range(z_size):
            anatomical = z_size - ordinal - 1
            _write_image(
                root / f"protocol_volume_{volume:04d}_slice_{ordinal:04d}.nrrd",
                np.full((1, 2, 3), volume * 10 + anatomical, dtype=np.float32),
                origin=(0.0, 0.0, float(anatomical)),
            )


def test_load_4d_list_axis_nrrd_and_ras_to_lps(tmp_path: Path) -> None:
    values = np.arange(3 * 2 * 3 * 4, dtype=np.int16).reshape(3, 2, 3, 4)
    path = tmp_path / "dicom.nrrd"
    _write_raw_4d_nrrd(path, values)

    loaded = three_way.load_dicom_nrrd(path)

    assert loaded.array_tzyx.shape == (3, 2, 3, 4)
    np.testing.assert_array_equal(loaded.array_tzyx, values)
    assert loaded.header["fields"]["space"] == "right-anterior-superior"
    assert loaded.header["ras_to_lps_verification"]["verified"] is True
    assert loaded.image.GetOrigin() == pytest.approx((-10.0, -20.0, 30.0))
    np.testing.assert_allclose(three_way.ras_to_lps_point((10, 20, 30)), (-10, -20, 30))
    np.testing.assert_allclose(
        three_way.ras_to_lps_direction(np.eye(3)), np.diag((-1.0, -1.0, 1.0))
    )


def test_volume_offsets_require_complete_unambiguous_coverage() -> None:
    established = three_way.determine_volume_offset([1, 2, 3], 3)
    assert established == {
        "status": "established",
        "fire_to_dicom_offset": -1,
        "evidence": "complete contiguous FIRE coverage equals the DICOM time dimension",
    }
    assert three_way.determine_volume_offset([1, 2], 4)["status"] == "indeterminate"
    assert three_way.determine_volume_offset([0, 0, 1], 3)["reason"].startswith("duplicate")


def test_fire_slice_order_is_derived_from_physical_origins(tmp_path: Path) -> None:
    _make_fire_tree(tmp_path, z_size=4)

    dataset = three_way.discover_fire_dataset(tmp_path, "synthetic")

    assert dataset.slice_ordinal_to_anatomical == {0: 3, 1: 2, 2: 1, 3: 0}
    reconstructed = dataset.volume_array(1)
    assert [int(reconstructed[z, 0, 0]) for z in range(4)] == [10, 11, 12, 13]


def test_ambiguous_slice_mapping_is_rejected(tmp_path: Path) -> None:
    _make_fire_tree(tmp_path, z_size=2)
    duplicate = tmp_path / "protocol_volume_0001_slice_0001.nrrd"
    _write_image(duplicate, np.zeros((1, 2, 3), dtype=np.float32), origin=(0.0, 0.0, 1.0))

    with pytest.raises(three_way.ComparisonError, match="one-to-one mapping"):
        three_way.discover_fire_dataset(tmp_path, "synthetic")


def test_axis_permutation_and_reflection_are_losslessly_compared() -> None:
    target = sitk.GetImageFromArray(np.arange(24, dtype=np.int16).reshape(2, 3, 4))
    source = sitk.PermuteAxes(target, [1, 0, 2])
    source = sitk.Flip(source, [True, False, False])

    row, details = three_way.compare_loaded(
        "test", "volume:0000", _loaded(target), _loaded(source), 0.001, 0.99, None
    )

    assert details["physical"]["axis_relation"]["source_axis_to_target_axis"] == [1, 0, 2]
    assert details["physical"]["axis_relation"]["reflected_source_axes"] == [0]
    assert row["lossless_reorientation_feasible"] is True
    assert row["numeric_values_equal"] is True


def test_equal_values_with_different_dtypes_are_reported_separately() -> None:
    values = np.arange(24, dtype=np.int16).reshape(2, 3, 4)
    integer_image = sitk.GetImageFromArray(values)
    float_image = sitk.GetImageFromArray(values.astype(np.float32))

    row, _details = three_way.compare_loaded(
        "test", "volume:0000", _loaded(integer_image), _loaded(float_image),
        0.001, 0.99, None,
    )

    assert row["numeric_values_equal"] is True
    assert row["exact_voxel_equal_including_dtype"] is False
    assert row["byte_equal_after_common_cast"] is True


def test_affine_intensity_scaling_is_diagnostic_not_primary_normalization() -> None:
    source_values = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    source = sitk.GetImageFromArray(source_values)
    target = sitk.GetImageFromArray(2.0 * source_values + 3.0)
    row, details = three_way.compare_loaded(
        "test", "volume:0000", _loaded(target), _loaded(source), 0.001, 0.99, None
    )
    aggregate = three_way.aggregate_pair("test", [row], [details])

    assert row["numeric_values_equal"] is False
    assert row["mean_absolute_error"] > 0
    fit = aggregate["diagnostic_affine_intensity_fit"]
    assert fit["scale"] == pytest.approx(2.0)
    assert fit["offset"] == pytest.approx(3.0)
    assert fit["residual_rmse"] == pytest.approx(0.0, abs=1e-6)


def test_half_voxel_origin_is_spatial_failure_despite_perfect_correlation() -> None:
    values = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    target = sitk.GetImageFromArray(values)
    source = sitk.GetImageFromArray(values)
    source.SetOrigin((0.5, 0.5, 0.0))

    row, _details = three_way.compare_loaded(
        "test", "volume:0000", _loaded(target), _loaded(source), 0.001, 0.99, None
    )

    assert row["spatial_status"] == "SPATIAL_FAIL_GRID"
    assert row["pixel_comparison_basis"] == "index_only_signed_axis_diagnostic"
    assert row["pearson_correlation"] == pytest.approx(1.0)
    assert row["max_coordinate_discrepancy_mm"] == pytest.approx(np.sqrt(0.5))


def test_missing_fire_records_fail_before_comparison(tmp_path: Path) -> None:
    fire_root = tmp_path / "fire"
    fire_root.mkdir()
    _make_fire_tree(fire_root, z_size=2, volumes=2)
    fire = three_way.discover_fire_dataset(fire_root, "synthetic")
    dicom_path = tmp_path / "dicom.nrrd"
    _write_raw_4d_nrrd(dicom_path, np.zeros((3, 2, 2, 3), dtype=np.int16))
    dicom = three_way.load_dicom_nrrd(dicom_path)

    with pytest.raises(three_way.ComparisonError, match="correspondence incomplete"):
        three_way.analyze_fire_dicom_pair(
            "test", fire, dicom, True, 0.001, 0.99
        )


def test_report_files_are_generated(tmp_path: Path) -> None:
    metric = {
        "pearson_correlation": 1.0,
        "mean_absolute_error": 0.0,
        "rmse": 0.0,
    }
    pair = {
        "spatial_status": "SPATIAL_PASS_EXACT",
        "pixel_metrics": metric,
        "pixel_comparison_basis_counts": {"lossless_physical_reorientation": 1},
        "all_numeric_values_equal": True,
        "all_exact_voxel_equal_including_dtype": True,
        "all_byte_equal_after_common_cast": True,
        "mean_corresponding_coordinate_discrepancy_mm": 0.0,
        "maximum_corresponding_coordinate_discrepancy_mm": 0.0,
        "mean_physical_support_iou_percentage": 100.0,
        "lossless_reorientation_feasible_count": 1,
        "pair_count": 1,
    }
    summary = {
        "inputs": {"native_fire": "A", "dicom_nrrd": "B", "converted_fire": "C"},
        "spatial_tolerance_mm": 0.001,
        "discovery": {
            "native": {"record_count": 1},
            "converted": {"record_count": 1},
            "dicom": {"stored_sizes_xyzt": [2, 2, 1, 1]},
        },
        "correspondence": {
            "temporal": {"fire_to_dicom_offset": 0},
            "native_slice_ordinal_to_anatomical_z": {0: 0},
        },
        "intensity_distributions": {
            name: {
                "dtypes": ["<f4"], "minimum": 0.0, "maximum": 1.0,
                "mean": 0.5, "standard_deviation": 0.5,
                "integer_valued_percentage": 100.0,
            }
            for name in ("native", "dicom", "converted")
        },
        "pairwise_comparisons": {name: dict(pair) for name in three_way.PAIR_NAMES},
        "interpretation": {
            "half_voxel_geometry_diagnostic": {
                "expected_in_plane_half_voxel_diagonal_mm": 0.70710678,
            }
        },
    }
    row = {"comparison": "native_vs_dicom", "identifier": "volume:0000"}

    three_way.write_outputs(tmp_path, summary, [row])

    assert (tmp_path / "comparison_report.md").is_file()
    assert (tmp_path / "pairwise_comparison.csv").is_file()
    parsed = json.loads((tmp_path / "comparison_summary.json").read_text())
    assert parsed["spatial_tolerance_mm"] == 0.001

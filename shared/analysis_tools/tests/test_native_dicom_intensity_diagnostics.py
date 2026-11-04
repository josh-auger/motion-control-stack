from __future__ import annotations

import numpy as np
import pydicom
import pytest
from pydicom.sequence import Sequence

from shared.analysis_tools import diagnose_native_dicom_intensities as diagnostic


def test_integer_valued_float32_is_distinguished_from_storage_dtype() -> None:
    summary = diagnostic.value_distribution(np.arange(8, dtype=np.float32))

    assert summary["storage_dtype"] == "<f4"
    assert summary["integer_valued_percentage"] == 100.0


def test_fractional_native_values_are_measured() -> None:
    summary = diagnostic.value_distribution(
        np.asarray([0.0, 1.25, 2.5, 3.75], dtype=np.float32)
    )

    assert summary["integer_valued_percentage"] == 25.0
    assert summary["distance_to_nearest_integer_percentiles"]["100"] == 0.5


def test_pure_rounding_model() -> None:
    native = np.asarray([1.2, 1.6, 3.49, 3.51], dtype=np.float32)
    dicom = np.asarray([1, 2, 3, 4], dtype=np.int16)

    models = diagnostic.quantization_models(native, dicom)

    assert models["round"]["exact_agreement_percentage"] == 100.0
    assert models["round"]["rmse"] == 0.0


def test_truncation_is_evaluated_independently() -> None:
    native = np.asarray([1.9, 2.9, 3.1], dtype=np.float32)
    dicom = np.asarray([1, 2, 3], dtype=np.int16)

    models = diagnostic.quantization_models(native, dicom)

    assert models["truncate"]["exact_agreement_percentage"] == 100.0
    assert models["round"]["exact_agreement_percentage"] < 100.0


def test_linear_scaling_fit_recovers_coefficients() -> None:
    source = np.arange(20, dtype=np.float32)
    target = 2.5 * source + 3.0

    fit = diagnostic.fit_linear_relationship(target, source)

    assert fit["available"] is True
    assert fit["scale"] == pytest.approx(2.5)
    assert fit["intercept"] == pytest.approx(3.0)
    assert fit["r_squared"] == pytest.approx(1.0)
    assert fit["residual_rmse"] == pytest.approx(0.0, abs=1e-12)


def test_linear_scaling_followed_by_quantization() -> None:
    dicom = np.arange(20, dtype=np.int16)
    native = (1.75 * dicom + 4.25).astype(np.float32)
    fit = diagnostic.fit_linear_relationship(native, dicom)

    models = diagnostic.scaling_quantization_models(native, dicom, fit)

    assert models["round"]["exact_agreement_percentage"] == 100.0
    assert models["round"]["rmse"] == 0.0


def test_nonlinear_relationship_leaves_residual_error() -> None:
    source = np.linspace(0, 10, 101, dtype=np.float32)
    target = source * source

    fit = diagnostic.fit_linear_relationship(target, source)

    assert fit["available"] is True
    assert fit["r_squared"] < 1.0
    assert fit["residual_rmse"] > 1.0


def test_per_volume_scaling_drift_is_detected() -> None:
    dicom = np.tile(np.arange(12, dtype=np.float32).reshape(1, 1, 3, 4), (3, 1, 1, 1))
    native = np.empty_like(dicom)
    for volume, scale in enumerate((1.0, 1.5, 2.0)):
        native[volume] = scale * dicom[volume] + volume

    rows = diagnostic.fit_per_volume(native, dicom)
    stability = diagnostic.summarize_fit_stability(rows)

    assert [row["scale"] for row in rows] == pytest.approx([1.0, 1.5, 2.0])
    assert [row["intercept"] for row in rows] == pytest.approx([0.0, 1.0, 2.0])
    assert stability["scale"]["minimum"] == pytest.approx(1.0)
    assert stability["scale"]["maximum"] == pytest.approx(2.0)


def test_spatially_varying_scaling_is_not_explained_by_one_fit() -> None:
    source = np.arange(1, 17, dtype=np.float32).reshape(2, 2, 4)
    target = source.copy()
    target[1] *= 2.0

    global_fit = diagnostic.fit_linear_relationship(target, source)
    first = diagnostic.fit_linear_relationship(target[0], source[0])
    second = diagnostic.fit_linear_relationship(target[1], source[1])

    assert global_fit["residual_rmse"] > 0
    assert first["scale"] == pytest.approx(1.0)
    assert second["scale"] == pytest.approx(2.0)


def test_dicom_zero_floor_error_is_separated_from_positive_voxels() -> None:
    native = np.asarray([0.0, 10.0, 20.0, 30.0], dtype=np.float32)
    dicom = np.asarray([0, 0, 20, 30], dtype=np.int16)

    result = diagnostic.zero_floor_diagnostics(native, dicom)

    assert result["dicom_zero_voxel_percentage"] == 50.0
    assert result["raw_error_where_dicom_zero"]["rmse"] == pytest.approx(np.sqrt(50.0))
    assert result["raw_error_where_dicom_positive"]["rmse"] == 0.0
    assert result["dicom_zero_squared_error_percentage"] == 100.0


def test_constant_source_is_reported_as_unavailable() -> None:
    fit = diagnostic.fit_linear_relationship(
        np.arange(8, dtype=np.float32), np.ones(8, dtype=np.float32)
    )

    assert fit["available"] is False
    assert "zero variance" in fit["reason"]


def test_empty_or_invalid_mask_is_rejected() -> None:
    values = np.arange(4, dtype=np.float32)

    with pytest.raises(diagnostic.DiagnosticError, match="No finite voxels"):
        diagnostic.fit_linear_relationship(values, values, np.zeros(4, dtype=bool))
    with pytest.raises(diagnostic.DiagnosticError, match="mask shape"):
        diagnostic.fit_linear_relationship(values, values, np.ones(3, dtype=bool))


def test_enhanced_dicom_pixel_transform_metadata() -> None:
    dataset = pydicom.Dataset()
    frame_groups = []
    for _ in range(2):
        group = pydicom.Dataset()
        transform = pydicom.Dataset()
        transform.RescaleSlope = "1.25"
        transform.RescaleIntercept = "-3"
        transform.RescaleType = "US"
        group.PixelValueTransformationSequence = Sequence([transform])
        frame_groups.append(group)
    dataset.PerFrameFunctionalGroupsSequence = Sequence(frame_groups)

    metadata = diagnostic.extract_dicom_scaling_metadata(dataset)

    assert metadata["slopes"] == [1.25, 1.25]
    assert metadata["intercepts"] == [-3.0, -3.0]
    assert metadata["rescale_types"] == ["US", "US"]
    assert metadata["real_world_value_mapping_present"] is False

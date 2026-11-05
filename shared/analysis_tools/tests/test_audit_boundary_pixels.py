from __future__ import annotations

from pathlib import Path

import numpy as np
import pydicom
import pytest
import SimpleITK as sitk
from pydicom.sequence import Sequence

from shared.analysis_tools.audit_boundary_pixels import (
    AuditError,
    align_dicom_volume,
    boundary_masks,
    dicom_volume_geometry,
    index_enhanced_dataset,
    pair_metrics,
    squared_error_attribution,
    zero_mask_classes,
)
from shared.analysis_tools.compare_motion_control_images import LoadedImage, _sitk_geometry


def loaded(image: sitk.Image) -> LoadedImage:
    array = sitk.GetArrayFromImage(image)
    return LoadedImage(array, _sitk_geometry(image, array.shape), image)


def frame_group(
    temporal: int,
    anatomical_slice: int,
    position: tuple[float, float, float],
) -> pydicom.Dataset:
    group = pydicom.Dataset()
    content = pydicom.Dataset()
    content.TemporalPositionIndex = temporal
    content.InStackPositionNumber = anatomical_slice
    group.FrameContentSequence = Sequence([content])

    plane_position = pydicom.Dataset()
    plane_position.ImagePositionPatient = list(position)
    group.PlanePositionSequence = Sequence([plane_position])

    plane_orientation = pydicom.Dataset()
    plane_orientation.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    group.PlaneOrientationSequence = Sequence([plane_orientation])

    measures = pydicom.Dataset()
    measures.PixelSpacing = [1, 1]
    measures.SpacingBetweenSlices = 1
    measures.SliceThickness = 1
    group.PixelMeasuresSequence = Sequence([measures])
    return group


def enhanced_dataset(groups: list[pydicom.Dataset]) -> pydicom.Dataset:
    dataset = pydicom.Dataset()
    dataset.NumberOfFrames = len(groups)
    dataset.Rows = 2
    dataset.Columns = 3
    dataset.PerFrameFunctionalGroupsSequence = Sequence(groups)
    return dataset


def test_boundary_masks_follow_xy_coordinates_and_prevent_corner_double_counting():
    masks = boundary_masks((3, 4))
    assert np.flatnonzero(masks["first_x"]) == pytest.approx([0, 4, 8])
    assert np.flatnonzero(masks["last_x"]) == pytest.approx([3, 7, 11])
    assert np.flatnonzero(masks["first_y"]) == pytest.approx([0, 1, 2, 3])
    assert np.flatnonzero(masks["last_y"]) == pytest.approx([8, 9, 10, 11])
    assert np.count_nonzero(masks["all_outer_boundary"]) == 2 * 3 + 2 * 4 - 4
    assert np.count_nonzero(masks["all_outer_two_pixel_boundary"]) == 12


def test_exact_zero_mask_classification():
    native = np.array([[0.0, 2.0], [0.0, 4.0]])
    converted = np.array([[0.0, 0.0], [3.0, 5.0]])
    classes = zero_mask_classes(native, converted)
    assert classes["both_zero"].tolist() == [[True, False], [False, False]]
    assert classes["converted_zero_native_nonzero"].tolist() == [[False, True], [False, False]]
    assert classes["native_zero_converted_nonzero"].tolist() == [[False, False], [True, False]]
    assert classes["both_nonzero"].tolist() == [[False, False], [False, True]]


def test_squared_error_attribution_is_exhaustive_without_boundary_overlap():
    native = np.zeros((4, 4), dtype=float)
    converted = np.zeros((4, 4), dtype=float)
    native[0, 0] = 2.0  # Boundary converted-zero/native-nonzero: SSE 4.
    native[1, 1] = 3.0  # Interior converted-zero/native-nonzero: SSE 9.
    converted[1, 2] = 4.0  # Interior native-zero/converted-nonzero: SSE 16.
    native[2, 1], converted[2, 1] = 7.0, 2.0  # Interior both nonzero: SSE 25.
    result = squared_error_attribution(native, converted)
    assert result["total_squared_error"] == 54.0
    assert result["boundary_squared_error"] == 4.0
    assert result["interior_converted_zero_native_nonzero_squared_error"] == 9.0
    assert result["interior_native_zero_converted_nonzero_squared_error"] == 16.0
    assert result["interior_both_nonzero_squared_error"] == 25.0
    assert result["zero_mask_disagreement_squared_error"] == 29.0
    assert result["converted_zero_native_nonzero_squared_error"] == 13.0


def test_all_zero_and_constant_images_have_defined_errors_and_no_correlation():
    zero = np.zeros((3, 4), dtype=float)
    metrics = pair_metrics(zero, zero)
    assert metrics["mae"] == 0.0
    assert metrics["rmse"] == 0.0
    assert metrics["pearson_correlation"] is None

    first = np.ones((3, 4), dtype=float)
    second = np.full((3, 4), 2.0)
    constant = pair_metrics(first, second)
    assert constant["mae"] == 1.0
    assert constant["rmse"] == 1.0
    assert constant["pearson_correlation"] is None

    read_only_mask = np.broadcast_to(np.array([[True, False, True, False]]), zero.shape)
    masked = pair_metrics(zero, zero, read_only_mask)
    assert masked["pixel_count"] == 6


def test_enhanced_dicom_metadata_orders_frames_without_file_or_frame_order_assumption():
    dataset = enhanced_dataset([
        frame_group(7, 2, (0, 0, 1)),
        frame_group(7, 1, (0, 0, 0)),
    ])
    record = index_enhanced_dataset(dataset, Path("enhanced.dcm"))
    assert record.volume == 6
    assert record.temporal_position_index == 7
    assert record.frame_indices_by_anatomical_slice == (1, 0)

    decoded = np.stack((np.full((2, 3), 20), np.full((2, 3), 10)))
    ordered = decoded[np.asarray(record.frame_indices_by_anatomical_slice)]
    assert ordered[:, 0, 0].tolist() == [10, 20]
    image, residual = dicom_volume_geometry(
        dataset, record.frame_indices_by_anatomical_slice, ordered
    )
    assert image.GetSize() == (3, 2, 2)
    assert residual == pytest.approx(0.0)


@pytest.mark.parametrize(
    "groups, message",
    [
        ([frame_group(1, 1, (0, 0, 0)), frame_group(1, 1, (0, 0, 1))], "duplicate"),
        ([frame_group(1, 1, (0, 0, 0)), frame_group(2, 2, (0, 0, 1))], "multiple temporal"),
    ],
)
def test_ambiguous_enhanced_dicom_metadata_is_rejected(groups, message):
    with pytest.raises(AuditError, match=message):
        index_enhanced_dataset(enhanced_dataset(groups), Path("ambiguous.dcm"))


def test_missing_enhanced_dicom_frame_metadata_is_rejected():
    group = frame_group(1, 1, (0, 0, 0))
    del group.FrameContentSequence[0].InStackPositionNumber
    with pytest.raises(AuditError, match="lacks TemporalPositionIndex or InStackPositionNumber"):
        index_enhanced_dataset(enhanced_dataset([group]), Path("missing.dcm"))


def test_dicom_to_native_alignment_handles_lossless_axis_reflection():
    target_array = np.arange(2 * 3 * 4, dtype=np.uint16).reshape(2, 3, 4)
    native_image = sitk.GetImageFromArray(target_array)
    source_array = np.flip(target_array, axis=2)
    source_image = sitk.GetImageFromArray(source_array)
    source_image.SetOrigin((3.0, 0.0, 0.0))
    source_image.SetDirection((-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

    class FakeSeries:
        def load_volume(self, volume):
            assert volume == 4
            return source_array, source_image, 0.0

    aligned, relation = align_dicom_volume(FakeSeries(), 4, loaded(native_image), 0.001)
    assert relation["source_axis_signs"] == [-1, 1, 1]
    assert np.array_equal(aligned, target_array)

from __future__ import annotations

import numpy as np
import pytest

from shared.analysis_tools import audit_fire_geometry as audit


def test_even_identity_center_conventions_and_offset_sign() -> None:
    position = np.zeros(3)
    size = (72, 72, 1)
    spacing = (3.0, 3.0, 3.0)
    directions = np.eye(3)

    fire = audit.origin_from_center(position, size, spacing, directions)
    n_over_two = audit.origin_from_center(
        position,
        size,
        spacing,
        directions,
        convention="fov_midpoint",
        inplane_only=True,
    )

    np.testing.assert_allclose(fire, (-106.5, -106.5, 0.0))
    np.testing.assert_allclose(n_over_two, (-108.0, -108.0, 0.0))
    np.testing.assert_allclose(fire - n_over_two, (1.5, 1.5, 0.0))
    assert np.linalg.norm(fire - n_over_two) == pytest.approx(3.0 / np.sqrt(2.0))


def test_odd_dimensions_remain_half_voxel_apart_algebraically() -> None:
    fire = audit.origin_from_center(
        (0, 0, 0), (71, 71, 1), (3, 3, 3), np.eye(3)
    )
    n_over_two = audit.origin_from_center(
        (0, 0, 0),
        (71, 71, 1),
        (3, 3, 3),
        np.eye(3),
        convention="fov_midpoint",
        inplane_only=True,
    )

    np.testing.assert_allclose(fire, (-105.0, -105.0, 0.0))
    np.testing.assert_allclose(n_over_two, (-106.5, -106.5, 0.0))
    np.testing.assert_allclose(fire - n_over_two, (1.5, 1.5, 0.0))


def test_oblique_orientation_uses_physical_index_vectors() -> None:
    angle = np.deg2rad(30.0)
    directions = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    spacing = np.array((2.0, 4.0, 5.0))
    expected = 0.5 * (
        directions[:, 0] * spacing[0] + directions[:, 1] * spacing[1]
    )

    actual = audit.inplane_half_voxel_vector(directions, spacing)

    np.testing.assert_allclose(actual, expected)


def test_axis_reflection_changes_the_corresponding_offset_sign() -> None:
    directions = np.diag((1.0, -1.0, 1.0))

    half = audit.inplane_half_voxel_vector(directions, (3.0, 3.0, 3.0))

    np.testing.assert_allclose(half, (1.5, -1.5, 0.0))


def test_ras_to_lps_converts_points_and_each_vector_column() -> None:
    point = np.array((10.0, -20.0, 30.0))
    directions = np.array(
        [[1.0, 0.1, 0.2], [0.3, 1.0, 0.4], [0.5, 0.6, 1.0]]
    )

    np.testing.assert_allclose(audit.ras_to_lps_point(point), (-10.0, 20.0, 30.0))
    np.testing.assert_allclose(
        audit.ras_to_lps_directions(directions),
        np.diag((-1.0, -1.0, 1.0)) @ directions,
    )


def test_invalid_center_convention_is_rejected() -> None:
    with pytest.raises(ValueError, match="unknown center convention"):
        audit.origin_from_center(
            (0, 0, 0), (2, 2, 1), (1, 1, 1), np.eye(3), convention="unknown"
        )


@pytest.mark.skipif(not audit.DEFAULT_MRD.exists(), reason="REF capture is not available")
def test_real_ref_header_independently_reconstructs_fire_origin() -> None:
    native = audit.read_mrd_geometry(audit.DEFAULT_MRD, 10)
    observed = audit.parse_nrrd_header(audit.DEFAULT_FIRE_NRRD)
    direction = native["fire_meta_unit_directions_lps"]
    spacing = native["spacing_xyz_mm"]
    size = native["size_xyz"]
    position = native["position_lps"]

    # Independent expanded equation rather than a call back into origin_from_center.
    expected = position.copy()
    for axis in range(3):
        expected -= direction[:, axis] * spacing[axis] * (size[axis] - 1) / 2.0

    np.testing.assert_allclose(observed["origin"], expected, atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(
        observed["scaled_directions"], direction * spacing, atol=1e-12, rtol=0.0
    )


@pytest.mark.skipif(
    not audit.DEFAULT_DICOM_DIR.exists(), reason="enhanced REF DICOM is not available"
)
def test_real_ref_dicom_ipp_exhibits_inplane_n_over_two_relation() -> None:
    native = audit.read_mrd_geometry(audit.DEFAULT_MRD, 10)
    dicom = audit.read_enhanced_dicom_geometry(audit.DEFAULT_DICOM_DIR)
    steps = native["fire_meta_unit_directions_lps"] * native["spacing_xyz_mm"]
    offset_indices = np.linalg.solve(
        steps,
        native["position_lps"] - dicom["image_position_patient_lps"],
    )

    # DICOM DS values are rounded in the public tags, hence the 2e-4 tolerance.
    np.testing.assert_allclose(offset_indices, (36.0, 36.0, 0.0), atol=2e-4, rtol=0.0)


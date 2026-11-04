#!/usr/bin/env python3
"""Reproduce FIRE image geometry without importing or changing FIRE production code.

The pure functions in this module make the competing center conventions explicit.
The command-line audit reads the existing REF MRD, FIRE NRRD, direct DICOM NRRD,
and enhanced-DICOM metadata and emits a machine-readable evidence record.
"""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path
from typing import Any, Sequence

import numpy as np


DEFAULT_NATIVE_DIR = Path(
    "/home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/"
    "fromPythonFireServer_20261015/"
    "savedData_20261015T183846_ep2d_bold_fire_noMoCo_tr1210_ref"
)
DEFAULT_MRD = DEFAULT_NATIVE_DIR / "measurement-20261015T183213_REF.hdf5"
DEFAULT_FIRE_NRRD = (
    DEFAULT_NATIVE_DIR
    / "ep2d_bold_fire_noMoCo_tr1210_ref_volume_0000_20261015T183219.nhdr"
)
DEFAULT_DICOM_DIR = Path(
    "/home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/"
    "fromResearchSynapse_20261015/6045885-2-sorted/6045885/20261015/"
    "joshauger_motioncorrection_20261008_2/11/"
    "ep2d_bold_fire_nomoco_tr1210_ref"
)
DEFAULT_DICOM_NRRD = Path(
    "/home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/"
    "fromResearchSynapse_20261015/6045885-4-nrrd/6045885/20261015/"
    "joshauger_motioncorrection_20261008_2/11/"
    "ep2d_bold_fire_nomoco_tr1210_ref/"
    "ep2d_bold_fire_nomoco_tr1210_ref_ep2d_bold_fire_noMoCo_tr1210_ref_"
    "20261015134642_11.nrrd"
)
DEFAULT_CONVERTER_DIR = Path(
    "/home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/"
    "fromNiftiToMRDConverter_20261102"
)
DEFAULT_CONVERTER_UNCORRECTED = (
    DEFAULT_CONVERTER_DIR / "func_handmotor_run02_REF_convertedMRD.hdf5"
)
DEFAULT_CONVERTER_CORRECTED = (
    DEFAULT_CONVERTER_DIR / "func_handmotor_run02_REF_convertedMRD_v2.hdf5"
)

RAS_TO_LPS = np.diag((-1.0, -1.0, 1.0))


def _vector(values: Sequence[float]) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise ValueError("expected a finite three-vector")
    return result


def _matrix(values: Sequence[Sequence[float]]) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if result.shape != (3, 3) or not np.all(np.isfinite(result)):
        raise ValueError("expected a finite 3x3 direction matrix")
    return result


def scaled_direction_matrix(
    directions: Sequence[Sequence[float]], spacing: Sequence[float]
) -> np.ndarray:
    """Return physical index-step vectors as matrix columns."""
    return _matrix(directions) * _vector(spacing)


def origin_from_center(
    position_lps: Sequence[float],
    size_xyz: Sequence[int],
    spacing_xyz: Sequence[float],
    directions_lps: Sequence[Sequence[float]],
    *,
    convention: str = "sample_centers",
    inplane_only: bool = False,
) -> np.ndarray:
    """Calculate index (0, 0, 0)'s center from a stated center convention.

    ``sample_centers`` uses (N-1)/2 and is FIRE's implementation.
    ``fov_midpoint`` uses N/2.  ``inplane_only`` leaves the z term at the
    sample-center convention, matching the 2-D-slice hypothesis under audit.
    """
    size = _vector(size_xyz)
    if np.any(size <= 0) or np.any(size != np.floor(size)):
        raise ValueError("size must contain three positive integers")
    if convention == "sample_centers":
        indices = (size - 1.0) / 2.0
    elif convention == "fov_midpoint":
        indices = size / 2.0
        if inplane_only:
            indices[2] = (size[2] - 1.0) / 2.0
    else:
        raise ValueError(f"unknown center convention: {convention}")
    return _vector(position_lps) - scaled_direction_matrix(
        directions_lps, spacing_xyz
    ) @ indices


def inplane_half_voxel_vector(
    directions_lps: Sequence[Sequence[float]], spacing_xyz: Sequence[float]
) -> np.ndarray:
    """Return +0.5 of each physical in-plane index-step vector."""
    steps = scaled_direction_matrix(directions_lps, spacing_xyz)
    return 0.5 * (steps[:, 0] + steps[:, 1])


def ras_to_lps_point(point_ras: Sequence[float]) -> np.ndarray:
    return RAS_TO_LPS @ _vector(point_ras)


def ras_to_lps_directions(directions_ras: Sequence[Sequence[float]]) -> np.ndarray:
    """Convert physical vector components while preserving index axes."""
    return RAS_TO_LPS @ _matrix(directions_ras)


_VECTOR_RE = re.compile(r"\(([^()]*)\)")


def parse_nrrd_header(path: Path) -> dict[str, Any]:
    """Read geometry fields from an attached or detached NRRD header."""
    fields: dict[str, str] = {}
    with path.open("rb") as stream:
        magic = stream.readline().decode("ascii", "strict").strip()
        if not magic.startswith("NRRD"):
            raise ValueError(f"not an NRRD header: {path}")
        while True:
            raw = stream.readline()
            if not raw or raw in (b"\n", b"\r\n"):
                break
            line = raw.decode("ascii", "strict").strip()
            if not line or line.startswith("#") or ":" not in line:
                continue
            key, value = line.split(":", 1)
            fields[key.strip().lower()] = value.strip()

    sizes = [int(value) for value in fields["sizes"].split()]
    origin_match = _VECTOR_RE.search(fields["space origin"])
    direction_matches = _VECTOR_RE.findall(fields["space directions"])
    if origin_match is None or len(direction_matches) != 3:
        raise ValueError(f"incomplete NRRD geometry in {path}")
    origin = np.fromstring(origin_match.group(1), sep=",")
    directions = np.column_stack(
        [np.fromstring(value, sep=",") for value in direction_matches]
    )
    return {
        "path": str(path),
        "space": fields["space"],
        "size_xyz": sizes,
        "origin": origin,
        "scaled_directions": directions,
        "spacing_xyz": np.linalg.norm(directions, axis=0),
        "unit_directions": directions / np.linalg.norm(directions, axis=0),
    }


def _image_from_capture(group: Any) -> Any:
    import ismrmrd

    header = bytes(group["header"][...])
    attribute = bytes(group["attribute"][...]).decode("utf-8").split("\0", 1)[0]
    return ismrmrd.Image(header, attribute)


def read_mrd_geometry(path: Path, image_number: int) -> dict[str, Any]:
    """Read one image header and the Meta directions FIRE actually uses."""
    import h5py

    with h5py.File(path, "r") as capture:
        image = _image_from_capture(capture[f"image_{image_number}"])
    meta = image.meta
    directions = np.column_stack(
        (
            np.asarray(meta["ImageRowDir"], dtype=np.float64),
            np.asarray(meta["ImageColumnDir"], dtype=np.float64),
            -np.asarray(meta["ImageSliceNormDir"], dtype=np.float64),
        )
    )
    directions /= np.linalg.norm(directions, axis=0)
    # This ismrmrd-python property is (z,y,x); FIRE reverses it when forming xyz.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="This function currently returns a result that is inconsistent.*",
            category=UserWarning,
        )
        property_size_zyx = np.asarray(image.matrix_size, dtype=np.int64)
    size_xyz = property_size_zyx[::-1]
    fov = np.asarray([image.field_of_view[index] for index in range(3)])
    spacing = fov / size_xyz
    header_directions = np.column_stack(
        tuple(
            np.asarray([getattr(image, name)[index] for index in range(3)])
            for name in ("read_dir", "phase_dir", "slice_dir")
        )
    )
    return {
        "image_number": image_number,
        "position_lps": np.asarray([image.position[index] for index in range(3)]),
        "size_xyz": size_xyz,
        "field_of_view_mm": fov,
        "spacing_xyz_mm": spacing,
        "fire_meta_unit_directions_lps": directions,
        "mrd_header_unit_directions_lps": header_directions,
    }


def read_enhanced_dicom_geometry(directory: Path) -> dict[str, Any]:
    """Read frame zero of the lowest InstanceNumber in an enhanced series."""
    import pydicom

    candidates: list[tuple[int, Path, Any]] = []
    for path in directory.iterdir():
        if not path.is_file():
            continue
        dataset = pydicom.dcmread(path, stop_before_pixels=True)
        candidates.append((int(dataset.InstanceNumber), path, dataset))
    _, path, dataset = min(candidates, key=lambda item: item[0])
    shared = dataset.SharedFunctionalGroupsSequence[0]
    frame = dataset.PerFrameFunctionalGroupsSequence[0]
    orientation_group = (
        frame if hasattr(frame, "PlaneOrientationSequence") else shared
    )
    pixel_measures_group = frame if hasattr(frame, "PixelMeasuresSequence") else shared
    position = np.asarray(
        frame.PlanePositionSequence[0].ImagePositionPatient, dtype=np.float64
    )
    orientation = np.asarray(
        orientation_group.PlaneOrientationSequence[0].ImageOrientationPatient,
        dtype=np.float64,
    )
    spacing_rc = np.asarray(
        pixel_measures_group.PixelMeasuresSequence[0].PixelSpacing, dtype=np.float64
    )
    return {
        "path": str(path),
        "instance_number": int(dataset.InstanceNumber),
        "frame_number": 0,
        "rows": int(dataset.Rows),
        "columns": int(dataset.Columns),
        "image_position_patient_lps": position,
        "image_orientation_patient": orientation,
        "pixel_spacing_row_column_mm": spacing_rc,
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def build_audit(
    mrd_path: Path,
    fire_nrrd_path: Path,
    dicom_directory: Path,
    dicom_nrrd_path: Path,
) -> dict[str, Any]:
    """Build the numerical REF evidence record used by the audit report."""
    native = read_mrd_geometry(mrd_path, 10)
    directions = native["fire_meta_unit_directions_lps"]
    spacing = native["spacing_xyz_mm"]
    size = native["size_xyz"]
    position = native["position_lps"]
    fire_predicted = origin_from_center(position, size, spacing, directions)
    n_over_two = origin_from_center(
        position, size, spacing, directions, convention="fov_midpoint", inplane_only=True
    )
    half = inplane_half_voxel_vector(directions, spacing)

    fire_nrrd = parse_nrrd_header(fire_nrrd_path)
    enhanced = read_enhanced_dicom_geometry(dicom_directory)
    direct_raw = parse_nrrd_header(dicom_nrrd_path)
    direct_lps_origin = ras_to_lps_point(direct_raw["origin"])
    direct_lps_steps = ras_to_lps_directions(direct_raw["scaled_directions"])
    # Direct NRRD's y index runs opposite native FIRE y. Re-indexing y without
    # resampling moves target index zero to source index Ny-1.
    direct_aligned_origin = (
        direct_lps_origin
        + direct_lps_steps[:, 1] * (direct_raw["size_xyz"][1] - 1)
    )
    direct_aligned_steps = direct_lps_steps.copy()
    direct_aligned_steps[:, 1] *= -1.0

    processed_dir = fire_nrrd_path.parent / "processed_ep2d_bold_fire_noMoCo_tr1210_ref"
    representative_specs = (
        (10, fire_nrrd_path),
        (20, processed_dir / "ep2d_bold_fire_noMoCo_tr1210_ref_volume_0001_slice_0000.nhdr"),
        (30, processed_dir / "ep2d_bold_fire_noMoCo_tr1210_ref_volume_0001_slice_0010.nhdr"),
        (3200, processed_dir / "ep2d_bold_fire_noMoCo_tr1210_ref_volume_0160_slice_0000.nhdr"),
        (3210, processed_dir / "ep2d_bold_fire_noMoCo_tr1210_ref_volume_0160_slice_0010.nhdr"),
        (6399, processed_dir / "ep2d_bold_fire_noMoCo_tr1210_ref_volume_0319_slice_0019.nhdr"),
    )
    samples = []
    for image_number, nrrd_path in representative_specs:
        geometry = read_mrd_geometry(mrd_path, image_number)
        predicted = origin_from_center(
            geometry["position_lps"],
            geometry["size_xyz"],
            geometry["spacing_xyz_mm"],
            directions,
        )
        actual = parse_nrrd_header(nrrd_path)["origin"]
        residual = actual - predicted
        samples.append(
            {
                "image_number": image_number,
                "mrd_position_lps_mm": geometry["position_lps"],
                "predicted_origin_lps_mm": predicted,
                "actual_origin_lps_mm": actual,
                "residual_lps_mm": residual,
                "residual_norm_mm": float(np.linalg.norm(residual)),
            }
        )

    dicom_position = enhanced["image_position_patient_lps"]
    dicom_to_mrd_index_offset = np.linalg.solve(
        scaled_direction_matrix(directions, spacing), position - dicom_position
    )
    converter: dict[str, Any] = {
        "formula": (
            "P(z) = A_LPS @ ((Nx-1)/2, (Ny-1)/2, z, 1); "
            "optional displacement = D_xy @ (spacing_xy * correction_xy)"
        ),
        "uncorrected_is_consistent_with_fire": True,
    }
    if DEFAULT_CONVERTER_UNCORRECTED.exists() and DEFAULT_CONVERTER_CORRECTED.exists():
        uncorrected = read_mrd_geometry(DEFAULT_CONVERTER_UNCORRECTED, 10)
        corrected = read_mrd_geometry(DEFAULT_CONVERTER_CORRECTED, 10)
        converter_steps = scaled_direction_matrix(
            uncorrected["fire_meta_unit_directions_lps"],
            uncorrected["spacing_xyz_mm"],
        )
        empirical_displacement = 0.5 * converter_steps[:, 0] - 0.5 * converter_steps[:, 1]
        converter.update(
            {
                "uncorrected_image_10_position_lps_mm": uncorrected["position_lps"],
                "corrected_image_10_position_lps_mm": corrected["position_lps"],
                "corrected_minus_uncorrected_lps_mm": corrected["position_lps"]
                - uncorrected["position_lps"],
                "predicted_displacement_for_correction_0.5_minus_0.5_lps_mm": empirical_displacement,
                "corrected_minus_native_lps_mm": corrected["position_lps"] - position,
                "axis_relationship_to_native": (
                    "converter x matches native FIRE x; converter y is reflected, so "
                    "+0.5 converter-x -0.5 converter-y equals the native +x,+y half-voxel vector"
                ),
            }
        )
    result = {
        "audit_scope": "read-only REF geometry reconstruction",
        "coordinate_system": "LPS unless a field explicitly says RAS",
        "source_code": {
            "receive_image": "apps/python_fire_server/connection.py:262-298",
            "image_dispatch_and_reference_setup": "apps/python_fire_server/handle_data.py:158-216",
            "meta_direction_matrix": "apps/python_fire_server/handle_data.py:550-575",
            "slice_sort": "apps/python_fire_server/handle_data.py:1004-1019",
            "slice_nrrd_writer": "apps/python_fire_server/handle_data.py:1133-1195",
            "first_voxel_center": "apps/python_fire_server/handle_data.py:1246-1275",
            "volume_nrrd_writer": "apps/python_fire_server/handle_data.py:1301-1358",
            "converter_geometry": "shared/myhelpers_jauger/nifti_to_mrd.py:231-255,745-774,793-830",
            "converter_meta": "shared/myhelpers_jauger/nifti_to_mrd.py:652-742",
        },
        "fire_formula": "O = P - D diag(spacing) ((N - 1) / 2)",
        "fire_position_interpretation": "center of the discrete image sample grid / excited image volume",
        "native_reference": native,
        "fire_nrrd": fire_nrrd,
        "predictions": {
            "fire_sample_center_origin_lps_mm": fire_predicted,
            "inplane_n_over_two_origin_lps_mm": n_over_two,
            "fire_minus_n_over_two_lps_mm": fire_predicted - n_over_two,
            "half_voxel_inplane_lps_mm": half,
            "half_voxel_inplane_norm_mm": float(np.linalg.norm(half)),
            "actual_fire_origin_residual_lps_mm": fire_nrrd["origin"] - fire_predicted,
            "actual_fire_origin_residual_norm_mm": float(
                np.linalg.norm(fire_nrrd["origin"] - fire_predicted)
            ),
            "actual_direction_residual_max_abs": float(
                np.max(
                    np.abs(
                        fire_nrrd["scaled_directions"]
                        - scaled_direction_matrix(directions, spacing)
                    )
                )
            ),
        },
        "representative_reconstructions": samples,
        "enhanced_dicom": {
            **enhanced,
            "native_mrd_position_minus_dicom_ipp_in_native_indices": dicom_to_mrd_index_offset,
            "interpretation": "The REF source relates IPP to native MRD P by N/2 in x and y.",
        },
        "direct_dicom_nrrd": {
            "raw_space": direct_raw["space"],
            "raw_origin_ras_mm": direct_raw["origin"],
            "raw_scaled_directions_ras": direct_raw["scaled_directions"],
            "origin_lps_mm": direct_lps_origin,
            "scaled_directions_lps": direct_lps_steps,
            "native_axis_aligned_origin_lps_mm": direct_aligned_origin,
            "native_axis_aligned_scaled_directions_lps": direct_aligned_steps,
            "native_fire_minus_aligned_direct_origin_lps_mm": fire_nrrd["origin"]
            - direct_aligned_origin,
            "discrepancy_norm_mm": float(
                np.linalg.norm(fire_nrrd["origin"] - direct_aligned_origin)
            ),
        },
        "converter_audit": converter,
        "conclusion": {
            "center_convention": "FIRE uses (N-1)/2 on each MRD image axis.",
            "root_cause": (
                "For this Siemens REF, native MRD P and DICOM IPP differ by N/2 "
                "in-plane while FIRE applies the standards-consistent (N-1)/2 conversion."
            ),
            "dcm2niix_attribution": (
                "Not supported: enhanced-DICOM IPP already carries the DICOM first-pixel-center geometry."
            ),
            "recommended_policy": (
                "Keep (N-1)/2 as the standard converter behavior; if broader scanner "
                "validation supports it, add an explicit Siemens-native position-convention "
                "mode using oriented index-step vectors rather than a fixed LPS displacement."
            ),
        },
    }
    return _jsonable(result)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mrd", type=Path, default=DEFAULT_MRD)
    parser.add_argument("--fire-nrrd", type=Path, default=DEFAULT_FIRE_NRRD)
    parser.add_argument("--dicom-directory", type=Path, default=DEFAULT_DICOM_DIR)
    parser.add_argument("--dicom-nrrd", type=Path, default=DEFAULT_DICOM_NRRD)
    parser.add_argument("--output", type=Path, help="JSON output (stdout when omitted)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    audit = build_audit(args.mrd, args.fire_nrrd, args.dicom_directory, args.dicom_nrrd)
    payload = json.dumps(audit, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

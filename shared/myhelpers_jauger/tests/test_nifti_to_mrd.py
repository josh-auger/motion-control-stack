from __future__ import annotations

import json
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import ismrmrd
import nibabel as nib
import numpy as np
import pytest


HELPER_DIR = Path(__file__).resolve().parents[1]
FIRE_SERVER_DIR = Path(__file__).resolve().parents[3] / "apps" / "python_fire_server"
sys.path.insert(0, str(HELPER_DIR))
import nifti_to_mrd as converter  # noqa: E402


MRD_NS = "http://www.ismrm.org/ISMRMRD"


def minimal_xml() -> bytes:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<ismrmrdHeader xmlns="{MRD_NS}">
  <version>1</version>
  <studyInformation><studyInstanceUID>reference</studyInstanceUID></studyInformation>
  <measurementInformation>
    <measurementID>reference</measurementID>
    <protocolName>reference</protocolName>
  </measurementInformation>
  <acquisitionSystemInformation>
    <systemVendor>ReferenceVendor</systemVendor>
    <systemModel>ReferenceModel</systemModel>
    <systemFieldStrength_T>3</systemFieldStrength_T>
  </acquisitionSystemInformation>
  <encoding>
    <encodedSpace>
      <matrixSize><x>2</x><y>2</y><z>1</z></matrixSize>
      <fieldOfView_mm><x>2</x><y>2</y><z>1</z></fieldOfView_mm>
    </encodedSpace>
    <reconSpace>
      <matrixSize><x>2</x><y>2</y><z>1</z></matrixSize>
      <fieldOfView_mm><x>2</x><y>2</y><z>1</z></fieldOfView_mm>
    </reconSpace>
    <encodingLimits>
      <slice><minimum>0</minimum><maximum>0</maximum><center>0</center></slice>
      <repetition><minimum>0</minimum><maximum>0</maximum><center>0</center></repetition>
    </encodingLimits>
    <trajectory>cartesian</trajectory>
  </encoding>
  <sequenceParameters>
    <TR>1000</TR><TE>30</TE><flipAngle_deg>45</flipAngle_deg>
  </sequenceParameters>
  <userParameters />
</ismrmrdHeader>""".encode()


def write_bytes(parent: h5py.Group, name: str, payload: bytes) -> None:
    parent.create_dataset(name, data=np.frombuffer(payload, dtype=np.uint8))


def make_reference(path: Path, series=(0,)) -> Path:
    xml = minimal_xml()
    # Ensure the synthetic fixture itself is a valid ISMRMRD header.
    ismrmrd.xsd.CreateFromDocument(xml.decode())
    with h5py.File(path, "w") as h5file:
        write_bytes(h5file, "Config File", bytes(1024))
        write_bytes(h5file, "Metadata XML", xml + bytes((0,)))
        for number, series_index in enumerate(series):
            header = ismrmrd.ImageHeader()
            header.version = 1
            header.flags = 15
            header.data_type = ismrmrd.DATATYPE_FLOAT
            header.matrix_size = (2, 2, 1)
            header.field_of_view = (2, 2, 1)
            header.channels = 1
            header.image_type = ismrmrd.IMTYPE_MAGNITUDE
            header.image_series_index = series_index
            group = h5file.create_group(f"image_{number}")
            write_bytes(group, "header", bytes(header))
            write_bytes(group, "attribute", ismrmrd.Meta().serialize().encode())
            write_bytes(group, "data", np.zeros((2, 2), np.float32).tobytes())
    return path


def default_metadata(z_size: int) -> dict:
    return {
        "RepetitionTime": 1.2,
        "EchoTime": 0.031,
        "FlipAngle": 52,
        "MagneticFieldStrength": 3,
        "Manufacturer": "TestVendor",
        "ManufacturersModelName": "TestModel",
        "ProtocolName": "synthetic_protocol",
        "SeriesDescription": "synthetic_series",
        "SliceThickness": 2.5,
        "SliceTiming": [0.1 * z for z in range(z_size)],
        "PhaseEncodingDirection": "j-",
    }


def make_input(
    directory: Path,
    data: np.ndarray,
    affine: np.ndarray | None = None,
    metadata: dict | None = None,
    slope: float | None = None,
    intercept: float | None = None,
) -> tuple[Path, Path]:
    if affine is None:
        affine = np.diag([2.0, 3.0, 2.5, 1.0])
    image = nib.Nifti1Image(data, affine)
    if slope is not None or intercept is not None:
        image.header.set_slope_inter(slope, intercept)
    nifti = directory / "scan.nii.gz"
    nib.save(image, nifti)
    sidecar = directory / "scan.json"
    sidecar.write_text(json.dumps(metadata if metadata is not None else default_metadata(data.shape[2])))
    return nifti, sidecar


def run_conversion(
    tmp_path: Path,
    data: np.ndarray,
    affine: np.ndarray | None = None,
    metadata: dict | None = None,
    use_reference: bool = False,
    **kwargs,
):
    nifti, sidecar = make_input(tmp_path, data, affine, metadata, **kwargs)
    output = tmp_path / "converted.hdf5"
    if use_reference:
        reference = make_reference(tmp_path / "reference.hdf5")
        result = converter.convert(nifti, reference, output, json_path=sidecar)
    else:
        result = converter.convert(nifti, output, json_path=sidecar)
    return result, nifti, output


def test_basic_4d_conversion_and_image_count(tmp_path):
    data = np.arange(4 * 5 * 3 * 2, dtype=np.int16).reshape(4, 5, 3, 2)
    result, _nifti, output = run_conversion(tmp_path, data)

    assert result.image_count == 6
    assert result.shape == data.shape
    assert result.dtype == np.dtype(np.float32)
    with h5py.File(output, "r") as h5file:
        assert set(("Config File", "Metadata XML")).issubset(h5file)
        assert converter._bytes(h5file["Config File"]).startswith(b"fire_emitter_dcm")
        assert len(converter.image_group_names(h5file)) == 6
        header, attribute, pixels = converter.read_image(h5file["image_0"])
        assert tuple(header.matrix_size) == (4, 5, 1)
        assert pixels.shape == (1, 1, 5, 4)
        meta = converter.validate_image_attribute(attribute, header)
        assert set(converter.REQUIRED_DIRECTION_ATTRIBUTES).issubset(meta.keys())
        assert meta["SlicePosLightMarker"] == [
            f"{float(value):.18f}" for value in header.position
        ]
        assert meta["RepetitionTime"] == "1200"
        assert meta["EchoTime"] == "31"


def test_slice_timing_controls_order_and_indices(tmp_path):
    data = np.zeros((2, 3, 4, 2), dtype=np.int16)
    metadata = default_metadata(4)
    metadata["SliceTiming"] = [0.3, 0.0, 0.2, 0.1]
    result, _nifti, output = run_conversion(tmp_path, data, metadata=metadata)

    assert result.slice_order == (1, 3, 2, 0)
    with h5py.File(output, "r") as h5file:
        indices = []
        for name in converter.image_group_names(h5file):
            header, _attribute, _pixels = converter.read_image(h5file[name])
            indices.append((header.repetition, header.slice, header.image_index))
    assert indices == [
        (0, 1, 0), (0, 3, 0), (0, 2, 0), (0, 0, 0),
        (1, 1, 0), (1, 3, 0), (1, 2, 0), (1, 0, 0),
    ]


def test_scaled_voxel_intensities_are_preserved(tmp_path):
    raw = np.arange(3 * 4 * 2 * 2, dtype=np.int16).reshape(3, 4, 2, 2)
    result, nifti, output = run_conversion(tmp_path, raw, slope=2.0, intercept=10.0)
    expected = np.asanyarray(nib.load(nifti).dataobj)
    restored = converter.reconstruct_4d(output, expected.shape)

    assert result.dtype == np.dtype(np.float64)
    np.testing.assert_array_equal(restored, expected)
    assert restored[0, 0, 0, 0] == 10


def test_mrd_readback_reconstructs_every_volume(tmp_path):
    rng = np.random.default_rng(4)
    data = rng.integers(-1000, 1000, size=(5, 4, 3, 4), dtype=np.int16)
    _result, _nifti, output = run_conversion(tmp_path, data)
    restored = converter.reconstruct_4d(output, data.shape)
    np.testing.assert_array_equal(restored, data)


def test_ras_to_lps_geometry(tmp_path):
    data = np.zeros((4, 6, 3, 1), dtype=np.int16)
    affine = np.array(
        [[2.0, 0, 0, 10.0], [0, 3.0, 0, 20.0], [0, 0, 4.0, 30.0], [0, 0, 0, 1.0]]
    )
    metadata = default_metadata(3)
    metadata["SliceThickness"] = 3.5
    _result, _nifti, output = run_conversion(tmp_path, data, affine, metadata)

    with h5py.File(output, "r") as h5file:
        header, _attribute, _pixels = converter.read_image(h5file["image_0"])
    np.testing.assert_allclose(header.read_dir, [-1, 0, 0], atol=1e-6)
    np.testing.assert_allclose(header.phase_dir, [0, -1, 0], atol=1e-6)
    np.testing.assert_allclose(header.slice_dir, [0, 0, 1], atol=1e-6)
    np.testing.assert_allclose(header.position, [-13, -27.5, 30], atol=1e-6)
    np.testing.assert_allclose(header.field_of_view, [8, 18, 3.5], atol=1e-6)


@pytest.mark.filterwarnings("ignore:builtin type Swig.*:DeprecationWarning")
def test_fire_coordinate_frame_extraction_uses_generated_attributes(tmp_path):
    data = np.zeros((4, 6, 3, 1), dtype=np.int16)
    angle = 0.31
    affine = np.array(
        [
            [2 * math.cos(angle), -3 * math.sin(angle), 0, 10],
            [2 * math.sin(angle), 3 * math.cos(angle), 0, 20],
            [0, 0, 4, 30],
            [0, 0, 0, 1],
        ]
    )
    _result, _nifti, output = run_conversion(tmp_path, data, affine)

    with h5py.File(output, "r") as h5file:
        header_bytes = converter._bytes(h5file["image_0/header"])
        attribute = converter._bytes(h5file["image_0/attribute"])
    image = ismrmrd.Image(header_bytes, attribute.decode("utf-8"))

    sys.path.insert(0, str(FIRE_SERVER_DIR))
    try:
        from handle_data import handleData
    finally:
        sys.path.pop(0)
    fire_matrix = handleData.__new__(handleData).get_image_direction_matrix(image)
    expected = np.column_stack(
        (image.read_dir, image.phase_dir, image.slice_dir)
    )
    np.testing.assert_allclose(fire_matrix, expected, atol=1e-6)


def test_missing_image_row_direction_reproduces_original_fire_failure(tmp_path):
    data = np.zeros((2, 2, 1, 1), dtype=np.int16)
    _result, _nifti, output = run_conversion(tmp_path, data)
    with h5py.File(output, "r") as h5file:
        header_bytes = converter._bytes(h5file["image_0/header"])
    header = ismrmrd.ImageHeader.from_buffer_copy(header_bytes)
    incomplete_meta = ismrmrd.Meta()
    incomplete_meta["Keep_image_geometry"] = 1
    incomplete_attribute = incomplete_meta.serialize().encode("utf-8")
    image = ismrmrd.Image(header_bytes, incomplete_attribute.decode("utf-8"))

    sys.path.insert(0, str(FIRE_SERVER_DIR))
    try:
        from handle_data import handleData
    finally:
        sys.path.pop(0)
    with pytest.raises(KeyError, match="ImageRowDir"):
        handleData.__new__(handleData).get_image_direction_matrix(image)
    with pytest.raises(converter.ConversionError, match="ImageRowDir"):
        converter.validate_image_attribute(incomplete_attribute, header)


@pytest.mark.parametrize(
    "linear",
    [
        np.array([[0, -2, 0], [3, 0, 0], [0, 0, 4.0]]),  # permutation/reflection
        np.array(
            [
                [2 * math.cos(0.35), -3 * math.sin(0.35), 0],
                [2 * math.sin(0.35), 3 * math.cos(0.35), 0],
                [0, 0, -4.0],
            ]
        ),  # oblique/reflected
    ],
)
def test_axis_permutations_reflections_and_oblique_affines(tmp_path, linear):
    data = np.zeros((3, 4, 2, 1), dtype=np.int16)
    affine = np.eye(4)
    affine[:3, :3] = linear
    affine[:3, 3] = (7, -8, 9)
    _result, _nifti, output = run_conversion(tmp_path, data, affine)
    expected = converter.derive_geometry(affine, data.shape, default_metadata(2))

    with h5py.File(output, "r") as h5file:
        for name in converter.image_group_names(h5file):
            header, _attribute, _pixels = converter.read_image(h5file[name])
            np.testing.assert_allclose(
                np.column_stack((header.read_dir, header.phase_dir, header.slice_dir)),
                expected.directions,
                atol=1e-6,
            )


def test_json_metadata_mapping(tmp_path):
    data = np.zeros((2, 3, 2, 4), dtype=np.int16)
    metadata = default_metadata(2)
    metadata["MultibandAccelerationFactor"] = 2
    _result, _nifti, output = run_conversion(tmp_path, data, metadata=metadata)

    with h5py.File(output, "r") as h5file:
        xml = converter._without_nul(converter._bytes(h5file["Metadata XML"]))
    root = ET.fromstring(xml)
    ns = {"m": MRD_NS}
    assert root.findtext("m:measurementInformation/m:protocolName", namespaces=ns) == "synthetic_protocol"
    assert root.findtext("m:acquisitionSystemInformation/m:systemVendor", namespaces=ns) == "TestVendor"
    assert root.findtext("m:acquisitionSystemInformation/m:systemModel", namespaces=ns) == "TestModel"
    assert float(root.findtext("m:sequenceParameters/m:TR", namespaces=ns)) == 1200
    assert float(root.findtext("m:sequenceParameters/m:TE", namespaces=ns)) == 31
    assert float(root.findtext("m:sequenceParameters/m:flipAngle_deg", namespaces=ns)) == 52
    assert root.findtext("m:encoding/m:encodingLimits/m:slice/m:maximum", namespaces=ns) == "1"
    assert root.findtext("m:encoding/m:encodingLimits/m:repetition/m:maximum", namespaces=ns) == "3"
    parsed = ismrmrd.xsd.CreateFromDocument(xml.decode("utf-8"))
    assert parsed.encoding[0].trajectory.value == "cartesian"
    assert parsed.userParameters.userParameterLong[0].name == "MultiBandFactor"
    assert parsed.userParameters.userParameterLong[0].value == 2


@pytest.mark.filterwarnings("ignore:builtin type Swig.*:DeprecationWarning")
def test_reference_free_xml_satisfies_fire_metadata_accessors(tmp_path):
    data = np.zeros((2, 2, 4, 1), dtype=np.int16)
    metadata = default_metadata(4)
    metadata["MultibandAccelerationFactor"] = 2
    _result, _nifti, output = run_conversion(tmp_path, data, metadata=metadata)

    with h5py.File(output, "r") as h5file:
        xml = converter._without_nul(converter._bytes(h5file["Metadata XML"]))
    parsed = ismrmrd.xsd.CreateFromDocument(xml.decode("utf-8"))

    sys.path.insert(0, str(FIRE_SERVER_DIR))
    try:
        from handle_data import handleData
    finally:
        sys.path.pop(0)
    handler = handleData.__new__(handleData)
    handler.metadata_object = parsed

    assert handler.get_protocol_name_from_metadata() == "synthetic_protocol"
    assert handler.get_smsfactor_from_metadata() == 2
    assert handler.get_nslices_per_volume_from_metadata() == 4
    spacing_parameters = {
        parameter.name: parameter.value
        for parameter in parsed.userParameters.userParameterDouble
    }
    assert spacing_parameters["SpacingBetweenSlices"] == pytest.approx(2.5)
    # FIRE's current accessor lowercases the name before comparing it to a
    # mixed-case literal, so it follows its documented volume-derived fallback.
    assert handler.get_spacing_between_slices_from_metadata() is None


def test_missing_optional_json_fields_uses_safe_defaults(tmp_path):
    data = np.zeros((2, 2, 2, 1), dtype=np.int16)
    result, _nifti, output = run_conversion(tmp_path, data, metadata={})
    assert result.slice_order == (0, 1)
    with h5py.File(output, "r") as h5file:
        header, _attribute, _pixels = converter.read_image(h5file["image_0"])
    assert header.acquisition_time_stamp == 0
    assert header.field_of_view[2] == pytest.approx(2.5)


def test_reference_free_invalid_metadata_is_rejected(tmp_path):
    data = np.zeros((2, 2, 2, 1), dtype=np.int16)
    metadata = default_metadata(2)
    metadata["MultibandAccelerationFactor"] = 0
    nifti, sidecar = make_input(tmp_path, data, metadata=metadata)

    with pytest.raises(converter.ConversionError, match="MultibandAccelerationFactor"):
        converter.convert(nifti, tmp_path / "out.hdf5", json_path=sidecar)


def test_cli_reference_argument_is_optional():
    args = converter.build_argument_parser().parse_args(
        ["--nifti", "scan.nii.gz", "--output", "converted.hdf5"]
    )
    assert args.reference_mrd is None


def test_legacy_reference_assisted_path_matches_reference_free_images(tmp_path):
    data = np.arange(3 * 4 * 2 * 2, dtype=np.int16).reshape(3, 4, 2, 2)
    nifti, sidecar = make_input(tmp_path, data)
    reference = make_reference(tmp_path / "reference.hdf5")
    reference_free = tmp_path / "reference_free.hdf5"
    assisted = tmp_path / "assisted.hdf5"

    converter.convert(nifti, reference_free, json_path=sidecar)
    converter.convert(nifti, reference, assisted, json_path=sidecar)

    with h5py.File(reference_free, "r") as free_file, h5py.File(assisted, "r") as assisted_file:
        assert converter.image_group_names(free_file) == converter.image_group_names(assisted_file)
        for name in converter.image_group_names(free_file):
            for field in ("header", "attribute", "data"):
                assert converter._bytes(free_file[f"{name}/{field}"]) == converter._bytes(
                    assisted_file[f"{name}/{field}"]
                )


def test_invalid_nifti_dimensions_are_rejected(tmp_path):
    data = np.zeros((2, 2, 2), dtype=np.int16)
    nifti = tmp_path / "three_dimensional.nii.gz"
    nib.save(nib.Nifti1Image(data, np.eye(4)), nifti)
    sidecar = tmp_path / "three_dimensional.json"
    sidecar.write_text("{}")
    reference = make_reference(tmp_path / "reference.hdf5")

    with pytest.raises(converter.ConversionError, match="Expected a 4D NIfTI"):
        converter.convert(nifti, reference, tmp_path / "out.hdf5", sidecar)


def test_incompatible_multiseries_reference_is_rejected(tmp_path):
    data = np.zeros((2, 2, 2, 1), dtype=np.int16)
    nifti, sidecar = make_input(tmp_path, data)
    reference = make_reference(tmp_path / "reference.hdf5", series=(0, 1))

    with pytest.raises(converter.ConversionError, match="exactly one image series"):
        converter.convert(nifti, reference, tmp_path / "out.hdf5", sidecar)


def test_output_overwrite_protection(tmp_path):
    data = np.zeros((2, 2, 1, 1), dtype=np.int16)
    nifti, sidecar = make_input(tmp_path, data)
    reference = make_reference(tmp_path / "reference.hdf5")
    output = tmp_path / "out.hdf5"
    output.write_bytes(b"do not replace")

    with pytest.raises(converter.ConversionError, match="--overwrite"):
        converter.convert(nifti, reference, output, sidecar)
    assert output.read_bytes() == b"do not replace"
    converter.convert(nifti, reference, output, sidecar, overwrite=True)
    assert h5py.is_hdf5(output)


def test_missing_and_invalid_inputs_and_json_autodiscovery(tmp_path):
    reference = make_reference(tmp_path / "reference.hdf5")
    with pytest.raises(converter.ConversionError, match="does not exist"):
        converter.convert(tmp_path / "missing.nii.gz", reference, tmp_path / "out.hdf5")

    data = np.zeros((2, 2, 1, 1), dtype=np.int16)
    nifti, sidecar = make_input(tmp_path, data)
    sidecar.write_text("not json")
    with pytest.raises(converter.ConversionError, match="Invalid JSON"):
        converter.convert(nifti, reference, tmp_path / "out.hdf5")

    sidecar.write_text("{}")
    result = converter.convert(nifti, reference, tmp_path / "out.hdf5")
    assert result.image_count == 1

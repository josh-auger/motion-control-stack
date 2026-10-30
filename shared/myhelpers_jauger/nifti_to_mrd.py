#!/usr/bin/env python3
"""Convert a 4D NIfTI image to the raw FIRE-server MRD HDF5 layout.

FIRE stores the scanner wire payload without using ismrmrd.Dataset: each image
is represented by ``image_N/{header,attribute,data}``, with the config and XML
header at the HDF5 root.  This module deliberately writes that same layout.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import uuid
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import h5py
import ismrmrd
import nibabel as nib
import numpy as np


LOGGER = logging.getLogger("nifti_to_mrd")
IMAGE_HEADER_SIZE = 198
MRD_NAMESPACE = "http://www.ismrm.org/ISMRMRD"
LPS_FROM_RAS = np.diag([-1.0, -1.0, 1.0])
DEFAULT_CONFIG_FILE = "fire_emitter_dcm"
REQUIRED_DIRECTION_ATTRIBUTES = (
    "ImageRowDir",
    "ImageColumnDir",
    "ImageSliceNormDir",
)


class ConversionError(RuntimeError):
    """Raised when an input or output cannot satisfy conversion assumptions."""


@dataclass(frozen=True)
class ReferenceInfo:
    config: bytes
    xml: bytes | None
    template_header: ismrmrd.ImageHeader
    image_series_index: int
    image_type: int
    data_dtype: np.dtype
    image_count: int


@dataclass(frozen=True)
class Geometry:
    directions: np.ndarray
    spacing: np.ndarray
    slice_thickness: float
    fov: np.ndarray
    affine_lps: np.ndarray


@dataclass(frozen=True)
class ConversionResult:
    output: Path
    image_count: int
    shape: tuple[int, int, int, int]
    dtype: np.dtype
    slice_order: tuple[int, ...]
    max_abs_error: float


def _bytes(dataset: h5py.Dataset) -> bytes:
    return np.asarray(dataset, dtype=np.uint8).tobytes()


def _without_nul(payload: bytes) -> bytes:
    return payload.split(bytes((0,)), 1)[0]


def _image_number(name: str) -> int:
    match = re.fullmatch(r"image_(\d+)", name)
    if match is None:
        raise ValueError(name)
    return int(match.group(1))


def image_group_names(h5file: h5py.File) -> list[str]:
    names = [name for name in h5file if re.fullmatch(r"image_\d+", name)]
    return sorted(names, key=_image_number)


def discover_json(nifti_path: Path) -> Path:
    name = nifti_path.name
    if name.lower().endswith(".nii.gz"):
        return nifti_path.with_name(name[:-7] + ".json")
    if name.lower().endswith(".nii"):
        return nifti_path.with_name(name[:-4] + ".json")
    raise ConversionError(f"NIfTI filename must end in .nii or .nii.gz: {nifti_path}")


def _require_readable_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise ConversionError(f"{label} does not exist or is not a file: {path}")
    if not os.access(path, os.R_OK):
        raise ConversionError(f"{label} is not readable: {path}")


def load_json(path: Path) -> dict[str, Any]:
    _require_readable_file(path, "JSON sidecar")
    try:
        with path.open("r", encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise ConversionError(f"Invalid JSON sidecar {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ConversionError(f"JSON sidecar must contain an object: {path}")
    return value


def load_nifti(path: Path) -> nib.spatialimages.SpatialImage:
    _require_readable_file(path, "NIfTI input")
    try:
        image = nib.load(str(path))
    except Exception as exc:
        raise ConversionError(f"Cannot read NIfTI input {path}: {exc}") from exc
    if len(image.shape) != 4:
        raise ConversionError(
            f"Expected a 4D NIfTI with shape (X,Y,Z,T), got {image.shape}"
        )
    if any(int(size) < 1 for size in image.shape):
        raise ConversionError(f"NIfTI dimensions must all be positive: {image.shape}")
    affine = np.asarray(image.affine, dtype=np.float64)
    if affine.shape != (4, 4) or not np.all(np.isfinite(affine)):
        raise ConversionError("NIfTI affine must be a finite 4x4 matrix")
    return image


def inspect_reference(path: Path) -> ReferenceInfo:
    """Validate a FIRE raw-layout reference and return reusable conventions."""
    _require_readable_file(path, "Reference MRD")
    try:
        with h5py.File(path, "r") as h5file:
            for required in ("Config File", "Metadata XML"):
                if required not in h5file or not isinstance(h5file[required], h5py.Dataset):
                    raise ConversionError(f"Reference MRD is missing root dataset {required!r}")
            names = image_group_names(h5file)
            if not names:
                raise ConversionError("Reference MRD has no image_N groups")

            xml = _without_nul(_bytes(h5file["Metadata XML"]))
            try:
                ET.fromstring(xml)
                ismrmrd.xsd.CreateFromDocument(xml.decode("utf-8"))
            except Exception as exc:
                raise ConversionError(f"Reference MRD metadata is not valid MRD XML: {exc}") from exc

            series: set[int] = set()
            template: ismrmrd.ImageHeader | None = None
            data_dtype: np.dtype | None = None
            image_type = int(ismrmrd.IMTYPE_MAGNITUDE)
            for name in names:
                group = h5file[name]
                if not isinstance(group, h5py.Group):
                    raise ConversionError(f"Reference object {name} is not an HDF5 group")
                if any(field not in group for field in ("header", "attribute", "data")):
                    raise ConversionError(f"Reference group {name} is incomplete")
                header_bytes = _bytes(group["header"])
                if len(header_bytes) != IMAGE_HEADER_SIZE:
                    raise ConversionError(
                        f"Reference {name}/header is {len(header_bytes)} bytes; expected 198"
                    )
                header = ismrmrd.ImageHeader.from_buffer_copy(header_bytes)
                dtype = np.dtype(ismrmrd.get_dtype_from_data_type(header.data_type))
                expected = (
                    int(header.channels)
                    * int(header.matrix_size[0])
                    * int(header.matrix_size[1])
                    * int(header.matrix_size[2])
                    * dtype.itemsize
                )
                if len(group["data"]) != expected:
                    raise ConversionError(
                        f"Reference {name}/data has {len(group['data'])} bytes; expected {expected}"
                    )
                series.add(int(header.image_series_index))
                if template is None:
                    template = header
                    data_dtype = dtype
                    image_type = int(header.image_type)
            if len(series) != 1:
                raise ConversionError(
                    "Reference MRD must contain exactly one image series; "
                    f"found {sorted(series)}"
                )
            assert template is not None and data_dtype is not None
            return ReferenceInfo(
                config=_bytes(h5file["Config File"]),
                xml=xml,
                template_header=template,
                image_series_index=next(iter(series)),
                image_type=image_type,
                data_dtype=data_dtype,
                image_count=len(names),
            )
    except OSError as exc:
        raise ConversionError(f"Cannot read reference MRD {path}: {exc}") from exc


def default_output_conventions() -> ReferenceInfo:
    """Return neutral FIRE conventions for reference-free conversion."""
    header = ismrmrd.ImageHeader()
    header.version = 1
    header.flags = 15
    header.image_type = ismrmrd.IMTYPE_MAGNITUDE
    header.image_series_index = 0
    config = DEFAULT_CONFIG_FILE.encode("utf-8")[:1024].ljust(1024, bytes((0,)))
    return ReferenceInfo(
        config=config,
        xml=None,
        template_header=header,
        image_series_index=0,
        image_type=int(ismrmrd.IMTYPE_MAGNITUDE),
        data_dtype=np.dtype(np.float32),
        image_count=0,
    )


def derive_geometry(affine_ras: np.ndarray, shape: Sequence[int], metadata: dict[str, Any]) -> Geometry:
    """Translate NIfTI voxel-center geometry from RAS+ coordinates to MRD LPS."""
    affine_lps = np.eye(4, dtype=np.float64)
    affine_lps[:3, :] = LPS_FROM_RAS @ np.asarray(affine_ras, dtype=np.float64)[:3, :]
    voxel_vectors = affine_lps[:3, :3]
    spacing = np.linalg.norm(voxel_vectors, axis=0)
    if not np.all(np.isfinite(spacing)) or np.any(spacing <= 0):
        raise ConversionError(f"NIfTI affine has invalid voxel spacing: {spacing}")
    directions = voxel_vectors / spacing
    gram = directions.T @ directions
    if not np.allclose(gram, np.eye(3), atol=1e-4, rtol=0.0):
        raise ConversionError(
            "NIfTI affine contains shear; MRD direction-cosine headers cannot represent it "
            "without resampling"
        )
    thickness_value = metadata.get("SliceThickness", spacing[2])
    try:
        slice_thickness = float(thickness_value)
    except (TypeError, ValueError) as exc:
        raise ConversionError(f"Invalid SliceThickness value: {thickness_value!r}") from exc
    if not np.isfinite(slice_thickness) or slice_thickness <= 0:
        raise ConversionError(f"SliceThickness must be a positive finite number: {slice_thickness}")
    x, y, z = (int(shape[index]) for index in range(3))
    fov = np.array((spacing[0] * x, spacing[1] * y, spacing[2] * z), dtype=np.float64)
    return Geometry(directions, spacing, slice_thickness, fov, affine_lps)


def slice_order_from_metadata(metadata: dict[str, Any], z_size: int) -> tuple[int, ...]:
    values = metadata.get("SliceTiming")
    if values is None:
        return tuple(range(z_size))
    if not isinstance(values, list) or len(values) != z_size:
        raise ConversionError(
            f"SliceTiming must contain exactly {z_size} entries; got "
            f"{len(values) if isinstance(values, list) else type(values).__name__}"
        )
    try:
        timings = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ConversionError("SliceTiming entries must be numeric") from exc
    if not np.all(np.isfinite(timings)) or np.any(timings < 0):
        raise ConversionError("SliceTiming entries must be finite and nonnegative")
    return tuple(sorted(range(z_size), key=lambda z: (float(timings[z]), z)))


def _scaled_dtype(image: nib.spatialimages.SpatialImage) -> np.dtype:
    # ArrayProxy reports the dtype that results after applying NIfTI slope/intercept.
    proxy = image.dataobj
    sample = np.asanyarray(proxy[0:1, 0:1, 0:1, 0:1])
    return np.dtype(sample.dtype).newbyteorder("=")


def _supported_mrd_dtype(dtype: np.dtype) -> np.dtype:
    dtype = np.dtype(dtype).newbyteorder("=")
    direct = {
        np.dtype(np.uint16), np.dtype(np.int16), np.dtype(np.uint32),
        np.dtype(np.int32), np.dtype(np.float32), np.dtype(np.float64),
        np.dtype(np.complex64), np.dtype(np.complex128),
    }
    if dtype in direct:
        return dtype
    if dtype.kind in "bu" and dtype.itemsize <= 2:
        return np.dtype(np.uint16)
    if dtype.kind == "i" and dtype.itemsize <= 2:
        return np.dtype(np.int16)
    raise ConversionError(f"NIfTI scaled dtype {dtype} has no lossless MRD representation")


def choose_output_dtype(source_dtype: np.dtype, reference_dtype: np.dtype) -> np.dtype:
    """Prefer the reference datatype only when NumPy considers the cast lossless."""
    source_dtype = np.dtype(source_dtype).newbyteorder("=")
    reference_dtype = np.dtype(reference_dtype).newbyteorder("=")
    supported_reference = {
        np.dtype(np.uint16), np.dtype(np.int16), np.dtype(np.uint32),
        np.dtype(np.int32), np.dtype(np.float32), np.dtype(np.float64),
        np.dtype(np.complex64), np.dtype(np.complex128),
    }
    if reference_dtype in supported_reference and np.can_cast(
        source_dtype, reference_dtype, casting="safe"
    ):
        return reference_dtype
    return _supported_mrd_dtype(source_dtype)


def _mrd_data_type(dtype: np.dtype) -> int:
    mapping = {
        np.dtype(np.uint16): ismrmrd.DATATYPE_USHORT,
        np.dtype(np.int16): ismrmrd.DATATYPE_SHORT,
        np.dtype(np.uint32): ismrmrd.DATATYPE_UINT,
        np.dtype(np.int32): ismrmrd.DATATYPE_INT,
        np.dtype(np.float32): ismrmrd.DATATYPE_FLOAT,
        np.dtype(np.float64): ismrmrd.DATATYPE_DOUBLE,
        np.dtype(np.complex64): ismrmrd.DATATYPE_CXFLOAT,
        np.dtype(np.complex128): ismrmrd.DATATYPE_CXDOUBLE,
    }
    try:
        return int(mapping[np.dtype(dtype).newbyteorder("=")])
    except KeyError as exc:
        raise ConversionError(f"Unsupported output dtype: {dtype}") from exc


def _namespace(root: ET.Element) -> str:
    if root.tag.startswith("{"):
        return root.tag[1:].split("}", 1)[0]
    return ""


def _q(namespace: str, name: str) -> str:
    return f"{{{namespace}}}{name}" if namespace else name


def _set_existing(root: ET.Element, namespace: str, path: str, value: Any) -> None:
    element = root.find("/".join(_q(namespace, part) for part in path.split("/")))
    if element is not None:
        element.text = str(value)


def _set_limit(root: ET.Element, namespace: str, name: str, maximum: int) -> None:
    base = f"encoding/encodingLimits/{name}"
    _set_existing(root, namespace, f"{base}/minimum", 0)
    _set_existing(root, namespace, f"{base}/maximum", maximum)
    _set_existing(root, namespace, f"{base}/center", 0)


def _optional_float(
    metadata: dict[str, Any],
    key: str,
    default: float | None = None,
    *,
    positive: bool = False,
) -> float | None:
    value = metadata.get(key)
    if value is None:
        return default
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ConversionError(f"JSON field {key!r} must be numeric") from exc
    if not np.isfinite(number) or (positive and number <= 0):
        qualifier = "positive and finite" if positive else "finite"
        raise ConversionError(f"JSON field {key!r} must be {qualifier}")
    return number


def _multiband_factor(metadata: dict[str, Any], z_size: int) -> int:
    value = metadata.get("MultibandAccelerationFactor")
    if value is not None:
        try:
            factor = int(value)
        except (TypeError, ValueError) as exc:
            raise ConversionError(
                "JSON field 'MultibandAccelerationFactor' must be a positive integer"
            ) from exc
        if factor < 1 or factor > z_size or factor != float(value):
            raise ConversionError(
                "JSON field 'MultibandAccelerationFactor' must be a positive integer"
            )
        return factor

    timings = metadata.get("SliceTiming")
    if isinstance(timings, list) and len(timings) == z_size:
        numeric = np.asarray(timings, dtype=np.float64)
        _unique, counts = np.unique(numeric, return_counts=True)
        if len(counts) and np.all(counts == counts[0]):
            return int(counts[0])
    return 1


def build_reference_free_xml(
    shape: Sequence[int],
    geometry: Geometry,
    metadata: dict[str, Any],
    nifti_name: str,
) -> bytes:
    """Construct a schema-backed MRD XML header without scanner metadata."""
    x, y, z, t = (int(value) for value in shape)
    source_name = (
        nifti_name[:-7]
        if nifti_name.lower().endswith(".nii.gz")
        else Path(nifti_name).stem
    )
    protocol = str(
        metadata.get("ProtocolName")
        or metadata.get("SeriesDescription")
        or source_name
        or "offline_nifti"
    )
    measurement_id = f"offline-nifti-{uuid.uuid5(uuid.NAMESPACE_URL, f'{source_name}|{tuple(shape)}')}"

    matrix = ismrmrd.xsd.matrixSizeType(x=x, y=y, z=1)
    fov = ismrmrd.xsd.fieldOfViewMm(
        x=float(geometry.fov[0]),
        y=float(geometry.fov[1]),
        z=float(geometry.slice_thickness),
    )
    encoded_space = ismrmrd.xsd.encodingSpaceType(
        matrixSize=matrix, fieldOfView_mm=fov
    )
    recon_space = ismrmrd.xsd.encodingSpaceType(
        matrixSize=ismrmrd.xsd.matrixSizeType(x=x, y=y, z=1),
        fieldOfView_mm=ismrmrd.xsd.fieldOfViewMm(
            x=float(geometry.fov[0]),
            y=float(geometry.fov[1]),
            z=float(geometry.slice_thickness),
        ),
    )

    def zero_limit() -> ismrmrd.xsd.limitType:
        return ismrmrd.xsd.limitType(minimum=0, maximum=0, center=0)

    limits = ismrmrd.xsd.encodingLimitsType(
        kspace_encoding_step_1=ismrmrd.xsd.limitType(
            minimum=0, maximum=y - 1, center=y // 2
        ),
        kspace_encoding_step_2=zero_limit(),
        average=zero_limit(),
        slice=ismrmrd.xsd.limitType(minimum=0, maximum=z - 1, center=0),
        contrast=zero_limit(),
        phase=zero_limit(),
        repetition=ismrmrd.xsd.limitType(minimum=0, maximum=t - 1, center=0),
        set=zero_limit(),
        segment=zero_limit(),
    )
    echo_train_length = metadata.get("EchoTrainLength")
    if echo_train_length is not None:
        original_echo_train_length = echo_train_length
        try:
            echo_train_length = int(original_echo_train_length)
        except (TypeError, ValueError) as exc:
            raise ConversionError("JSON field 'EchoTrainLength' must be an integer") from exc
        if echo_train_length < 1 or echo_train_length != float(original_echo_train_length):
            raise ConversionError("JSON field 'EchoTrainLength' must be a positive integer")
    encoding = ismrmrd.xsd.encodingType(
        encodedSpace=encoded_space,
        reconSpace=recon_space,
        encodingLimits=limits,
        trajectory=ismrmrd.xsd.trajectoryType.CARTESIAN,
        echoTrainLength=echo_train_length,
    )

    tr = _optional_float(metadata, "RepetitionTime", positive=True)
    te = _optional_float(metadata, "EchoTime", positive=True)
    flip = _optional_float(metadata, "FlipAngle")
    echo_spacing = _optional_float(metadata, "EffectiveEchoSpacing", positive=True)
    sequence_type = metadata.get("PulseSequenceName") or metadata.get("ScanningSequence")
    sequence = ismrmrd.xsd.sequenceParametersType(
        TR=[] if tr is None else [tr * 1000.0],
        TE=[] if te is None else [te * 1000.0],
        flipAngle_deg=[] if flip is None else [flip],
        sequence_type=None if sequence_type is None else str(sequence_type),
        echo_spacing=[] if echo_spacing is None else [echo_spacing * 1000.0],
    )

    multiband_factor = _multiband_factor(metadata, z)
    spacing_between_slices = _optional_float(
        metadata,
        "SpacingBetweenSlices",
        default=float(geometry.spacing[2]),
        positive=True,
    )
    strings = [
        ismrmrd.xsd.userParameterStringType(
            name="NIfTIToMRDConverter",
            value="shared/myhelpers_jauger/nifti_to_mrd.py",
        ),
        ismrmrd.xsd.userParameterStringType(name="SourceNIfTI", value=nifti_name),
    ]
    if metadata.get("PhaseEncodingDirection") is not None:
        strings.append(
            ismrmrd.xsd.userParameterStringType(
                name="NIfTIPhaseEncodingDirection",
                value=str(metadata["PhaseEncodingDirection"]),
            )
        )
    user_parameters = ismrmrd.xsd.userParametersType(
        userParameterLong=[
            ismrmrd.xsd.userParameterLongType(
                name="MultiBandFactor", value=multiband_factor
            )
        ],
        userParameterDouble=[
            ismrmrd.xsd.userParameterDoubleType(
                name="SpacingBetweenSlices", value=spacing_between_slices
            )
        ],
        userParameterString=strings,
    )

    field_strength = _optional_float(
        metadata, "MagneticFieldStrength", default=0.0, positive=True
    )
    system = ismrmrd.xsd.acquisitionSystemInformationType(
        systemVendor=str(metadata.get("Manufacturer") or "Unknown"),
        systemModel=str(metadata.get("ManufacturersModelName") or "Unknown"),
        systemFieldStrength_T=field_strength,
    )
    patient_position = None
    if metadata.get("PatientPosition") is not None:
        try:
            patient_position = ismrmrd.xsd.patientPositionType(
                str(metadata["PatientPosition"])
            )
        except ValueError:
            LOGGER.warning(
                "Ignoring unsupported PatientPosition value %r",
                metadata["PatientPosition"],
            )
    measurement = ismrmrd.xsd.measurementInformationType(
        measurementID=measurement_id,
        patientPosition=patient_position,
        protocolName=protocol,
        sequenceName=(
            None
            if metadata.get("PulseSequenceName") is None
            else str(metadata["PulseSequenceName"])
        ),
        seriesDescription=(
            None
            if metadata.get("SeriesDescription") is None
            else str(metadata["SeriesDescription"])
        ),
    )
    frequency = _optional_float(metadata, "ImagingFrequency", positive=True)
    experimental = None
    if frequency is not None:
        experimental = ismrmrd.xsd.experimentalConditionsType(
            H1resonanceFrequency_Hz=int(round(frequency * 1_000_000.0))
        )

    header = ismrmrd.xsd.ismrmrdHeader(
        version=1,
        measurementInformation=measurement,
        acquisitionSystemInformation=system,
        experimentalConditions=experimental,
        encoding=[encoding],
        sequenceParameters=sequence,
        userParameters=user_parameters,
    )
    xml = ismrmrd.xsd.ToXML(header).encode("utf-8")
    try:
        ismrmrd.xsd.CreateFromDocument(xml.decode("utf-8"))
    except Exception as exc:
        raise ConversionError(f"Generated MRD XML is invalid: {exc}") from exc
    return xml + bytes((0,))


def build_xml(
    reference_xml: bytes | None,
    shape: Sequence[int],
    geometry: Geometry,
    metadata: dict[str, Any],
    nifti_name: str,
) -> bytes:
    if reference_xml is None:
        return build_reference_free_xml(shape, geometry, metadata, nifti_name)
    try:
        root = ET.fromstring(_without_nul(reference_xml))
    except ET.ParseError as exc:
        raise ConversionError(f"Cannot parse reference MRD XML: {exc}") from exc
    namespace = _namespace(root)
    if namespace:
        ET.register_namespace("", namespace)
    encodings = root.findall(_q(namespace, "encoding"))
    if len(encodings) != 1:
        raise ConversionError(
            f"Reference XML must describe exactly one encoding; found {len(encodings)}"
        )
    x, y, z, t = (int(value) for value in shape)

    for space in ("encodedSpace", "reconSpace"):
        for axis, size in zip("xyz", (x, y, 1)):
            _set_existing(root, namespace, f"encoding/{space}/matrixSize/{axis}", size)
        for axis, value in zip("xyz", (geometry.fov[0], geometry.fov[1], geometry.slice_thickness)):
            _set_existing(root, namespace, f"encoding/{space}/fieldOfView_mm/{axis}", f"{value:.9g}")
    _set_limit(root, namespace, "slice", z - 1)
    _set_limit(root, namespace, "repetition", t - 1)

    protocol = metadata.get("ProtocolName") or metadata.get("SeriesDescription")
    if protocol:
        _set_existing(root, namespace, "measurementInformation/protocolName", protocol)
    generated_id = f"offline-nifti-{uuid.uuid4()}"
    _set_existing(root, namespace, "measurementInformation/measurementID", generated_id)
    _set_existing(root, namespace, "studyInformation/studyInstanceUID", generated_id)
    for key, path in (
        ("Manufacturer", "acquisitionSystemInformation/systemVendor"),
        ("ManufacturersModelName", "acquisitionSystemInformation/systemModel"),
        ("MagneticFieldStrength", "acquisitionSystemInformation/systemFieldStrength_T"),
    ):
        if metadata.get(key) is not None:
            _set_existing(root, namespace, path, metadata[key])
    # ISMRMRD sequence times are milliseconds; BIDS JSON times are seconds.
    for key, xml_name in (("RepetitionTime", "TR"), ("EchoTime", "TE")):
        if metadata.get(key) is not None:
            _set_existing(root, namespace, f"sequenceParameters/{xml_name}", float(metadata[key]) * 1000.0)
    if metadata.get("FlipAngle") is not None:
        _set_existing(root, namespace, "sequenceParameters/flipAngle_deg", metadata["FlipAngle"])

    # Store converter provenance in an existing userParameters block when possible.
    user_parameters = root.find(_q(namespace, "userParameters"))
    if user_parameters is not None:
        for name, value in (
            ("NIfTIToMRDConverter", "shared/myhelpers_jauger/nifti_to_mrd.py"),
            ("SourceNIfTI", nifti_name),
            ("NIfTIPhaseEncodingDirection", metadata.get("PhaseEncodingDirection", "unknown")),
        ):
            parameter = ET.SubElement(user_parameters, _q(namespace, "userParameterString"))
            ET.SubElement(parameter, _q(namespace, "name")).text = str(name)
            ET.SubElement(parameter, _q(namespace, "value")).text = str(value)

    xml = ET.tostring(root, encoding="utf-8", xml_declaration=True)
    try:
        ismrmrd.xsd.CreateFromDocument(xml.decode("utf-8"))
    except Exception as exc:
        raise ConversionError(f"Generated MRD XML is invalid: {exc}") from exc
    return xml + bytes((0,))


def _format_meta_vector(values: Sequence[float]) -> list[str]:
    return [f"{float(value):.18f}" for value in values]


def build_attribute(
    metadata: dict[str, Any], header: ismrmrd.ImageHeader
) -> bytes:
    """Build FIRE-compatible per-image Meta attributes.

    FIRE negates ``ImageSliceNormDir`` when constructing its image direction
    matrix.  Storing the negative header slice direction therefore makes the
    resulting FIRE matrix exactly ``[read_dir, phase_dir, slice_dir]``.
    """
    meta = ismrmrd.Meta()
    meta["DataRole"] = "Image"
    meta["ImageProcessingHistory"] = ["NIFTI_TO_MRD", "PYTHON"]
    meta["Keep_image_geometry"] = 1
    meta["SequenceDescriptionAdditional"] = str(
        metadata.get("SeriesDescription") or metadata.get("ProtocolName") or "NIFTI_TO_MRD"
    )
    if metadata.get("PhaseEncodingDirection") is not None:
        meta["NIfTIPhaseEncodingDirection"] = str(metadata["PhaseEncodingDirection"])
    meta["ImageRowDir"] = _format_meta_vector(header.read_dir)
    meta["ImageColumnDir"] = _format_meta_vector(header.phase_dir)
    meta["ImageSliceNormDir"] = _format_meta_vector(
        [-float(value) for value in header.slice_dir]
    )
    meta["SlicePosLightMarker"] = _format_meta_vector(header.position)
    if metadata.get("RepetitionTime") is not None:
        meta["RepetitionTime"] = f"{float(metadata['RepetitionTime']) * 1000.0:.12g}"
    if metadata.get("EchoTime") is not None:
        meta["EchoTime"] = f"{float(metadata['EchoTime']) * 1000.0:.12g}"
    sequence_description = metadata.get("SeriesDescription") or metadata.get(
        "ProtocolName"
    )
    if sequence_description is not None:
        meta["SequenceDescription"] = str(sequence_description)
    return meta.serialize().encode("utf-8")


def _meta_vector(meta: ismrmrd.Meta, key: str) -> np.ndarray:
    if key not in meta.keys():
        raise ConversionError(f"Image attribute payload is missing required field {key!r}")
    values = meta[key]
    if not isinstance(values, list) or len(values) != 3:
        raise ConversionError(f"Image attribute {key!r} must contain three values")
    try:
        vector = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ConversionError(f"Image attribute {key!r} must be numeric") from exc
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ConversionError(f"Image attribute {key!r} must contain three finite values")
    return vector


def validate_image_attribute(
    attribute: bytes, header: ismrmrd.ImageHeader
) -> ismrmrd.Meta:
    """Decode and validate the image geometry metadata required by FIRE."""
    try:
        payload = _without_nul(attribute).decode("utf-8")
        meta = ismrmrd.Meta.deserialize(payload)
    except Exception as exc:
        raise ConversionError(f"Image attribute payload is not valid ISMRMRD Meta XML: {exc}") from exc

    row = _meta_vector(meta, "ImageRowDir")
    column = _meta_vector(meta, "ImageColumnDir")
    slice_normal = _meta_vector(meta, "ImageSliceNormDir")
    expected = (
        np.asarray(header.read_dir, dtype=np.float64),
        np.asarray(header.phase_dir, dtype=np.float64),
        -np.asarray(header.slice_dir, dtype=np.float64),
    )
    for key, actual, wanted in zip(
        REQUIRED_DIRECTION_ATTRIBUTES, (row, column, slice_normal), expected
    ):
        if not np.allclose(actual, wanted, atol=1e-6, rtol=0.0):
            raise ConversionError(
                f"Image attribute {key!r} is inconsistent with the MRD image header"
            )

    # This is the exact matrix construction performed by FIRE's
    # handleData.get_image_direction_matrix().
    fire_matrix = np.column_stack((row, column, -slice_normal))
    header_matrix = np.column_stack(
        (header.read_dir, header.phase_dir, header.slice_dir)
    )
    if not np.allclose(fire_matrix, header_matrix, atol=1e-6, rtol=0.0):
        raise ConversionError("FIRE image direction matrix does not match the MRD header")
    position = _meta_vector(meta, "SlicePosLightMarker")
    if not np.allclose(position, header.position, atol=1e-4, rtol=0.0):
        raise ConversionError(
            "Image attribute 'SlicePosLightMarker' is inconsistent with the MRD header"
        )
    return meta


def _slice_position(geometry: Geometry, shape: Sequence[int], z_index: int) -> np.ndarray:
    center_index = np.array(
        ((int(shape[0]) - 1) / 2.0, (int(shape[1]) - 1) / 2.0, float(z_index), 1.0)
    )
    return (geometry.affine_lps @ center_index)[:3]


def _relative_timestamp(metadata: dict[str, Any], repetition: int, z_index: int) -> int:
    """Return a derived relative 2.5 ms tick count, or zero if timing is unavailable."""
    tr = metadata.get("RepetitionTime")
    if tr is None:
        return 0
    try:
        seconds = repetition * float(tr)
        timings = metadata.get("SliceTiming")
        if isinstance(timings, list):
            seconds += float(timings[z_index])
        ticks = round(seconds * 400.0)
    except (TypeError, ValueError, IndexError):
        return 0
    return int(max(0, min(ticks, np.iinfo(np.uint32).max)))


def build_header(
    reference: ReferenceInfo,
    shape: Sequence[int],
    geometry: Geometry,
    dtype: np.dtype,
    metadata: dict[str, Any],
    z_index: int,
    repetition: int,
) -> ismrmrd.ImageHeader:
    header = ismrmrd.ImageHeader()
    header.version = int(reference.template_header.version) or 1
    header.flags = int(reference.template_header.flags)
    header.data_type = _mrd_data_type(dtype)
    header.measurement_uid = 0  # Reference scanner UID is not valid for an offline conversion.
    header.matrix_size = (int(shape[0]), int(shape[1]), 1)
    header.field_of_view = (
        float(geometry.fov[0]), float(geometry.fov[1]), float(geometry.slice_thickness)
    )
    header.channels = 1
    position = _slice_position(geometry, shape, z_index)
    header.position = tuple(float(value) for value in position)
    header.read_dir = tuple(float(value) for value in geometry.directions[:, 0])
    header.phase_dir = tuple(float(value) for value in geometry.directions[:, 1])
    header.slice_dir = tuple(float(value) for value in geometry.directions[:, 2])
    header.average = 0
    header.slice = int(z_index)
    header.contrast = 0
    header.phase = 0
    header.repetition = int(repetition)
    header.set = 0
    header.acquisition_time_stamp = _relative_timestamp(metadata, repetition, z_index)
    header.image_type = int(reference.image_type) or int(ismrmrd.IMTYPE_MAGNITUDE)
    header.image_index = 0  # This is how incoming FIRE images are stored in the reference.
    header.image_series_index = int(reference.image_series_index)
    return header


def _write_bytes(group: h5py.Group, name: str, payload: bytes) -> None:
    group.create_dataset(name, data=np.frombuffer(payload, dtype=np.uint8))


def write_mrd(
    output: Path,
    image: nib.spatialimages.SpatialImage,
    source_data: np.ndarray,
    metadata: dict[str, Any],
    reference: ReferenceInfo,
    geometry: Geometry,
    dtype: np.dtype,
    order: Sequence[int],
    nifti_name: str,
) -> None:
    x, y, z_size, t_size = (int(value) for value in image.shape)
    xml = build_xml(reference.xml, image.shape, geometry, metadata, nifti_name)
    little_dtype = np.dtype(dtype).newbyteorder("<")
    with h5py.File(output, "w") as h5file:
        _write_bytes(h5file, "Config File", reference.config)
        _write_bytes(h5file, "Metadata XML", xml)
        image_number = 0
        for repetition in range(t_size):
            for z_index in order:
                header = build_header(
                    reference, image.shape, geometry, dtype, metadata, z_index, repetition
                )
                attribute = build_attribute(metadata, header)
                # NIfTI uses (X,Y); MRD bytes are (channel,z,Y,X), with X fastest.
                source_slice = source_data[:, :, z_index, repetition]
                serialized = np.ascontiguousarray(source_slice.T, dtype=little_dtype)
                group = h5file.create_group(f"image_{image_number}")
                _write_bytes(group, "header", bytes(header))
                _write_bytes(group, "attribute", attribute)
                _write_bytes(group, "data", serialized.tobytes(order="C"))
                image_number += 1
        h5file.flush()
    expected = z_size * t_size
    if image_number != expected:
        raise ConversionError(f"Internal image count error: wrote {image_number}, expected {expected}")


def read_image(group: h5py.Group) -> tuple[ismrmrd.ImageHeader, bytes, np.ndarray]:
    header_bytes = _bytes(group["header"])
    if len(header_bytes) != IMAGE_HEADER_SIZE:
        raise ConversionError(f"Invalid image header length: {len(header_bytes)}")
    header = ismrmrd.ImageHeader.from_buffer_copy(header_bytes)
    dtype = np.dtype(ismrmrd.get_dtype_from_data_type(header.data_type)).newbyteorder("<")
    x, y, z = (int(value) for value in header.matrix_size)
    channels = int(header.channels)
    raw = _bytes(group["data"])
    expected = channels * z * y * x * dtype.itemsize
    if len(raw) != expected:
        raise ConversionError(f"Image data is {len(raw)} bytes; expected {expected}")
    data = np.frombuffer(raw, dtype=dtype).reshape(channels, z, y, x)
    return header, _bytes(group["attribute"]), data


def reconstruct_4d(path: Path, shape: Sequence[int]) -> np.ndarray:
    x, y, z_size, t_size = (int(value) for value in shape)
    result: np.ndarray | None = None
    seen: set[tuple[int, int]] = set()
    with h5py.File(path, "r") as h5file:
        names = image_group_names(h5file)
        if len(names) != z_size * t_size:
            raise ConversionError(f"MRD contains {len(names)} images; expected {z_size * t_size}")
        for name in names:
            header, _attribute, data = read_image(h5file[name])
            if int(header.channels) != 1 or int(header.matrix_size[2]) != 1:
                raise ConversionError(f"{name} is not a single-channel 2D image")
            if tuple(int(v) for v in header.matrix_size[:2]) != (x, y):
                raise ConversionError(f"{name} has incompatible matrix size")
            key = (int(header.repetition), int(header.slice))
            if key in seen or not (0 <= key[0] < t_size and 0 <= key[1] < z_size):
                raise ConversionError(f"Invalid or duplicate repetition/slice index in {name}: {key}")
            if result is None:
                result = np.empty((x, y, z_size, t_size), dtype=data.dtype.newbyteorder("="))
            result[:, :, key[1], key[0]] = data[0, 0].T
            seen.add(key)
    assert result is not None
    return result


def validate_output(
    path: Path,
    image: nib.spatialimages.SpatialImage,
    source_data: np.ndarray,
    metadata: dict[str, Any],
    geometry: Geometry,
    order: Sequence[int],
    dtype: np.dtype,
) -> float:
    expected_count = int(image.shape[2]) * int(image.shape[3])
    max_abs_error = 0.0
    with h5py.File(path, "r") as h5file:
        for required in ("Config File", "Metadata XML"):
            if required not in h5file:
                raise ConversionError(f"Output validation failed: missing {required}")
        xml = _without_nul(_bytes(h5file["Metadata XML"]))
        try:
            ismrmrd.xsd.CreateFromDocument(xml.decode("utf-8"))
        except Exception as exc:
            raise ConversionError(f"Output XML validation failed: {exc}") from exc
        names = image_group_names(h5file)
        if len(names) != expected_count:
            raise ConversionError(
                f"Output has {len(names)} images; expected {expected_count}"
            )
        for image_number, name in enumerate(names):
            repetition = image_number // int(image.shape[2])
            expected_slice = int(order[image_number % int(image.shape[2])])
            header, attribute, data = read_image(h5file[name])
            if (int(header.repetition), int(header.slice)) != (repetition, expected_slice):
                raise ConversionError(
                    f"Output ordering error in {name}: got repetition/slice "
                    f"{header.repetition}/{header.slice}, expected {repetition}/{expected_slice}"
                )
            if np.dtype(ismrmrd.get_dtype_from_data_type(header.data_type)) != np.dtype(dtype):
                raise ConversionError(f"Output datatype mismatch in {name}")
            expected_position = _slice_position(geometry, image.shape, expected_slice)
            if not np.allclose(header.position, expected_position, atol=1e-4, rtol=0.0):
                raise ConversionError(f"Output geometry mismatch in {name}")
            header_directions = np.column_stack(
                (header.read_dir, header.phase_dir, header.slice_dir)
            )
            if not np.allclose(
                header_directions, geometry.directions, atol=1e-6, rtol=0.0
            ):
                raise ConversionError(f"Output direction-cosine mismatch in {name}")
            expected_fov = (
                float(geometry.fov[0]),
                float(geometry.fov[1]),
                float(geometry.slice_thickness),
            )
            if not np.allclose(
                header.field_of_view, expected_fov, atol=1e-4, rtol=0.0
            ):
                raise ConversionError(f"Output field-of-view mismatch in {name}")
            validate_image_attribute(attribute, header)
            source = source_data[:, :, expected_slice, repetition]
            restored = data[0, 0].T
            if not np.array_equal(restored, source):
                error = float(np.max(np.abs(restored.astype(np.complex128) - source.astype(np.complex128))))
                max_abs_error = max(max_abs_error, error)
                raise ConversionError(
                    f"Output intensity mismatch in {name}; maximum absolute error {error}"
                )
    return max_abs_error


def convert(
    nifti_path: Path | str,
    reference_path: Path | str | None = None,
    output_path: Path | str | None = None,
    json_path: Path | str | None = None,
    overwrite: bool = False,
) -> ConversionResult:
    # Preserve the original positional API: convert(nifti, reference, output).
    # Reference-free callers may use convert(nifti, output) or output_path=...
    if output_path is None:
        if reference_path is None:
            raise ConversionError("An output path is required")
        output_path = reference_path
        reference_path = None
    nifti_path = Path(nifti_path).expanduser().resolve()
    reference_path = (
        Path(reference_path).expanduser().resolve()
        if reference_path is not None
        else None
    )
    output_path = Path(output_path).expanduser().resolve()
    json_path = (
        Path(json_path).expanduser().resolve() if json_path is not None else discover_json(nifti_path)
    )
    if output_path.exists() and not overwrite:
        raise ConversionError(f"Output already exists (use --overwrite): {output_path}")
    if output_path.is_dir():
        raise ConversionError(f"Output path is a directory: {output_path}")
    if not output_path.parent.is_dir():
        raise ConversionError(f"Output directory does not exist: {output_path.parent}")
    if output_path in tuple(
        path for path in (nifti_path, reference_path, json_path) if path is not None
    ):
        raise ConversionError("Output path must differ from every input path")

    LOGGER.info("Reading NIfTI and metadata")
    image = load_nifti(nifti_path)
    metadata = load_json(json_path)
    if reference_path is None:
        LOGGER.info("Using reference-free FIRE metadata conventions")
        reference = default_output_conventions()
    else:
        LOGGER.info("Inspecting FIRE reference MRD")
        reference = inspect_reference(reference_path)
    geometry = derive_geometry(image.affine, image.shape, metadata)
    order = slice_order_from_metadata(metadata, int(image.shape[2]))
    source_dtype = _scaled_dtype(image)
    dtype = choose_output_dtype(source_dtype, reference.data_dtype)
    # Materialize once. Random slicing of a compressed .nii.gz proxy repeatedly
    # decompresses the file and is pathologically slow for long time series.
    source_data = np.asanyarray(image.dataobj)
    LOGGER.info(
        "Converting shape %s as %s (%d images)",
        image.shape,
        dtype,
        int(image.shape[2]) * int(image.shape[3]),
    )

    temporary = output_path.with_name(f".{output_path.name}.{uuid.uuid4().hex}.tmp")
    try:
        write_mrd(
            temporary, image, source_data, metadata, reference, geometry, dtype, order, nifti_path.name
        )
        LOGGER.info("Validating structure, geometry, and voxel values")
        max_abs_error = validate_output(
            temporary, image, source_data, metadata, geometry, order, dtype
        )
        os.replace(temporary, output_path)
    except BaseException:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            LOGGER.warning("Could not remove incomplete temporary output %s", temporary)
        raise
    return ConversionResult(
        output=output_path,
        image_count=int(image.shape[2]) * int(image.shape[3]),
        shape=tuple(int(value) for value in image.shape),
        dtype=dtype,
        slice_order=tuple(order),
        max_abs_error=max_abs_error,
    )


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert a 4D NIfTI and JSON sidecar to FIRE raw-layout MRD HDF5"
    )
    parser.add_argument("--nifti", required=True, type=Path, help="Input .nii or .nii.gz file")
    parser.add_argument(
        "--json", dest="json_path", type=Path,
        help="JSON sidecar (default: same basename as --nifti)",
    )
    parser.add_argument(
        "--reference-mrd",
        type=Path,
        help="Optional FIRE reference .hdf5 for legacy template-assisted conversion",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output .hdf5 path")
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing output")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = build_argument_parser().parse_args(argv)
    try:
        result = convert(
            args.nifti,
            reference_path=args.reference_mrd,
            output_path=args.output,
            json_path=args.json_path,
            overwrite=args.overwrite,
        )
    except (ConversionError, OSError, ValueError) as exc:
        LOGGER.error("%s", exc)
        return 2
    LOGGER.info(
        "Wrote %s: %d images, shape %s, dtype %s, max round-trip error %.3g",
        result.output,
        result.image_count,
        result.shape,
        result.dtype,
        result.max_abs_error,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

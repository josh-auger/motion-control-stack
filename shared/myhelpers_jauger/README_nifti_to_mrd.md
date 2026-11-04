# NIfTI-to-MRD converter

`nifti_to_mrd.py` converts one 4D NIfTI scan and its BIDS-style JSON
sidecar into the raw HDF5 layout written by `apps/python_fire_server`.
It is intended for offline analysis and future FIRE replay work; no FIRE
service needs to be running.

## Dependencies

The helper uses the repository's existing Python dependencies: `numpy`,
`nibabel`, `h5py`, and `ismrmrd`.

## Usage

Reference-free conversion is the default:

```bash
python shared/myhelpers_jauger/nifti_to_mrd.py \
  --nifti /path/to/scan.nii.gz \
  --json /path/to/scan.json \
  --output /path/to/converted.hdf5
```

`--json` may be omitted when the sidecar has the NIfTI basename. Existing
outputs are protected unless `--overwrite` is supplied.

The legacy reference-assisted path remains available:

```bash
python shared/myhelpers_jauger/nifti_to_mrd.py \
  --nifti /path/to/scan.nii.gz \
  --json /path/to/scan.json \
  --reference-mrd /path/to/FIRE-reference.hdf5 \
  --output /path/to/converted.hdf5
```

An optional reference must use the FIRE raw layout, contain one image series,
and provide root `Config File` and `Metadata XML` datasets plus complete
`image_N/{header,attribute,data}` groups.

The October 6 validation command is:

```bash
python shared/myhelpers_jauger/nifti_to_mrd.py \
  --nifti /home/jauger/Radiology_Research/Scan_data/20261006_MRN6045885_handmotor_1_baseline/fromResearchSynapse_DICOMS_20261006/6045885-3-converted/6045885/20261006/joshauger_motioncorrection_handmotor_sesh1_baseline/18/func-bold_task-rest480_run-01_slimmon/func-bold_task-rest480_run-01_slimmon_func-bold_task-rest480_run-01_SLIMMON_20261006143118_18.nii.gz \
  --output /home/jauger/Radiology_Research/Scan_data/20261006_MRN6045885_handmotor_1_baseline/dataConversionTesting_Nifti_to_MRD_rest480_run02/nifti_to_mrd_oct6_reference_free_v3.hdf5
```

## Input metadata

The sidecar must be a JSON object, but no individual JSON key is mandatory:
matrix size, voxel geometry, and datatype come from NIfTI. `ProtocolName` or
`SeriesDescription` is preferred; the NIfTI basename is the fallback protocol.
Useful optional fields include `RepetitionTime`, `EchoTime`, `FlipAngle`,
`SliceTiming`, `SliceThickness`, `SpacingBetweenSlices`,
`MultibandAccelerationFactor`, `MagneticFieldStrength`, `Manufacturer`,
`ManufacturersModelName`, `PatientPosition`, `PulseSequenceName`,
`EffectiveEchoSpacing`, `EchoTrainLength`, and `PhaseEncodingDirection`.
Provided numeric values must be finite and physically valid. Missing slice
timing produces ascending-Z ordering; missing multiband information defaults to
one unless a uniform factor can be derived from `SliceTiming`.

## Data, ordering, and geometry

The NIfTI is interpreted as `(X,Y,Z,T)`. Output images are time-major: every
slice of volume 0, then every slice of volume 1. Within a volume, valid JSON
`SliceTiming` values define acquisition order (stable by NIfTI Z index for
simultaneous slices); otherwise Z order is used. MRD `repetition` is `T`,
MRD `slice` is the NIfTI Z index, and `image_index` remains zero to match the
observed FIRE reference convention.

Each NIfTI `(X,Y)` slice is transposed only for MRD's byte layout
`(channel,z,Y,X)` and is transposed back on read. This is serialization, not
resampling. NIfTI slope/intercept scaling is applied by nibabel. The reference
datatype is used only when NumPy classifies that cast as lossless; otherwise a
lossless MRD datatype is selected. Conversion validation compares every stored
voxel to nibabel's scaled values exactly.

The affine's three voxel vectors are converted from NIfTI RAS+ to MRD LPS.
Their normalized vectors populate `read_dir`, `phase_dir`, and `slice_dir`.
Each image `position` uses the Siemens in-plane field-of-view midpoint
convention. If the NIfTI affine in LPS has first-voxel center `O` and oriented
voxel-step columns `v_x`, `v_y`, and `v_z`, the position for slice `z` is:

```text
P_MRD(z) = O + (N_x / 2) v_x + (N_y / 2) v_y + z v_z
```

Only the two in-plane dimensions use `N/2`. The through-plane coordinate remains
the individual NIfTI slice index, so the converter does not introduce a
half-slice displacement. In-plane FOV remains matrix size times affine spacing,
and header FOV Z remains JSON `SliceThickness` (or affine Z spacing).

This position convention replaces the former `(N-1)/2` calculation and the
empirical `--inplane-origin-correction-voxels 0.5 -0.5` workflow. That CLI
option has been removed: no additional correction should be supplied. Applying
the revised FIRE `N/2` interpretation to the generated position recovers the
original NIfTI first-voxel center exactly:

```text
O = P_MRD(z) - (N_x / 2) v_x - (N_y / 2) v_y - z v_z
```

The calculation uses affine voxel-step vectors directly, so permutations,
reflections (including the observed DICOM-derived Y reflection), and oblique
orthogonal orientations are handled without hard-coded correction signs or
resampling. Sheared affines remain unsupported because MRD direction cosines
cannot represent shear exactly.

This convention is intended to match the native Siemens MRD stream and the
revised FIRE origin calculation. It assumes one 2-D MRD image per NIfTI slice;
it should not be generalized to a centered 3-D MRD record without separately
defining the through-plane convention.

In reference-free mode the MRD XML is constructed with the installed ISMRMRD
schema classes. It contains reconstructed/encoded matrix and FOV, Cartesian
encoding, slice and repetition limits, protocol/system descriptions, available
sequence parameters, multiband factor, spacing between slices, and converter
provenance. BIDS times in seconds are converted to MRD milliseconds. The root
`Config File` is the 1,024-byte FIRE application default `fire_emitter_dcm`.
The optional legacy path instead updates the supplied reference XML and copies
its config payload.

Image attributes are generated directly and never copy scanner miniheaders.
When TR exists, `acquisition_time_stamp` is a relative value derived in standard
2.5 ms ticks from TR and SliceTiming; it is not a wall-clock acquisition time.
Measurement identifiers are explicitly synthetic, and image `measurement_uid`
is zero.

Every image attribute payload includes `ImageRowDir`, `ImageColumnDir`, and
`ImageSliceNormDir`, which FIRE requires when establishing its reference image
coordinate frame. Row and column directions match the MRD header `read_dir`
and `phase_dir`; `ImageSliceNormDir` is the negative of header `slice_dir`
because FIRE negates that metadata field when building its direction matrix.
`SlicePosLightMarker` matches header `position`. Repetition time, echo time,
and sequence description are populated from the JSON sidecar when available.
All vectors are serialized as three finite numeric strings using ISMRMRD Meta
XML, decoded again during output validation, and checked against each header.

## FIRE reference findings

`python_fire_server` reads the scanner's image message and writes the received
198-byte header, attribute payload, and data payload unchanged. It does not
alter incoming image pixels or headers. The October 6 reference has:

- `Config File` (1,024 `uint8` bytes) and `Metadata XML` (219,236 `uint8`
  bytes) at the root;
- 2,800 `image_N` groups, exactly 50 volumes x 56 slices, in one series;
- `header` (198 bytes), XML `attribute` (about 5 KB), and `data` (32,400
  bytes) in every image group;
- 90 x 90 x 1, one-channel float32 magnitude images;
- time-major ordering with SMS/acquisition-timing slice order;
- `image_index=0`, `image_series_index=0`, `slice=0..55`, and
  `repetition=0..49`;
- 216 x 216 x 2.4 mm per-image FOV and LPS scanner geometry.

The matching NIfTI is 90 x 90 x 56 x 50, int16, with 2.4 mm isotropic
spacing, TR 1 s, TE 29 ms, and 56 SliceTiming entries. Reference pixels and
the NIfTI are extremely highly correlated but are not numerically identical:
the reference contains scanner float reconstructions while the NIfTI contains
DICOM-exported integer values. The converter therefore validates against the
NIfTI, rather than modifying values to resemble the reference.

## Validation and limitations

Every conversion reopens its temporary output, schema-parses the MRD XML,
checks all image records, indices, ordering, datatype, positions, and exact
voxel values, then atomically publishes the requested output. Tests can be run
with:

```bash
python -m pytest -q shared/myhelpers_jauger/tests/test_nifti_to_mrd.py
```

Known limitations:

- only single-channel, 2D-per-record magnitude image series are produced;
- a reference with multiple image series or encodings is rejected;
- Siemens private scanner slice numbers and absolute acquisition timestamps
  are not recoverable from NIfTI and are not fabricated;
- private ICE miniheaders, coil descriptions, waveform data, patient/study
  identifiers, and scanner measurement dependencies are omitted;
- the raw FIRE HDF5 layout is intentionally not the standard compound-dataset
  layout used by `ismrmrd.Dataset`;
- replay requires the raw-layout FIRE sender used for the validated offline
  workflow; the generic compound-dataset `ismrmrd.Dataset` client is not a
  reader for this archival layout.

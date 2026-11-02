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
Each image `position` is the physical center of that NIfTI slice. In-plane FOV
is matrix size times affine spacing, and header FOV Z is JSON `SliceThickness`
(or affine Z spacing). Permutations, reflections, and oblique orthogonal
affines are supported without resampling. Sheared affines are rejected because
MRD image direction cosines cannot represent shear exactly.

## Optional in-plane origin correction

`--inplane-origin-correction-voxels X Y` applies a physical translation to
every MRD image position. `X` is measured along the converted image's
`read_dir`; `Y` is measured along its `phase_dir`. The displacement is

```text
X * spacing_x * read_dir + Y * spacing_y * phase_dir
```

in LPS millimeters. The default is `0 0`, which is byte-identical to omitting
the option and preserves the original geometry. For example:

```bash
python3 shared/myhelpers_jauger/nifti_to_mrd.py \
  --nifti /path/to/scan.nii.gz \
  --json /path/to/scan.json \
  --output /path/to/scan_halfvoxel_corrected.hdf5 \
  --inplane-origin-correction-voxels 0.5 -0.5
```

Only header `position` and its duplicate `SlicePosLightMarker` attribute are
translated. Pixel arrays, dimensions, spacing, direction vectors, FOV, slice
ordering, timing, indices, and other metadata are unchanged.

For the 2026-10-15 REF acquisition, the measured converted-to-native index
mapping had residual `(-0.5000025, -0.4998563, -0.0000105)` relative to the
ideal X-preserving/Y-reflecting mapping. Because the converted Y direction is
reflected relative to native, cancelling the two in-plane residuals requires
`(+0.5, -0.5)` in the converted basis, not `(+0.5, +0.5)`. With its oblique
3 mm directions, that correction is an LPS displacement of approximately
`(+1.587035, +1.395713, +0.182501)` mm.

The exact validated conversion command was:

```bash
python3 shared/myhelpers_jauger/nifti_to_mrd.py \
  --nifti /home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/fromResearchSynapse_20261015/6045885-3-converted/6045885/20261015/joshauger_motioncorrection_20261008_2/11/ep2d_bold_fire_nomoco_tr1210_ref/ep2d_bold_fire_nomoco_tr1210_ref_ep2d_bold_fire_noMoCo_tr1210_ref_20261015134642_11.nii.gz \
  --json /home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/fromResearchSynapse_20261015/6045885-3-converted/6045885/20261015/joshauger_motioncorrection_20261008_2/11/ep2d_bold_fire_nomoco_tr1210_ref/ep2d_bold_fire_nomoco_tr1210_ref_ep2d_bold_fire_noMoCo_tr1210_ref_20261015134642_11.json \
  --output /home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/fromNiftiToMRDConverter_20261102/ref_convertedMRD_halfvoxel_corrected.hdf5 \
  --inplane-origin-correction-voxels 0.5 -0.5
```

Direct comparison with `func_handmotor_run02_REF_convertedMRD.hdf5` verified
all 6,400 pixel payloads, config/XML, directions, FOV, ordering, and nonposition
header/attribute fields were identical. Every header changed only `position`,
and every attribute changed only `SlicePosLightMarker`. Float32 header storage
produced a maximum displacement error of `6.32e-6` mm.

A complete offline FIRE replay produced
`data/savedData_20261102T190917_ep2d_bold_fire_noMoCo_tr1210_ref`. The comparison
report is in `data/image_comparison_REF_native_vs_halfvoxel_corrected`:

| Metric | Original converted REF | Half-voxel corrected REF |
| --- | ---: | ---: |
| Matched image pairs | 6,381 | 6,381 |
| Mean coordinate discrepancy | 2.121014886 mm | 0.000431798 mm |
| Maximum coordinate discrepancy | 2.121023229 mm | 0.000436805 mm |
| Physical support overlap A / B / IoU | 98.6153% / 98.6149% / 97.2680% | 99.9990% / 99.9985% / 99.9975% |
| Physical-alignment Pearson correlation | 0.955001349 | 0.998077712 |
| Physical-alignment MAE / RMSE | 34.688992 / 73.687336 | 6.437691 / 15.493410 |

The remaining roughly `0.00043` mm residual exceeds the comparison utility's
strict `0.0001` mm lossless-grid tolerance, so it is still classified as
requiring resampling even though it approaches the index-only reflection
result. Intensity equality is not expected because the source NIfTI comes from
downstream enhanced DICOMs rather than the scanner reconstruction.

This correction is dataset-dependent. Do not automatically apply the REF
offsets to unrelated acquisitions; determine and verify their signs and
magnitudes in the relevant physical coordinate frames first.

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

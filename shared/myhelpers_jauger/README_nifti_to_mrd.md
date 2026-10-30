# NIfTI-to-MRD converter

`nifti_to_mrd.py` converts one 4D NIfTI scan and its BIDS-style JSON
sidecar into the raw HDF5 layout written by `apps/python_fire_server`.
It is intended for offline analysis and future FIRE replay work; no FIRE
service needs to be running.

## Dependencies

The helper uses the repository's existing Python dependencies: `numpy`,
`nibabel`, `h5py`, and `ismrmrd`.

## Usage

```bash
python shared/myhelpers_jauger/nifti_to_mrd.py \
  --nifti /path/to/scan.nii.gz \
  --json /path/to/scan.json \
  --reference-mrd /path/to/FIRE-reference.hdf5 \
  --output /path/to/converted.hdf5
```

`--json` may be omitted when the sidecar has the NIfTI basename. Existing
outputs are protected unless `--overwrite` is supplied. The reference is
required and must be a FIRE raw-layout HDF5 file with one image series,
root `Config File` and `Metadata XML` datasets, and at least one complete
`image_N/{header,attribute,data}` group.

The October 6 validation command is:

```bash
python shared/myhelpers_jauger/nifti_to_mrd.py \
  --nifti /home/jauger/Radiology_Research/Scan_data/20261006_MRN6045885_handmotor_1_baseline/fromResearchSynapse_DICOMS_20261006/6045885-3-converted/6045885/20261006/joshauger_motioncorrection_handmotor_sesh1_baseline/18/func-bold_task-rest480_run-01_slimmon/func-bold_task-rest480_run-01_slimmon_func-bold_task-rest480_run-01_SLIMMON_20261006143118_18.nii.gz \
  --reference-mrd /home/jauger/Radiology_Research/Scan_data/20261006_MRN6045885_handmotor_1_baseline/fromPythonFireServer_20261006/savedData_20261006T191516_func-bold_task-rest480_run-01_SLIMMON/measurement-20261006T191403.hdf5 \
  --output /tmp/nifti_to_mrd_oct6_validation.hdf5
```

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

The copied XML structure is updated for reconstructed matrix/FOV, slice and
repetition limits, protocol, scanner fields, TR, TE, and flip angle. BIDS times
in seconds are converted to the MRD XML convention of milliseconds. Image
attributes are newly generated and do not copy scanner miniheaders. When TR
exists, `acquisition_time_stamp` is a relative value derived in standard
2.5 ms ticks from TR and SliceTiming; it is not a wall-clock acquisition time.
Scanner measurement/study identifiers are replaced with synthetic offline
identifiers, and image `measurement_uid` is zero.

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
- the raw FIRE HDF5 layout is intentionally not the standard compound-dataset
  layout used by `ismrmrd.Dataset`;
- existing replay utilities may still require a transport adapter that emits
  these stored image records as MRD wire messages.

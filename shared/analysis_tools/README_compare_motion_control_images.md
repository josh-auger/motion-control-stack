# Motion-control image comparison

`compare_motion_control_images.py` performs both a direct array comparison and
a SimpleITK physical-space comparison of images saved by two FIRE-server
acquisitions. It is intended for native-versus-converted MRD replay validation
and can also compare other `savedData_*` directory pairs. It does not need a
scanner, server, CUDA, registration program, or source MRD file.

## What is compared

The FIRE-server writes each received image payload unchanged to
`<protocol>_volume_####_slice_####.raw`. The first reference volume is exposed
through one detached 3D NRRD header; later received slices have individual
detached NRRD headers. The default discovery combines those representations so
every received voxel is compared once:

- `*_volume_####_YYYYMMDDTHHMMSS.nhdr`: reference 3D volume;
- `*_volume_####_slice_####.nhdr`: subsequent received slices;
- equivalent supported image headers with the same identifiers.

The volume and slice indices form the matching key. Protocol prefixes,
directory placement, timestamps, and lexicographic file position are not used
to pair images. Duplicate keys are ambiguous and cause a coverage failure.

Derived files containing markers such as `_upsampled`, `registered`,
`resampled`, or `_mask` are excluded from the primary comparison and listed in
the reports. Other supported image files without a reliable volume or
volume/slice identity are also reported as excluded. `.raw` payloads are not
independently discovered because their NRRD headers define shape, type, and
geometry.

NRRD, MetaImage, and NIfTI files are read with SimpleITK. The NumPy array is
obtained with `sitk.GetArrayFromImage()`. SimpleITK indexes an image as
`(x, y, z, ...)`, while the NumPy array axes are reversed as `(..., z, y, x)`;
the implementation keeps this conversion explicit whenever it applies an
index permutation. A direct raw-NRRD fallback is retained only for data such as
a 4D NRRD with a nonspatial axis that the installed ITK NRRD reader cannot
represent. Physical-space validation is unavailable for that fallback case.

## Equality rules

For compatible shapes, the utility reports exact element equality, dtype
equality, unequal voxel count and percentage, maximum absolute difference,
MAE, RMSE, mean signed difference (`A - B`), and Pearson correlation when it is
defined. A one-bit floating-point difference is a mismatch. Co-located NaNs
are treated as the same numerical value; NaN and signed infinity counts are
reported separately. Array-byte equality is reported independently, so cases
such as equal signed-zero values with different representations remain visible.

`exact_voxel_equal` requires equal shapes, equal dtypes, and exact numerical
values. Byte equality is diagnostic and is not substituted for voxel equality.
Statistics use only locations finite in both arrays; nonfinite mismatches still
count as unequal voxels.

Spatial metadata is evaluated separately using image dimensions, voxel
spacing, physical origin, direction cosine matrix, orientation code, image
center, center corners, voxel-support corners, physical bounding box, and
available reader metadata. Unrelated header comments and timestamps do not
affect voxel equality.

## Physical-coordinate conventions

SimpleITK normalizes supported image geometry to LPS physical coordinates in
millimeters. An image index maps to a physical point as:

```text
p = origin + direction @ diag(spacing) @ index
```

The utility validates this mapping through SimpleITK's
`TransformIndexToPhysicalPoint()` and
`TransformContinuousIndexToPhysicalPoint()`. Equal dimensions do not imply
equal physical coverage. Reports include image-center and corner coordinates,
physical support overlap, bounding-box differences, and discrepancies between
coordinates implied by the best discrete axis relation.

The physical comparison has three deliberately separate levels:

1. **Direct array comparison** compares the original NumPy arrays without any
   axis or intensity change.
2. **Lossless physical reorientation** detects a signed axis permutation,
   identifies the nearest integer source-to-target index offset, and verifies
   corresponding physical coordinates. Exact equivalence and equivalence
   within `--spatial-tolerance-mm` are reported separately. When the residual
   is accepted by the configured tolerance, the discrete transpose/reflection
   is preferred over interpolation. No pixel values are changed.
3. **Interpolated physical comparison** is used when the grids share
   SimpleITK's LPS coordinate frame but do not have corresponding voxel
   centers. The converted image is resampled onto the native grid with an
   explicit identity physical transform and linear interpolation. Statistics
   use only target voxel centers within the source support. Physical support
   overlap and sampled-center overlap are reported separately.

An orientation label such as LPS versus LAS is not by itself proof that a
NumPy Y reflection is correct. The full direction, spacing, origin, index
offset, sizes, and physical point mapping must agree. The report includes an
index-only signed-axis diagnostic for investigating such cases, but labels it
as nonphysical whenever the required integer-grid mapping fails.

After physical alignment, a diagnostic least-squares model
`native = scale × converted + offset` is reported with residual MAE and RMSE.
It never modifies source intensities or the primary comparison. Intensity
distributions also report the fraction of integer-valued samples to expose
quantization-like behavior.

Legacy direct-array categories are retained for downstream compatibility:

- `PASS_EXACT`: complete matching, exact voxels, and identical geometry;
- `PASS_VOXELS_ONLY`: complete matching and exact voxels, with geometry differences;
- `FAIL_VOXELS`: at least one matched array differs in shape, dtype, or voxel values;
- `FAIL_COVERAGE`: missing, extra, ambiguous, or zero matched identities;
- `FAIL_GEOMETRY`: at least one matched pair has different spatial geometry.

Spatial and pixel validation use independent statuses.

Spatial statuses:

- `SPATIAL_PASS_EXACT`: dimensions, spacing, physical coverage, and the
  discrete voxel-center mapping agree to strict numerical precision;
- `SPATIAL_PASS_TOLERANCE`: the same structural checks pass and corresponding
  physical coordinates are within `--spatial-tolerance-mm`, but are not exact;
- `SPATIAL_FAIL_GRID`: dimensions, spacing, direction relationship, or
  voxel-center mapping is incompatible or exceeds the tolerance;
- `SPATIAL_FAIL_COVERAGE`: the physical fields of view do not agree within the
  tolerance or have insufficient overlap;
- `SPATIAL_INDETERMINATE`: SimpleITK cannot represent the needed physical
  geometry.

The spatial tolerance applies to coordinate residuals, not to dimensions or
materially different voxel spacing. Actual residuals are always reported.
Intensity correlation never changes a spatial status.

Pixel statuses:

- `PIXEL_PASS_EXACT`: values and dtypes are exact either directly or after a
  verified lossless reindexing;
- `PIXEL_SIMILAR`: values are nonidentical but the primary aligned comparison
  meets `--pixel-correlation-threshold`;
- `PIXEL_DIFFERENT`: comparable values do not meet that criterion;
- `PIXEL_INCOMPARABLE`: no meaningful direct or aligned statistics exist.

An interpolated comparison is never called exact. The fitted scale and offset
remain diagnostic and do not affect primary metrics or pixel values.
Per-pair pixel statuses are retained in the CSV; the overall `PIXEL_SIMILAR`
decision uses the aggregate aligned correlation across all valid overlap
voxels, so low-information individual slices do not dominate the dataset-level
assessment.

High post-alignment correlation is evidence of similar image content, not
proof that the original physical headers were correct.

When spatial correspondence passes but pixels are merely similar or different,
the metadata support comparing retrospective motion estimates in a consistent
physical coordinate frame. That conclusion does not guarantee identical
registration outputs and does not independently prove absolute scanner-space
accuracy.

## Command line

Only the input directories are required:

```bash
python3 shared/analysis_tools/compare_motion_control_images.py \
  --dir-a /path/to/savedData_A \
  --dir-b /path/to/savedData_B
```

Without `--output-dir`, reports are written under the current working directory
as `comparison_results_<dir-a-name>_vs_<dir-b-name>`. An output directory equal
to or inside an input directory is rejected, preventing reports from changing
future discovery results.

Options:

- `--output-dir PATH`: report destination;
- `--label-a TEXT`, `--label-b TEXT`: human-readable dataset labels;
- `--image-pattern GLOB`: filter supported candidates by relative path or basename;
- `--spatial-tolerance-mm FLOAT`: coordinate acceptance tolerance, default
  `0.001` mm;
- `--pixel-correlation-threshold FLOAT`: threshold for `PIXEL_SIMILAR`, default
  `0.99`;
- `--validation-target spatial|pixels|both`: controls only the process exit
  code; all diagnostics are still calculated. Coverage failures always fail;
- `--verbose`: emit periodic comparison progress.

With `--validation-target spatial`, exit code 1 means coverage or spatial
failure. With `pixels`, it means coverage or `PIXEL_DIFFERENT`/
`PIXEL_INCOMPARABLE`. With `both` (the default), either required assessment can
fail the process. Input and usage errors return 2.

REF example from the repository root:

```bash
python3 shared/analysis_tools/compare_motion_control_images.py \
  --dir-a data/savedData_20261102T154607_ep2d_bold_fire_noMoCo_tr1210_ref \
  --dir-b data/savedData_20261102T190917_ep2d_bold_fire_noMoCo_tr1210_ref \
  --label-a "Native REF MRD" \
  --label-b "Corrected converted REF MRD" \
  --spatial-tolerance-mm 0.001 \
  --pixel-correlation-threshold 0.99 \
  --validation-target both \
  --output-dir data/image_comparison_REF_spatial_and_pixel
```

## Reports

The output directory contains:

- `comparison_summary.json`: inputs, discovery and coverage details, excluded
  files, aggregate metrics, geometry status, categories, and validation flags;
- `image_comparison.csv`: one row per unambiguous matched image identity;
- `comparison_report.md`: concise human-readable findings and the first
  significant discrepancies;
- `figures/`: at most two representative/worst-case panels containing the
  native image, converted storage ordering, physically aligned converted
  image, difference map, and diagnostic intensity fit.

## Limitations

- Default matching requires FIRE-style `volume` and optional `slice` indices.
  Files without a reliable identifier are reported, not positionally paired.
- Image I/O and spatial interpretation depend on SimpleITK/ITK support. The
  direct NRRD fallback has no physical validation.
- Signed permutations/reflections are detected with a small numerical matrix
  tolerance. Spatial acceptance uses the separately reported configurable
  physical-coordinate tolerance.
- General rotated-grid physical overlap uses an axis-aligned bounding-box
  approximation; signed-axis relations use exact support bounds in target
  continuous-index space.
- Linear resampling is diagnostic and is never classified as exact or
  lossless. Discrete label images would require nearest-neighbor interpolation;
  this utility assumes continuous MRI intensities.
- Differences observed after replay cannot by themselves identify which stage
  introduced them. In particular, upstream DICOM-to-NIfTI conversion may be a
  contributor and the MRD converter should not be assumed to be the cause.

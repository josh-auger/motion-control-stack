# Three-way REF image comparison

`compare_three_way_mrd_images.py` compares three representations of the same
FIRE reference acquisition:

1. scanner-native FIRE-server NRRDs;
2. the direct 4D DICOM-to-NRRD export from `dcm2niix`;
3. FIRE-server NRRDs produced by replaying the corrected NIfTI-derived MRD.

The analysis is read-only with respect to all three input datasets. It writes a
JSON summary, a per-record CSV, a Markdown report, and four diagnostic figures.

## Usage

From the repository root:

```bash
python3 shared/analysis_tools/compare_three_way_mrd_images.py \
  --native-fire /path/to/native/savedData_REF \
  --dicom-nrrd /path/to/direct_dicom_export.nrrd \
  --converted-fire /path/to/convertedMRD/savedData_REF \
  --output-dir data/image_comparison_REF_three_way \
  --spatial-tolerance-mm 0.001
```

For the supplied REF experiment:

```bash
python3 shared/analysis_tools/compare_three_way_mrd_images.py \
  --native-fire /home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/fromPythonFireServer_20261015/savedData_20261015T183846_ep2d_bold_fire_noMoCo_tr1210_ref \
  --dicom-nrrd /home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/fromResearchSynapse_20261015/6045885-4-nrrd/6045885/20261015/joshauger_motioncorrection_20261008_2/11/ep2d_bold_fire_nomoco_tr1210_ref/ep2d_bold_fire_nomoco_tr1210_ref_ep2d_bold_fire_noMoCo_tr1210_ref_20261015134642_11.nrrd \
  --converted-fire /home/jauger/Radiology_Research/Scan_data/20261015_MRN6045885_handmotor_3_pMoCo/fromRetroSVR_20261102/savedData_20261102T193244_ep2d_bold_fire_noMoCo_tr1210_ref_convertedMRD \
  --output-dir data/image_comparison_REF_three_way \
  --spatial-tolerance-mm 0.001
```

The tool refuses to write into a nonempty output directory unless
`--overwrite` is supplied. It never edits input files.

## Correspondence rules

- FIRE records are matched by parsed volume and slice identifiers, never by
  directory position.
- The temporal offset is accepted only when complete, contiguous FIRE coverage
  matches the DICOM time dimension. Partial or duplicate coverage is reported
  as indeterminate.
- FIRE slice ordinals are mapped to anatomical Z indices from physical slice
  origins in the volume-0 reference geometry. Reversed or interleaved orders
  are therefore handled explicitly.
- The DICOM NRRD's declared RAS coordinates are exposed by SimpleITK in LPS;
  this conversion negates physical X and Y.
- Pixel comparisons use lossless reorientation only when physical geometry
  passes the configured tolerance. If it does not, signed-axis index metrics
  are labeled as index-only diagnostics rather than physical equivalence.

## Outputs

```text
data/image_comparison_REF_three_way/
├── comparison_report.md
├── comparison_summary.json
├── pairwise_comparison.csv
└── figures/
    ├── figure1_three_way_anatomy.png
    ├── figure2_pairwise_differences.png
    ├── figure3_intensity_scatterplots.png
    └── figure4_spatial_geometry.png
```

The JSON report preserves input paths, mappings, geometry, tolerances,
alignment methods, aggregate statistics, intensity distributions, and the
diagnostic affine intensity fits. The CSV contains one row for every pairwise
FIRE record comparison.

## Interpretation cautions

The tool keeps spatial validation independent from pixel similarity. High
correlation does not establish physical or temporal correspondence. It also
does not normalize intensity values for the primary comparison. Any affine
intensity fit is diagnostic only.

Both the direct NRRD and the NIfTI used by the conversion pathway originate in
`dcm2niix`, so their agreement establishes preservation downstream of that
shared decoder; it does not independently isolate `dcm2niix`'s upstream DICOM
scaling behavior.

## Native-versus-DICOM intensity diagnostics

After generating the three-way outputs, run the focused intensity investigation
with the same inputs:

```bash
python3 shared/analysis_tools/diagnose_native_dicom_intensities.py \
  --native-fire /path/to/native/savedData_REF \
  --dicom-nrrd /path/to/direct_dicom_export.nrrd \
  --converted-fire /path/to/convertedMRD/savedData_REF \
  --output-dir data/image_comparison_REF_three_way
```

The utility automatically discovers a unique `measurement-*.hdf5` capture in
each FIRE directory and derives the corresponding `-2-sorted` enhanced-DICOM
series from the `-4-nrrd` path. Explicit `--native-mrd`, `--converted-mrd`, and
`--dicom-series-dir` arguments are available when automatic discovery is not
appropriate.

It appends a separate investigation section to `comparison_report.md` and
adds:

```text
intensity_diagnostics.json
intensity_per_volume.csv
intensity_representative_slices.csv
figures/figure5_native_dicom_intensity_relationship.png
figures/figure6_intensity_residual_distributions.png
figures/figure7_scaling_stability.png
figures/figure8_representative_intensity_profiles.png
```

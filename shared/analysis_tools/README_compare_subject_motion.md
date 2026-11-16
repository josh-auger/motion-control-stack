# Retrospective six-subject NOCO motion comparison

The primary motion-change metric is **consecutive slice-group displacement (mm)**
at the native temporal resolution of every registration estimate.

Run from the repository root:

```bash
python3 shared/analysis_tools/compare_subject_motion.py
python3 -m pytest shared/analysis_tools/tests/test_compare_subject_motion.py shared/analysis_tools/test_analysis_tools.py -q
```

Run tests after the analyzer has exited, so artifact-consistency checks read a
completed export. Add `--junitxml=data/GBS_six_subject_comparison/validation_tests.xml`
to save the automated validation results beside the report. The artifact audit
is skipped when the generated report is absent.

The default input root is `data/`; the output is
`data/GBS_six_subject_comparison/`. Override these using `--data-root` and
`--output-dir`. Discovery requires six unique `savedData_*_GBS<digits>`
directories; repeated subjects are ambiguous and rejected. For synthetic or other
cohorts, `--expected-subjects N` changes that discovery check. The current acquisition
geometry requires 14 groups. Original replay directories are read-only inputs; an
output inside, equal to, or containing any input directory is refused.

Dependencies are already in the repository requirements: NumPy, SimpleITK,
Matplotlib, h5py, and ismrmrd, plus Python's standard library. No production
services are imported or changed. Existing `analysis_common`, registration-status,
and registration-runtime helpers are reused. Pure monitor arithmetic is reproduced
in the analysis module and checked against the existing function bodies in tests.

## Observations and derived measurements

* Registration outcomes/finalization: `registration_status.json`, schema 1.
* Rigid transforms: exact status-linked `.tfm` names, searched at the results root
  and in `processed_*/transforms/`. Conflicting duplicate files are refused.
* Acquisition geometry/expected counts: encoding limits, TR, and `SliceTiming` in
  the metadata JSON. A unique reference NRRD header confirms physical LPS space.
* Slice groups: pointer-file members are arrival ordinals. They are not spatial
  slice indices into the metadata's `SliceTiming` object.
* Acquisition time: MRD image-header `acquisition_time_stamp` in 400 Hz ticks,
  relative to the first received slice. The first image of each group supplies its
  time; coverage of all four images is checked. Metadata-derived times are used
  only when the MRD file is absent, with a warning. Missing/invalid stored times
  remain missing. For these replays the converter derived MRD times from the source
  NIfTI acquisition metadata; these are not measurements of scanner wall time.
* Runtime: queue-item monotonic registration end minus start, in milliseconds.
  All actual calls are included, retaining startup overhead inside the call.
  Identity reference setup is excluded. Runtime does not establish acquisition time.
* Motion CSV: an independent saved ledger for parameter/displacement cross-checks;
  its recorded group displacement is exported separately from recomputed metrics.

All non-reference volumes retain 14 group estimates. Volume 0 has one recorded
identity reference, not 14 estimated transforms. Accordingly the six current runs
each have 480 volumes, 6,706 expected calls, and 6,707 status/motion observations.
Expected opportunities are derived from metadata rather than the number of
successful records, so absent status records cannot disappear from the denominator.

Translations are saved ITK center-relative translation parameters in mm; rotation
trajectories are saved Euler parameters in degrees, preserving `ComputeZYX`.
Reference-relative translation magnitude is their Euclidean norm. The primary
rotation magnitude is the Euclidean norm of the Euler angles, matching existing
plotting helpers; a separate matrix axis-angle magnitude is also exported. These
reference-relative quantities are distinct from displacement between observations.

## Consecutive slice-group displacement formulation

Use the saved ITK mapping `T(p) = R(p-c) + c + t`, in physical LPS coordinates.
The relative mapping is **current composed with inverse(previous)**:

```text
Delta = T_current o inverse(T_previous)
R_delta = R_current R_previous^T
c_delta = c_previous
t_delta = R_current(-R_previous^T t_previous + c_previous - c_current)
          + t_current + c_current - c_previous
```

As in `apps/motion_monitor/monitor_directory.py`, decompose Delta with the default
SimpleITK Euler3D convention (`ComputeZYX=False`, matrix order Rz Rx Ry).
For its Euler angles `(x,y,z)` in radians, retain the monitor's exact formulation:

```text
q = [-1 + cos(x)cos(y) + cos(x)cos(z) + cos(y)cos(z)
     + sin(x)sin(y)sin(z)] / 2
theta = acos(clip(q, -1, 1))
D = ||t_delta||_2 + r sqrt((1-cos(theta))^2 + sin(theta)^2)
r = 50 mm
```

Clipping protects against floating-point roundoff. This deliberately retains the
repository's legacy Euler-based rotational proxy, including its mixed-angle term;
it is not conventional Power FD or a matrix-geodesic angle.

**Consecutive slice-group displacement** applies D to adjacent expected group
slots in acquisition sequence, retaining all 14 estimates per non-reference
volume and every transition across volume boundaries. The recorded identity for
volume 0 supplies the initial baseline, timed at the last reference slice group;
its own displacement is undefined and excluded from displacement statistics.
For each complete subject, this gives **6,706 valid consecutive-group transitions**.
There is no subsampling or averaging. This metric describes change between
registration estimates and is distinct from absolute reference-relative motion.

Missing, skipped, failed, and unreadable transforms stay missing. The first
displacement is undefined. A missing transform invalidates both adjacent pairs;
gaps are never bridged. Volume completeness is reported independently.
Reference identity is also excluded from reference-relative magnitude distributions.

Motion events use **strictly D > 0.3 mm**, a descriptive threshold matching the
repository launcher, configurable with `--motion-event-threshold-mm`. Counts and
percentages use valid transitions; `event_runs` counts consecutive exceedances,
with missing observations splitting runs. Group events additionally report the
number of affected volumes. These are not image-quality classifications.
The HTML displacement summary table displays **Affected volumes (>0.3 mm)**,
using the existing `group_event_affected_volumes` value. It counts unique acquired
volumes containing at least one consecutive slice-group displacement exceeding
0.30 mm; each volume counts once. Cross-volume transitions are included and
attributed to the volume containing the current slice group. The `event_runs`
statistic remains unchanged in the CSV/JSON exports and embedded summaries.
Statistics use finite observations only, population SD, and linearly interpolated
percentiles. `missing_count` includes the initial undefined displacement.
Cross-subject distributions of subject means weight subjects equally.

## Exceedance by slice-group position

The **High-motion exceedance by slice-group position** section follows the
consecutive-displacement summary table. It groups the existing displacement
measurements by the current estimate's explicit volume and slice-group metadata;
it does not recompute motion or change registration ordering. Consecutive metadata
pairs are checked before classification. Group 0 is the cross-volume transition
from the preceding volume's final estimate to the current volume's first estimate;
the initial such transition uses the recorded volume-0 identity. Groups 1–13 are
intra-volume transitions. No within-reference observations are invented.

For each group, exceedance percentage is `100 * count(D > 0.30 mm) / count(valid D)`.
Missing values do not enter the denominator; a group with no valid measurements
has a missing rate, mean, and P95. Each complete subject contributes 479 valid
measurements per group, totaling 6,706. Mean and P95 are calculated over all valid
measurements in the group, including those below threshold. Exceedances indicate
elevated estimated motion between registration measurements, not necessarily
image corruption. This is a descriptive positional comparison.

## Artifacts

* `report.html`: self-contained HTML with six embedded PNGs, tables, methods,
  warnings, and provenance. No external JavaScript, fonts, or image dependencies.
* `figure_1_motion_trajectories.png`: six aligned subject rows, translation and
  rotation columns, shared acquisition-time and per-unit axes; every group shown.
* `figure_2_group_displacement.png`: six vertically aligned subject panels showing
  every consecutive slice-group displacement against subject-specific acquisition
  time, with a common vertical scale and the descriptive threshold.
* `figure_3_motion_distributions.png`: consecutive slice-group displacement
  boxplots with mean diamonds and a table of mean, median, P95, maximum, and
  percentage of valid transitions exceeding the threshold.
* `figure_4_registration_performance.png`: registered/failed/skipped/missing counts.
* `figure_5_registration_runtime.png`: full-call distributions, startup retained.
* `slice_group_exceedance_by_position.png`: one 2×3 figure in GBS024, GBS029,
  GBS061, GBS067, GBS069, GBS071 order, with 14 bars per subject, group 0 highlighted,
  and a common percentage scale. The original five figures are unchanged.
* `slice_group_exceedance_by_position.csv`: supporting `subject`,
  `slice_group_index`, `transition_type` (`cross_volume` or `intra_volume`),
  `valid_transition_count`, `above_threshold_count`, `exceedance_pct`,
  `mean_displacement_mm`, and `p95_displacement_mm` columns.
* `summary.json` (schema 2): metadata, completeness, group statistics, runtimes,
  cross-subject comparisons, provenance, assumptions, and warnings.
* `subject_summary.csv`: metadata and processing completeness.
* `motion_statistics.csv`: group motion and runtime statistics, including
  mean, median, P95, maximum and event frequencies.
* `group_motion.csv` / `group_motion.json`: every expected group slot, outcome,
  time/source, transform source/type/center/convention, motion, displacement,
  runtime, and recorded CSV displacement.
* `validation.json` and `validation_tests.xml`: retained validation records,
  refreshed after running the focused tests and checking source preservation.

The analyzer removes the three known superseded sampling artifacts from its output
directory after a successful report write. It leaves validation records and
unrelated files intact. No original replay directory is modified.

Missing numeric values are CSV blanks and JSON nulls; plotting preserves gaps.
Figures 2 and 3 use a common full-range symlog scale when a large outlier would
compress the typical observations (maximum > five times max(threshold, 0.3 mm)).
The scale is linear through max(threshold, 0.3 mm), logarithmic above that value,
and explicitly labeled. All outliers remain plotted.
An incompatible subject remains in summary tables with explicit warnings rather
than being silently omitted. Discovery ambiguity/count errors stop the command.

This is a descriptive comparison of independently acquired **NOCO** scans.
Recorded registration success does not establish registration accuracy.
The report makes no claim about prospective correction effectiveness or
motion-corrected image quality.

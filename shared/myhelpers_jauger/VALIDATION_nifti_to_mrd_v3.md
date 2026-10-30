# Reference-free NIfTI-to-MRD validation (V3)

## Reference-dependency audit

The V2 reference supplied five kinds of information. None require an external
scan once their role is made explicit:

| Former dependency | Classification | V3 source |
| --- | --- | --- |
| Raw HDF5 hierarchy and byte datasets | Required structural metadata | Fixed FIRE layout written directly with `h5py` |
| XML matrix, FOV, limits, protocol, and sequence fields | Derivable | NIfTI shape/affine and JSON sidecar |
| Header version/flags, magnitude type, series 0, float32 preference | Defaultable | Documented FIRE defaults |
| 1,024-byte config filename | Defaultable | `fire_emitter_dcm` plus null padding |
| Scanner UIDs, ICE miniheaders, coils, dependencies, absolute times | Unavailable/private | Omitted; never fabricated |

The minimum offline replay contract is: a parseable MRD XML header containing
one encoding, protocol name, slice/repetition limits, and a user-parameter
multiband factor; one 1,024-byte config payload; ordered
`image_N/{header,attribute,data}` byte datasets; and the FIRE direction Meta
attributes described in the V2 report. Queue and motion-monitor processing does
not read the MRD XML or HDF5 directly; it consumes outputs published by FIRE.

## Reference-free XML

V3 builds `ismrmrd.xsd.ismrmrdHeader` objects and serializes them with
`ismrmrd.xsd.ToXML`. It populates:

- encoded and reconstructed `90 x 90 x 1` matrix and image FOV;
- Cartesian trajectory and encoding limits for 56 slices and 50 repetitions;
- protocol, vendor, model, field strength, and optional patient position;
- TR, TE, flip angle, sequence name, echo spacing, and echo-train length when
  present;
- `MultiBandFactor=4`, taken from JSON (or uniformly derivable SliceTiming);
- `SpacingBetweenSlices=2.4` mm and converter/source provenance.

The measurement identifier is deterministic and explicitly synthetic. Patient,
study, series, frame-of-reference, and scanner measurement UIDs are omitted.
Missing optional scanner-system values use `Unknown` and field strength `0.0`
only where FIRE's metadata-save path requires the element to be present.

## Generated artifact

```text
/home/jauger/Radiology_Research/Scan_data/20261006_MRN6045885_handmotor_1_baseline/dataConversionTesting_Nifti_to_MRD_rest480_run02/nifti_to_mrd_oct6_reference_free_v3.hdf5
```

Prior FIRE, V1, and V2 files were not overwritten.

## V2 versus V3

| Property | V2 reference-assisted | V3 reference-free |
| --- | ---: | ---: |
| File size | 100,487,584 bytes | 100,271,222 bytes |
| Root objects | 2,802 | 2,802 |
| Image records | 2,800 | 2,800 |
| XML bytes | 219,621 | 3,259 |
| Config payload | `fire_emitter_dcm` | byte-identical |
| Image headers | validated V2 | byte-identical for all 2,800 |
| Attributes | validated V2 | byte-identical for all 2,800 |
| Pixel payloads | validated V2 | byte-identical for all 2,800 |

The XML size reduction is expected: V3 omits scanner-private protocol blobs,
coil labels, UIDs, dependency records, and unrelated acquisition fields. Both
V2 and V3 contain the same interpretation-critical values: protocol, Siemens
Prisma description, 3 T field strength, matrix/FOV, Cartesian trajectory,
slice/repetition limits, TR 1000 ms, TE 29 ms, flip 52 degrees, and multiband
factor 4. V3 additionally records `SpacingBetweenSlices=2.4` mm directly from
JSON.

## Validation results

- Exact NIfTI round trip for all 22,680,000 voxels; maximum error `0`
- Correct 50-volume x 56-slice organization and all 2,800 index pairs
- V2/V3 image header, attribute, and pixel bytes identical
- V3 XML parsed by `ismrmrd.xsd.CreateFromDocument`
- FIRE protocol, multiband factor, and slice-count metadata accessors returned
  the expected values from reference-free XML
- `SpacingBetweenSlices` parsed as 2.4 mm; FIRE's current accessor has a
  case-comparison defect and intentionally falls back to deriving spacing from
  the first complete volume
- All FIRE-required Meta attributes validated during conversion
- FIRE `get_image_direction_matrix()` succeeded for images 0, 1, 55, 56,
  1399, and 2799 and matched their header direction matrices within `1e-6`
- Raw hierarchy, config, group names, and image payloads match the successfully
  replayed V2 file, preserving the offline sender's file-reading contract
- Automated synthetic tests: 19 passed

No full socket replay was performed for V3, so end-to-end replay success is not
claimed. The remaining integration check is a complete V3 replay through all
50 volumes. Scanner-private metadata remains unavailable, but it is not used by
the validated FIRE reference-volume path.

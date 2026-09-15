from pathlib import Path

import numpy as np
import SimpleITK as sitk

from apps.motion_monitor.running_tsnr import (
    RunningTSNR,
    atomic_save_tsnr_image,
    create_tsnr_mosaic,
)


# -------------------------------------------------------------------------
# USER CONFIGURATION
# -------------------------------------------------------------------------

DATA_DIR = Path(
    "/home/jauger/Radiology_Research/Scan_data/20251023_MRN6045885_moco_scan/fromPythonFireServer_20251023/scan002_savedData_20251023T170713_fetal1_SMS1_TR2.5"
)

OUTPUT_DIR = Path(f"{DATA_DIR}/tsnr_test_results")

DISPLAY_MAX = 100.0

CHECKPOINTS = {
    2,
    5,
    10,
    20,
    50,
}


# -------------------------------------------------------------------------
# INPUT DISCOVERY
# -------------------------------------------------------------------------

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

nhdr_files = sorted(
    path
    for path in DATA_DIR.glob("*.nhdr")
    if "_upsampled" not in path.name
)

print(f"Found {len(nhdr_files)} candidate NRRD headers.")

if len(nhdr_files) < 2:
    raise RuntimeError(
        "Need at least two 3D volumes to calculate TSNR."
    )

for path in nhdr_files[:10]:
    print(" ", path.name)


# -------------------------------------------------------------------------
# RUNNING TSNR
# -------------------------------------------------------------------------

running = RunningTSNR()

# Only retained for offline validation against NumPy.
# The real-time implementation does NOT retain prior volumes.
all_volumes = []

reference_shape = None
reference_size = None
reference_spacing = None
reference_direction = None


for index, path in enumerate(nhdr_files):

    print(
        f"\nReading {index + 1}/{len(nhdr_files)}: "
        f"{path.name}"
    )

    image = sitk.ReadImage(str(path))
    volume = sitk.GetArrayFromImage(image)

    print("  SITK size:     ", image.GetSize())
    print("  NumPy shape:   ", volume.shape)
    print("  spacing:       ", image.GetSpacing())

    if volume.ndim != 3:
        raise RuntimeError(
            f"{path.name} is not a 3D volume. "
            f"Shape = {volume.shape}"
        )

    if reference_shape is None:

        reference_shape = volume.shape
        reference_size = image.GetSize()
        reference_spacing = image.GetSpacing()
        reference_direction = image.GetDirection()

        print("\nReference geometry:")
        print("  shape:     ", reference_shape)
        print("  size:      ", reference_size)
        print("  spacing:   ", reference_spacing)
        print("  direction: ", reference_direction)

    else:

        if volume.shape != reference_shape:
            raise RuntimeError(
                f"Shape mismatch for {path.name}: "
                f"{volume.shape} != {reference_shape}"
            )

        if image.GetSize() != reference_size:
            raise RuntimeError(
                f"Image size mismatch for {path.name}"
            )

        if not np.allclose(
            image.GetSpacing(),
            reference_spacing,
        ):
            raise RuntimeError(
                f"Spacing mismatch for {path.name}"
            )

        if not np.allclose(
            image.GetDirection(),
            reference_direction,
        ):
            raise RuntimeError(
                f"Direction mismatch for {path.name}"
            )

    volume = volume.astype(np.float64)

    running.update(volume)

    # Keep data only for validation against direct NumPy TSNR.
    all_volumes.append(volume)

    print(f"  TSNR sample count: {running.count}")

    if running.count in CHECKPOINTS:

        tsnr = running.get_tsnr()

        mosaic = create_tsnr_mosaic(
            tsnr,
            output_shape=(
                volume.shape[1],
                volume.shape[2],
            ),
            display_max=DISPLAY_MAX,
        )

        output_path = (
            OUTPUT_DIR
            / f"tsnr_n{running.count:04d}.jpg"
        )

        atomic_save_tsnr_image(
            mosaic,
            output_path,
        )

        print(
            f"  Saved checkpoint: {output_path}"
        )


# -------------------------------------------------------------------------
# SAVE FINAL TSNR MOSAIC
# -------------------------------------------------------------------------

final_tsnr = running.get_tsnr()

final_mosaic = create_tsnr_mosaic(
    final_tsnr,
    output_shape=(
        reference_shape[1],
        reference_shape[2],
    ),
    display_max=DISPLAY_MAX,
)

final_path = OUTPUT_DIR / "tsnr_final.jpg"

atomic_save_tsnr_image(
    final_mosaic,
    final_path,
)

print(
    f"\nSaved final TSNR mosaic: {final_path}"
)


# -------------------------------------------------------------------------
# OFFLINE NUMPY VALIDATION
# -------------------------------------------------------------------------

print("\nBuilding conventional NumPy TSNR...")

stack = np.stack(
    all_volumes,
    axis=0,
)

expected_mean = np.mean(
    stack,
    axis=0,
)

expected_std = np.std(
    stack,
    axis=0,
    ddof=1,
)

expected_tsnr = np.divide(
    expected_mean,
    expected_std,
    out=np.zeros_like(expected_mean),
    where=expected_std > 0,
)

expected_tsnr = np.nan_to_num(
    expected_tsnr,
    nan=0.0,
    posinf=0.0,
    neginf=0.0,
)

np.testing.assert_allclose(
    final_tsnr,
    expected_tsnr,
    rtol=1e-10,
    atol=1e-10,
)

print("\nPASS: RunningTSNR matches direct NumPy TSNR.")
print(f"Volumes processed: {running.count}")
print(f"Final TSNR shape:  {final_tsnr.shape}")
print(f"Final mosaic shape: {final_mosaic.shape}")
"""Incremental temporal SNR calculation and scanner-sized display helpers."""

from __future__ import annotations

import os
import tempfile
from os import PathLike

import cv2
import numpy as np


TSNR_AXIAL_POSITIONS = (0.20, 0.40, 0.60, 0.80)


class RunningTSNR:
    """Maintain voxelwise running statistics without retaining prior volumes.

    Input arrays use the SimpleITK NumPy convention ``(z, y, x)``. Statistics
    are accumulated as float64 using Welford's online algorithm.
    """

    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        """Clear all accumulated statistics."""
        self.count = 0
        self.mean: np.ndarray | None = None
        self.M2: np.ndarray | None = None

    def update(self, volume: np.ndarray) -> None:
        """Add one 3D volume to the running statistics.

        Raises:
            ValueError: If ``volume`` is not 3D, is complex-valued, or does
                not match the shape established by the first volume.
            TypeError: If ``volume`` does not contain numeric data.
        """
        array = np.asarray(volume)
        if array.ndim != 3:
            raise ValueError(f"Expected a 3D volume, got shape {array.shape}.")
        if not np.issubdtype(array.dtype, np.number):
            raise TypeError(f"Expected numeric volume data, got dtype {array.dtype}.")
        if np.iscomplexobj(array):
            raise ValueError("Complex-valued volumes are not supported; supply magnitude data.")

        if self.mean is not None and array.shape != self.mean.shape:
            raise ValueError(
                f"Volume shape {array.shape} does not match initial shape "
                f"{self.mean.shape}."
            )

        values = array.astype(np.float64, copy=False)
        if self.count == 0:
            self.mean = values.copy()
            self.M2 = np.zeros(values.shape, dtype=np.float64)
            self.count = 1
            return

        # Shape validation above and count > 0 guarantee these are arrays.
        assert self.mean is not None
        assert self.M2 is not None
        self.count += 1
        with np.errstate(invalid="ignore", over="ignore"):
            delta = values - self.mean
            self.mean += delta / self.count
            delta2 = values - self.mean
            self.M2 += delta * delta2

    def get_variance(self) -> np.ndarray:
        """Return voxelwise sample variance, or zeros after only one volume."""
        if self.mean is None or self.M2 is None:
            raise RuntimeError("No volumes have been accumulated.")
        if self.count < 2:
            return np.zeros_like(self.mean)

        with np.errstate(invalid="ignore", over="ignore"):
            variance = self.M2 / (self.count - 1)
        return np.nan_to_num(variance, nan=0.0, posinf=0.0, neginf=0.0)

    def get_std(self) -> np.ndarray:
        """Return voxelwise sample standard deviation."""
        with np.errstate(invalid="ignore"):
            return np.sqrt(self.get_variance())

    def get_tsnr(self) -> np.ndarray:
        """Return voxelwise temporal SNR, replacing undefined values with zero."""
        if self.mean is None:
            raise RuntimeError("No volumes have been accumulated.")
        if self.count < 2:
            return np.zeros_like(self.mean)

        std = self.get_std()
        tsnr = np.zeros_like(self.mean)
        valid = (std > 0) & np.isfinite(std) & np.isfinite(self.mean)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            np.divide(self.mean, std, out=tsnr, where=valid)
        return np.nan_to_num(tsnr, nan=0.0, posinf=0.0, neginf=0.0)


def create_tsnr_mosaic(
    tsnr_volume: np.ndarray,
    output_shape: tuple[int, int] | None = None,
    display_max: float = 100.0,
) -> np.ndarray:
    """Create an exact-size uint8 2x2 mosaic of axial TSNR slices.

    ``tsnr_volume`` is expected in ``(z, y, x)`` order. The y-axis is flipped
    for display, matching the existing MRI display helper used by this
    repository. ``output_shape`` is ``(height, width)`` and defaults to the
    dimensions of one input slice.
    """
    volume = np.asarray(tsnr_volume)
    display_slices = get_tsnr_display_slices(volume, display_max=display_max)

    if output_shape is None:
        height, width = int(volume.shape[1]), int(volume.shape[2])
    else:
        if len(output_shape) != 2:
            raise ValueError("output_shape must contain (height, width).")
        height, width = (int(output_shape[0]), int(output_shape[1]))
        if height < 2 or width < 2:
            raise ValueError("output_shape height and width must both be at least 2.")

    scaled_slices = np.rint(display_slices * (255.0 / display_max)).astype(np.uint8)
    row_heights = (height // 2, height - height // 2)
    column_widths = (width // 2, width - width // 2)
    mosaic = np.zeros((height, width), dtype=np.uint8)

    for tile_number, display_slice in enumerate(scaled_slices):
        row, column = divmod(tile_number, 2)
        tile_height = row_heights[row]
        tile_width = column_widths[column]
        interpolation = (
            cv2.INTER_AREA
            if tile_height < volume.shape[1] or tile_width < volume.shape[2]
            else cv2.INTER_LINEAR
        )
        tile = cv2.resize(
            display_slice,
            (tile_width, tile_height),
            interpolation=interpolation,
        )
        y_start = 0 if row == 0 else row_heights[0]
        x_start = 0 if column == 0 else column_widths[0]
        mosaic[y_start:y_start + tile_height, x_start:x_start + tile_width] = tile

    return mosaic


def get_tsnr_axial_indices(depth: int) -> tuple[int, int, int, int]:
    """Map the established normalized axial locations to valid indices."""
    depth = int(depth)
    if depth <= 0:
        raise ValueError("TSNR volume depth must be greater than zero.")
    return tuple(
        min(int(position * depth), depth - 1)
        for position in TSNR_AXIAL_POSITIONS
    )


def get_tsnr_display_slices(
    tsnr_volume: np.ndarray,
    display_max: float = 100.0,
) -> np.ndarray:
    """Return four clipped, y-flipped axial slices in TSNR display units."""
    volume = np.asarray(tsnr_volume)
    if volume.ndim != 3:
        raise ValueError(f"Expected a 3D TSNR volume, got shape {volume.shape}.")
    if any(d == 0 for d in volume.shape):
        raise ValueError(f"TSNR volume dimensions must be nonzero, got {volume.shape}.")
    if not np.isfinite(display_max) or display_max <= 0:
        raise ValueError("display_max must be a finite value greater than zero.")

    finite_volume = np.nan_to_num(
        volume.astype(np.float64, copy=False),
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )
    display_volume = np.flip(
        np.clip(finite_volume, 0.0, display_max),
        axis=1,
    )
    return display_volume[list(get_tsnr_axial_indices(volume.shape[0]))]


def calculate_mean_tsnr(tsnr_volume: np.ndarray) -> float:
    """Return the mean over finite positive voxels in the complete TSNR map."""
    volume = np.asarray(tsnr_volume, dtype=np.float64)
    if volume.ndim != 3:
        raise ValueError(f"Expected a 3D TSNR volume, got shape {volume.shape}.")
    # Without a brain mask, only finite positive voxels contribute so zero-valued
    # background and undefined statistics do not dominate the displayed mean.
    valid = np.isfinite(volume) & (volume > 0)
    if not np.any(valid):
        return 0.0
    return float(np.mean(volume[valid]))


def atomic_save_tsnr_image(
    mosaic: np.ndarray,
    output_path: str | PathLike[str],
    jpeg_quality: int = 95,
) -> None:
    """Encode a grayscale mosaic and atomically replace ``output_path``."""
    image = np.asarray(mosaic)
    if image.ndim != 2:
        raise ValueError(f"Expected a 2D grayscale mosaic, got shape {image.shape}.")
    if not 0 <= jpeg_quality <= 100:
        raise ValueError("jpeg_quality must be between 0 and 100.")

    if image.dtype != np.uint8:
        image = np.nan_to_num(image, nan=0.0, posinf=255.0, neginf=0.0)
        image = np.clip(image, 0, 255).astype(np.uint8)

    success, encoded = cv2.imencode(
        ".jpg",
        image,
        [int(cv2.IMWRITE_JPEG_QUALITY), int(jpeg_quality)],
    )
    if not success:
        raise OSError("OpenCV failed to encode the TSNR mosaic as JPEG.")

    final_path = os.path.abspath(os.fspath(output_path))
    output_dir = os.path.dirname(final_path)
    basename = os.path.basename(final_path)
    if not basename:
        raise ValueError("output_path must include a filename.")

    temporary_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_dir,
            prefix=f".{basename}.",
            suffix=".tmp.jpg",
            delete=False,
        ) as temporary_file:
            temporary_path = temporary_file.name
            temporary_file.write(encoded.tobytes())
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
        os.replace(temporary_path, final_path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            try:
                os.remove(temporary_path)
            except FileNotFoundError:
                pass

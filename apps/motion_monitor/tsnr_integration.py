"""Runtime orchestration helpers for motion-monitor TSNR integration."""

from __future__ import annotations

import json
import logging
import os
import re
import time
from collections.abc import Callable
from os import PathLike

import numpy as np
import SimpleITK as sitk

try:
    from .running_tsnr import (
        RunningTSNR,
        atomic_save_tsnr_image,
        create_tsnr_mosaic,
    )
except ImportError:  # Support execution from the motion_monitor application directory.
    from running_tsnr import (  # type: ignore[no-redef]
        RunningTSNR,
        atomic_save_tsnr_image,
        create_tsnr_mosaic,
    )


_POINTER_RE = re.compile(
    r"^(?P<protocol>.+)_volume_(?P<volume>\d{4})_group_(?P<group>\d{4})\.txt$"
)


class TSNRVolumeProcessor:
    """Attempt idempotent, non-blocking TSNR processing for pointer events."""

    def __init__(
        self,
        input_dir: str | PathLike[str],
        output_path: str | PathLike[str],
        *,
        min_samples: int = 20,
        display_max: float = 100.0,
        idle_retry_interval: float = 1.0,
        accumulator: RunningTSNR | None = None,
        mosaic_factory: Callable[..., np.ndarray] = create_tsnr_mosaic,
        image_saver: Callable[[np.ndarray, str | PathLike[str]], None] = atomic_save_tsnr_image,
        logger: logging.Logger | None = None,
    ) -> None:
        if min_samples < 2:
            raise ValueError("min_samples must be at least 2.")
        if idle_retry_interval <= 0:
            raise ValueError("idle_retry_interval must be greater than zero.")

        self.input_dir = os.path.abspath(os.fspath(input_dir))
        self.output_path = os.path.abspath(os.fspath(output_path))
        self.min_samples = min_samples
        self.display_max = display_max
        self.idle_retry_interval = idle_retry_interval
        self.accumulator = accumulator if accumulator is not None else RunningTSNR()
        self.pending_volumes: set[int] = set()
        self.processed_volumes: set[int] = set()
        self.image_ready_volumes: dict[int, np.ndarray] = {}
        self.motion_ready_volumes: dict[int, dict[str, int]] = {}
        self._pointer_paths: dict[int, str] = {}
        self._protocol_name: str | None = None
        self._last_idle_retry = 0.0
        self._mosaic_factory = mosaic_factory
        self._image_saver = image_saver
        self._logger = logger if logger is not None else logging.getLogger(__name__)

    def reset(self) -> None:
        """Reset per-acquisition TSNR state without deleting the last JPEG."""
        self.accumulator.reset()
        self.pending_volumes.clear()
        self.processed_volumes.clear()
        self.image_ready_volumes.clear()
        self.motion_ready_volumes.clear()
        self._pointer_paths.clear()
        self._protocol_name = None
        self._last_idle_retry = 0.0

    def handle_pointer(
        self,
        pointer_path: str | PathLike[str],
        expected_slice_count: int | None,
    ) -> None:
        """Record a stable group pointer and make one bounded retry pass."""
        pointer_path = os.path.abspath(os.fspath(pointer_path))
        match = _POINTER_RE.fullmatch(os.path.basename(pointer_path))
        if match is None:
            self._logger.warning(
                "TSNR: ignored pointer with unexpected name: %s",
                os.path.basename(pointer_path),
            )
            return

        protocol_name = match.group("protocol")
        volume_number = int(match.group("volume"))
        if self._protocol_name is None:
            self._protocol_name = protocol_name
        elif protocol_name != self._protocol_name:
            self._logger.warning(
                "TSNR: ignored pointer for protocol %s while processing %s",
                protocol_name,
                self._protocol_name,
            )
            return

        if volume_number not in self.processed_volumes:
            self.pending_volumes.add(volume_number)
            self._pointer_paths[volume_number] = pointer_path

        # The reference remains the unconditional initialization sample.
        if volume_number == 0:
            self.motion_ready_volumes.setdefault(
                0, {"registered": 0, "skipped": 0, "failed": 0}
            )

        self.retry_pending(expected_slice_count)

    def handle_registration_status(
        self,
        status_path: str | PathLike[str],
    ) -> bool:
        """Consume all newly finalized volumes from a durable status snapshot."""
        try:
            with open(status_path, "r", encoding="utf-8") as status_file:
                document = json.load(status_file)
            if document.get("schema_version") != 1:
                raise ValueError("unsupported or missing schema_version")
            volumes = document.get("volumes")
            if not isinstance(volumes, dict):
                raise ValueError("volumes must be an object")

            newly_ready: list[int] = []
            for volume_text, volume_record in volumes.items():
                if not isinstance(volume_text, str) or not volume_text.isdigit():
                    raise ValueError(f"invalid volume key: {volume_text!r}")
                volume_number = int(volume_text)
                if not isinstance(volume_record, dict):
                    raise ValueError(f"volume {volume_number} record must be an object")
                if volume_record.get("registration_finalized") is not True:
                    continue
                groups = volume_record.get("groups")
                if not isinstance(groups, dict) or not groups:
                    raise ValueError(f"volume {volume_number} groups must be nonempty")
                counts = {"registered": 0, "skipped": 0, "failed": 0}
                for group_text, group_record in groups.items():
                    if not isinstance(group_text, str) or not group_text.isdigit():
                        raise ValueError(
                            f"volume {volume_number} has invalid group key {group_text!r}"
                        )
                    if not isinstance(group_record, dict):
                        raise ValueError(
                            f"volume {volume_number} group {group_text} must be an object"
                        )
                    status = group_record.get("status")
                    if status not in counts:
                        raise ValueError(
                            f"volume {volume_number} group {group_text} has invalid status"
                        )
                    if status == "registered" and not isinstance(
                        group_record.get("transform_filename"), str
                    ):
                        raise ValueError(
                            f"volume {volume_number} group {group_text} lacks transform_filename"
                        )
                    counts[status] += 1

                if volume_number not in self.motion_ready_volumes:
                    self.motion_ready_volumes[volume_number] = counts
                    newly_ready.append(volume_number)
                    self._logger.info(
                        "TSNR: volume %d motion ready (%d registered, %d skipped, %d failed)",
                        volume_number,
                        counts["registered"],
                        counts["skipped"],
                        counts["failed"],
                    )
        except Exception as error:
            self._logger.warning("TSNR: registration status unavailable: %s", error)
            return False

        for volume_number in sorted(newly_ready):
            self._maybe_process_volume(volume_number)
        return True

    def retry_pending(self, expected_slice_count: int | None) -> None:
        """Attempt each pending volume once, oldest volume number first."""
        for volume_number in sorted(self.pending_volumes):
            if volume_number in self.processed_volumes:
                self.pending_volumes.discard(volume_number)
                continue

            if volume_number not in self.image_ready_volumes:
                try:
                    volume = self._load_volume(volume_number, expected_slice_count)
                    if volume is None:
                        continue
                    if (
                        self.accumulator.mean is not None
                        and volume.shape != self.accumulator.mean.shape
                    ):
                        raise ValueError(
                            f"assembled shape {volume.shape} does not match reference "
                            f"shape {self.accumulator.mean.shape}"
                        )
                    self.image_ready_volumes[volume_number] = volume
                    self._logger.info("TSNR: volume %d image ready", volume_number)
                except Exception as error:
                    self._logger.warning(
                        "TSNR: volume %d remains pending: %s", volume_number, error
                    )
                    continue

            self._maybe_process_volume(volume_number)

    def _maybe_process_volume(self, volume_number: int) -> bool:
        """Accumulate exactly once after image and registration are both ready."""
        if volume_number in self.processed_volumes:
            return False
        volume = self.image_ready_volumes.get(volume_number)
        if volume is None or volume_number not in self.motion_ready_volumes:
            return False

        try:
            self.accumulator.update(volume)
        except Exception as error:
            self._logger.warning(
                "TSNR: volume %d remains pending: %s", volume_number, error
            )
            return False

        # Success ordering is deliberate: update first, then mark processed.
        self.processed_volumes.add(volume_number)
        self.pending_volumes.discard(volume_number)
        self._pointer_paths.pop(volume_number, None)
        self.image_ready_volumes.pop(volume_number, None)
        item_name = "reference" if volume_number == 0 else f"volume {volume_number}"
        self._logger.info(
            "TSNR: %s accumulated (n=%d)", item_name, self.accumulator.count
        )

        if volume_number > 0:
            try:
                self._retire_volume_files(volume_number)
            except Exception as error:
                # Retirement is secondary housekeeping. Never make a successfully
                # accumulated volume eligible for tSNR processing again.
                self._logger.warning(
                    "TSNR: volume %d retirement failed: %s", volume_number, error
                )

        # Display errors must not make an accumulated volume eligible again.
        self._update_display(volume.shape[1:])
        return True

    def _retire_volume_files(self, volume_number: int) -> None:
        """Move one consumed non-reference volume out of the acquisition root."""
        if volume_number <= 0:
            return
        if self._protocol_name is None:
            self._logger.warning(
                "TSNR: volume %d retirement skipped (protocol unavailable)",
                volume_number,
            )
            return

        escaped_protocol = re.escape(self._protocol_name)
        volume_text = f"{volume_number:04d}"
        pointer_pattern = re.compile(
            rf"^{escaped_protocol}_volume_{volume_text}_group_(?P<index>\d{{4}})\.txt$"
        )
        header_pattern = re.compile(
            rf"^{escaped_protocol}_volume_{volume_text}_slice_(?P<index>\d{{4}})\.nhdr$"
        )
        raw_pattern = re.compile(
            rf"^{escaped_protocol}_volume_{volume_text}_slice_(?P<index>\d{{4}})\.raw$"
        )

        try:
            names = os.listdir(self.input_dir)
        except OSError as error:
            self._logger.warning(
                "TSNR: volume %d retirement scan failed: %s", volume_number, error
            )
            return

        pointers: list[str] = []
        headers: dict[int, str] = {}
        raw_files: dict[int, str] = {}
        for name in names:
            source_path = os.path.join(self.input_dir, name)
            if not os.path.isfile(source_path):
                continue
            pointer_match = pointer_pattern.fullmatch(name)
            if pointer_match is not None:
                pointers.append(name)
                continue
            header_match = header_pattern.fullmatch(name)
            if header_match is not None:
                headers[int(header_match.group("index"))] = name
                continue
            raw_match = raw_pattern.fullmatch(name)
            if raw_match is not None:
                raw_files[int(raw_match.group("index"))] = name

        complete_indices = sorted(set(headers) & set(raw_files))
        if not pointers and not complete_indices:
            self._logger.warning(
                "TSNR: volume %d retirement found no complete image-work files",
                volume_number,
            )
            return

        processed_dir = os.path.join(
            self.input_dir, f"processed_{self._protocol_name}"
        )
        try:
            os.makedirs(processed_dir, exist_ok=True)
        except OSError as error:
            self._logger.warning(
                "TSNR: volume %d retirement directory unavailable: %s",
                volume_number,
                error,
            )
            return

        moved = {"pointer": 0, "header": 0, "raw": 0}

        def destination_exists(name: str) -> bool:
            destination_path = os.path.join(processed_dir, name)
            if not os.path.lexists(destination_path):
                return False
            self._logger.warning(
                "TSNR: retirement collision for %s; source left active", name
            )
            return True

        def move_file(name: str, file_type: str) -> bool:
            source_path = os.path.join(self.input_dir, name)
            destination_path = os.path.join(processed_dir, name)
            if destination_exists(name):
                return False
            try:
                os.rename(source_path, destination_path)
            except OSError as error:
                self._logger.warning(
                    "TSNR: failed to retire %s: %s", name, error
                )
                return False
            moved[file_type] += 1
            return True

        # Publish each detached NRRD into the processed directory raw-first so
        # a moved header never intentionally precedes its data file. Separate
        # renames are not transactional, so roll the raw file back when the
        # associated header move fails and rollback remains possible.
        for slice_index in complete_indices:
            header_name = headers[slice_index]
            raw_name = raw_files[slice_index]
            if destination_exists(header_name) or destination_exists(raw_name):
                continue
            if not move_file(raw_name, "raw"):
                continue
            if move_file(header_name, "header"):
                continue

            moved_raw_path = os.path.join(processed_dir, raw_name)
            original_raw_path = os.path.join(self.input_dir, raw_name)
            if os.path.lexists(original_raw_path):
                self._logger.warning(
                    "TSNR: could not roll back retired raw file %s; source exists",
                    raw_name,
                )
                continue
            try:
                os.rename(moved_raw_path, original_raw_path)
            except OSError as error:
                self._logger.warning(
                    "TSNR: failed to roll back retired raw file %s: %s",
                    raw_name,
                    error,
                )
            else:
                moved["raw"] -= 1

        for missing_index in sorted(set(headers) ^ set(raw_files)):
            unmatched_name = headers.get(missing_index) or raw_files[missing_index]
            self._logger.warning(
                "TSNR: incomplete detached NRRD pair left active: %s", unmatched_name
            )

        # Pointer publication into the archive comes after its image pairs.
        for pointer_name in sorted(pointers):
            move_file(pointer_name, "pointer")

        self._logger.info(
            "Retired volume %d: %d pointers, %d headers, %d raw files -> %s/",
            volume_number,
            moved["pointer"],
            moved["header"],
            moved["raw"],
            os.path.basename(processed_dir),
        )

    def retry_pending_if_due(
        self,
        expected_slice_count: int | None,
        *,
        now: float | None = None,
    ) -> bool:
        """Perform at most one retry pass when the idle throttle is due."""
        if not self.pending_volumes:
            return False
        current_time = time.monotonic() if now is None else now
        if current_time - self._last_idle_retry < self.idle_retry_interval:
            return False
        self._last_idle_retry = current_time
        self.retry_pending(expected_slice_count)
        return True

    def _load_volume(
        self,
        volume_number: int,
        expected_slice_count: int | None,
    ) -> np.ndarray | None:
        if volume_number == 0:
            return self._load_reference_volume()
        if 0 not in self.processed_volumes:
            self._logger.info(
                "TSNR: volume %d pending (reference volume not yet accumulated)",
                volume_number,
            )
            return None
        if expected_slice_count is None or expected_slice_count <= 0:
            self._logger.info(
                "TSNR: volume %d pending (slice count metadata unavailable)",
                volume_number,
            )
            return None
        return self._load_slice_volume(volume_number, expected_slice_count)

    def _load_reference_volume(self) -> np.ndarray | None:
        pointer_path = self._pointer_paths.get(0)
        if pointer_path is None:
            self._logger.info("TSNR: volume 0 pending (reference pointer unavailable)")
            return None

        try:
            with open(pointer_path, "r", encoding="utf-8") as pointer_file:
                listed_names = [line.strip() for line in pointer_file if line.strip()]
        except OSError as error:
            raise OSError(f"failed to read reference pointer: {error}") from error

        if len(listed_names) != 1:
            raise ValueError(
                f"reference pointer must list exactly one header; found {len(listed_names)}"
            )
        listed_name = listed_names[0]
        if os.path.basename(listed_name) != listed_name:
            raise ValueError("reference pointer must contain a header basename")
        expected_prefix = f"{self._protocol_name}_volume_0000_"
        if (
            not listed_name.startswith(expected_prefix)
            or not listed_name.endswith(".nhdr")
            or listed_name.endswith("_upsampled.nhdr")
        ):
            raise ValueError(f"unexpected reference header name: {listed_name}")

        reference_path = os.path.join(self.input_dir, listed_name)
        reference_image = sitk.ReadImage(reference_path)
        reference_array = sitk.GetArrayFromImage(reference_image)
        if reference_array.ndim != 3:
            raise ValueError(
                f"reference must be a 3D volume; got shape {reference_array.shape}"
            )
        return reference_array

    def _load_slice_volume(
        self,
        volume_number: int,
        expected_slice_count: int,
    ) -> np.ndarray | None:
        assert self._protocol_name is not None
        prefix = f"{self._protocol_name}_volume_{volume_number:04d}_slice_"
        suffix = ".nhdr"
        slice_paths_by_index: dict[int, list[str]] = {}
        unexpected_names: list[str] = []

        for name in os.listdir(self.input_dir):
            if not name.startswith(prefix) or not name.endswith(suffix):
                continue
            index_text = name[len(prefix):-len(suffix)]
            if not re.fullmatch(r"\d{4}", index_text):
                unexpected_names.append(name)
                continue
            slice_paths_by_index.setdefault(int(index_text), []).append(
                os.path.join(self.input_dir, name)
            )

        expected_indices = set(range(expected_slice_count))
        discovered_indices = set(slice_paths_by_index)
        duplicate_indices = sorted(
            index for index, paths in slice_paths_by_index.items() if len(paths) != 1
        )
        if (
            discovered_indices != expected_indices
            or duplicate_indices
            or unexpected_names
        ):
            missing = sorted(expected_indices - discovered_indices)
            unexpected = sorted(discovered_indices - expected_indices)
            self._logger.info(
                "TSNR: volume %d pending (%d/%d slices; missing=%d; "
                "unexpected=%d; duplicates=%d; malformed=%d)",
                volume_number,
                len(discovered_indices & expected_indices),
                expected_slice_count,
                len(missing),
                len(unexpected),
                len(duplicate_indices),
                len(unexpected_names),
            )
            return None

        loaded: list[tuple[float, int, np.ndarray]] = []
        in_plane_shape: tuple[int, int] | None = None
        in_plane_spacing: tuple[float, float] | None = None
        in_plane_direction: np.ndarray | None = None

        for slice_index in range(expected_slice_count):
            slice_path = slice_paths_by_index[slice_index][0]
            try:
                image = sitk.ReadImage(slice_path)
            except Exception as error:
                raise OSError(
                    f"failed to read slice {slice_index}: {error}"
                ) from error

            array = sitk.GetArrayFromImage(image)
            if image.GetDimension() != 3 or array.ndim != 3 or array.shape[0] != 1:
                raise ValueError(
                    f"slice {slice_index} is not one 3D single-plane image: "
                    f"dimension={image.GetDimension()}, shape={array.shape}"
                )
            plane = array[0]
            spacing = tuple(float(value) for value in image.GetSpacing()[:2])
            direction = np.asarray(image.GetDirection(), dtype=np.float64).reshape(3, 3)
            plane_direction = direction[:, :2]
            if not np.isfinite(spacing).all() or any(value <= 0 for value in spacing):
                raise ValueError(
                    f"slice {slice_index} has invalid in-plane spacing {spacing}"
                )
            if not np.isfinite(direction).all():
                raise ValueError(f"slice {slice_index} has non-finite direction values")
            if not np.allclose(
                plane_direction.T @ plane_direction,
                np.eye(2),
                rtol=1e-5,
                atol=1e-6,
            ):
                raise ValueError(
                    f"slice {slice_index} has invalid in-plane direction geometry"
                )

            if in_plane_shape is None:
                in_plane_shape = plane.shape
                in_plane_spacing = spacing
                in_plane_direction = plane_direction
            else:
                if plane.shape != in_plane_shape:
                    raise ValueError(
                        f"slice {slice_index} shape {plane.shape} does not match "
                        f"{in_plane_shape}"
                    )
                if not np.allclose(spacing, in_plane_spacing, rtol=1e-5, atol=1e-6):
                    raise ValueError(
                        f"slice {slice_index} in-plane spacing {spacing} does not "
                        f"match {in_plane_spacing}"
                    )
                if not np.allclose(
                    plane_direction, in_plane_direction, rtol=1e-5, atol=1e-6
                ):
                    raise ValueError(
                        f"slice {slice_index} in-plane direction is inconsistent"
                    )

            center_index = [(size - 1) / 2 for size in image.GetSize()]
            center_z = float(
                image.TransformContinuousIndexToPhysicalPoint(center_index)[2]
            )
            if not np.isfinite(center_z):
                raise ValueError(f"slice {slice_index} has a non-finite center z")
            loaded.append((center_z, slice_index, plane))

        loaded.sort(key=lambda item: (item[0], item[1]))
        return np.stack([item[2] for item in loaded], axis=0)

    def _update_display(self, output_shape: tuple[int, int]) -> None:
        if self.accumulator.count < self.min_samples:
            return
        try:
            tsnr_volume = self.accumulator.get_tsnr()
            mosaic = self._mosaic_factory(
                tsnr_volume,
                output_shape=output_shape,
                display_max=self.display_max,
            )
            self._image_saver(mosaic, self.output_path)
        except Exception as error:
            self._logger.warning("TSNR: display update failed: %s", error)
            return
        self._logger.info("TSNR: updated display (n=%d)", self.accumulator.count)

#!/usr/bin/env bash
# Mount an ext* chroot image read-only.
# Usage: ./mount_chroot_ro.sh IMAGE.img
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Syntax: $0 IMAGE.img" >&2
  exit 2
fi

IMG="$1"
[[ -f "$IMG" ]] || { echo "Image not found: $IMG" >&2; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 1; }; }
need sudo; need losetup; need mount; need mktemp

MNT="$(mktemp -d -t chroot_ro.XXXXXXXX)"
# Create a read-only loop device
LOOP="$(sudo losetup --read-only --find --show "$IMG")"

# Mount read-only (noexec,nodev,nosuid for safety)
sudo mount -o ro,noexec,nodev,nosuid "$LOOP" "$MNT"

echo "Mounted read-only:"
echo "  Image     : $IMG"
echo "  Loop dev  : $LOOP"
echo "  Mountpoint: $MNT"
echo
echo "When done, unmount with:"
echo "  sudo umount \"$MNT\" && sudo losetup -d \"$LOOP\" && rmdir \"$MNT\""

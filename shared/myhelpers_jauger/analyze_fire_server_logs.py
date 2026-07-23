import re
import argparse
from datetime import datetime, timezone
import matplotlib.pyplot as plt
import numpy as np

def extract_slices_per_volume(logfile_path):
    """
    Extract 'nSlices per volume' from logfile.
    """
    pattern = re.compile(r"nSlices per volume\s*=\s*(\d+)")
    with open(logfile_path, "r") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                return int(match.group(1))

    raise ValueError("Could not find 'nSlices per volume' in log file.")


def extract_receive_timestamps(logfile_path):
    """
    Parse logfile and extract datetime objects for matching lines.
    """
    receive_timestamps = []
    search_string = "Received MRD_MESSAGE_ISMRMRD_IMAGE (1022)"
    pattern = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3})")
    with open(logfile_path, "r") as f:
        for line in f:
            if search_string in line:
                match = pattern.match(line)
                if match:
                    dt_str = match.group(1)
                    dt = datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S,%f")
                    dt = dt.replace(tzinfo=timezone.utc)
                    receive_timestamps.append(dt)

    return receive_timestamps


def extract_acquisition_times(logfile_path):
    """
    Parse logfile and extract scanner acquisition timestamps (ticks since midnight, at 400 Hz).
    """
    acquisition_times = []
    pattern = re.compile(r"Image acquisition time\s*:\s*(\d+)")
    with open(logfile_path, "r") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                acquisition_times.append(int(match.group(1)))

    return acquisition_times


def plot_deltas_with_volumes(deltas, slices_per_volume, title, ylabel):
    plt.figure(figsize=(14, 5))
    # Plot deltas
    plt.plot(
        deltas,
        marker='o',
        linestyle='-',
        markersize=5,
        alpha=0.6
    )

    # Mark out volumes
    volume_indices = [k * slices_per_volume - 1 for k in range(1, (len(deltas) // slices_per_volume) + 1)]
    for idx in volume_indices:
        plt.axvline(
            x=idx,
            linestyle='--',
            linewidth=0.9,
            color='red',
            alpha=0.6
        )

    plt.xlim(0, len(deltas))
    # plt.xlim(50, 60)
    plt.ylim(0.065, 0.095)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Image index i", fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.grid(True)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Build a registration-grade multi-resolution pyramid.")
    parser.add_argument("--logfile", required=True, help="Log file from fire-server.")
    args = parser.parse_args()

    logfile = args.logfile
    receive_timestamps = extract_receive_timestamps(logfile)
    acquisition_times = extract_acquisition_times(logfile)

    if len(receive_timestamps) < 2:
        print("Not enough receive timestamps found.")
        return
    
    receive_deltas = np.array([
        (receive_timestamps[i] - receive_timestamps[i - 1]).total_seconds()
        for i in range(1, len(receive_timestamps))
    ])

    acquisition_deltas = None
    if len(acquisition_times) == 0:
        print("\nNo acquisition timestamps found. Skipping acquisition timing analysis.")
    elif len(acquisition_times) != len(receive_timestamps):
        print(
            f"\nWarning: Found {len(acquisition_times)} acquisition timestamps "
            f"but {len(receive_timestamps)} receive timestamps. "
            "Skipping acquisition timing analysis."
        )

    else:
        acquisition_deltas = np.diff(acquisition_times) / 400   # ticks since midnight conversion to ms


    slices_per_volume = extract_slices_per_volume(logfile)
    num_slices = len(receive_timestamps)
    num_volumes = num_slices // slices_per_volume
    remainder = num_slices % slices_per_volume
    print(f"Total slices: {num_slices}")
    print(f"Slices per volume: {slices_per_volume}")
    print(f"Total full volumes: {num_volumes}")
    if remainder != 0:
        print(f"Warning: {remainder} slices do not complete a full volume.")


    print("\nReceive timing deltas:")
    print(f"  Mean : {np.mean(receive_deltas):.6f} s")
    print(f"  Std  : {np.std(receive_deltas):.6f} s")
    print(f"  Min  : {np.min(receive_deltas):.6f} s")
    print(f"  Max  : {np.max(receive_deltas):.6f} s")

    if acquisition_deltas is not None:
        print("\nAcquisition timing deltas:")
        print(f"  Mean : {np.mean(acquisition_deltas):.6f} s")
        print(f"  Std  : {np.std(acquisition_deltas):.6f} s")
        print(f"  Min  : {np.min(acquisition_deltas):.6f} s")
        print(f"  Max  : {np.max(acquisition_deltas):.6f} s")


    plot_deltas_with_volumes(receive_deltas, slices_per_volume, "Time delay between receipt of image slice i and i+1", "Delta time (sec)")
    
    if acquisition_deltas is not None:
        plot_deltas_with_volumes(acquisition_deltas, slices_per_volume, "Time delay between acquisition of image slice i and i+1", "Delta time (sec)")


if __name__ == "__main__":
    main()
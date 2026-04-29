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


def extract_timestamps(logfile_path, search_string):
    """
    Parse logfile and extract datetime objects for matching lines.
    """
    timestamps = []
    pattern = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3})")
    with open(logfile_path, "r") as f:
        for line in f:
            if search_string in line:
                match = pattern.match(line)
                if match:
                    dt_str = match.group(1)
                    dt = datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S,%f")
                    dt = dt.replace(tzinfo=timezone.utc)
                    timestamps.append(dt)

    return timestamps


def compute_deltas(timestamps):
    return [timestamps[i] - timestamps[i - 1] for i in range(1, len(timestamps))]


def plot_deltas_with_volumes(deltas, slices_per_volume):
    plt.figure()
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
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Image index i", fontsize=14)
    plt.ylabel("Delta time (sec)", fontsize=14)
    plt.title("Time delay between receipt of image slice i and i+1", fontsize=16)
    plt.grid(True)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Build a registration-grade multi-resolution pyramid.")
    parser.add_argument("--logfile", required=True, help="Log file from fire-server.")
    args = parser.parse_args()

    # logfile = "/home/jauger/GitHubRepos/python-fire-server-jauger/received_data/savedData_20260428T141617_func-bold_task-rest480_run-01_SLIMMON/python-fire-server-jauger_20260428_141606.log"
    logfile = args.logfile
    search_string = "Received MRD_MESSAGE_ISMRMRD_IMAGE (1022)"

    slices_per_volume = extract_slices_per_volume(logfile)

    timestamps = extract_timestamps(logfile, search_string)
    if len(timestamps) < 2:
        print("Not enough matching entries found.")
        return

    num_slices = len(timestamps)
    num_volumes = num_slices // slices_per_volume
    remainder = num_slices % slices_per_volume
    print(f"Total slices: {num_slices}")
    print(f"Slices per volume: {slices_per_volume}")
    print(f"Total full volumes: {num_volumes}")
    if remainder != 0:
        print(f"Warning: {remainder} slices do not complete a full volume.")

    unix_times = [dt.timestamp() for dt in timestamps]
    deltas = compute_deltas(unix_times)

    plot_deltas_with_volumes(deltas, slices_per_volume)

    print(f"Mean time delay between image slices: {np.mean(deltas)}")
    print(f"Std dev time delay between image slices: {np.std(deltas)}")


if __name__ == "__main__":
    main()
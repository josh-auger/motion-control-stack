#!/usr/bin/env python3

# Title: process_queue_directory.py

# Description:
# Monitors a specified input directory for image files being written into the directory. New image files are identified
# and sent to sms-mi-reg for mutual information-based image registration.

# Created on: March 2025
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import os
import re
import sys
import csv
import glob
import time
import logging
import argparse
import subprocess
import SimpleITK as sitk
import numpy as np
import pandas as pd
from datetime import datetime
import json


def setup_logging(log_dir):
    """Configure logging to save logs with a timestamped filename in /working/."""
    log_filename = os.path.join(f"{log_dir}/log_local_queue_processor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ]
    )
    logging.info(f"Logging will be saved to : {log_filename}")


def reset_logging(log_dir):
    """Reset logging configuration to save logs to a new logfile with a fresh timestamp."""
    log_filename = os.path.join(f"{log_dir}/log_local_queue_processor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    # Remove all existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Reconfigure logging with new file
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ]
    )
    logging.info(f"Log reset. Logging will be saved to : {log_filename}")


def read_nrrd_image(filepath):
    """Read a NRRD image file and return the SimpleITK image object."""
    reader = sitk.ImageFileReader()
    reader.SetFileName(filepath)
    try:
        image = reader.Execute()
        logging.info(f"Read in {filepath}")
        return image
    except RuntimeError as e:
        logging.error(f"Failed to read {filepath}\nError: {str(e)}")
        logging.error("Terminating program and killing local-queue-processor container.")
        sys.exit(1)     # Exit program with non-zero status code


def create_identityTransformFile(reference_volume_filepath, indexcount, volcount, groupcount):
    """Create an identity transform file with center-of-rotation at the spatial center of the reference volume."""
    logging.info(f"Creating identity transform about center of reference volume {volcount:04d}...")
    transform_filename = f"alignTransform_{indexcount:04d}_{volcount:04d}-{groupcount:04d}_identity.tfm"
    transform_filepath = os.path.join("/working/", transform_filename)

    refImage = read_nrrd_image(reference_volume_filepath)

    image_origin = refImage.TransformContinuousIndexToPhysicalPoint([0 for index in refImage.GetSize()])
    logging.info('\timage_origin is : ' + str(image_origin))
    image_center = refImage.TransformContinuousIndexToPhysicalPoint([(index - 1) / 2.0 for index in refImage.GetSize()])
    logging.info('\timage_center is : ' + str(image_center))
    image_far = refImage.TransformContinuousIndexToPhysicalPoint([(index - 1) for index in refImage.GetSize()])
    logging.info('\timage_far is : ' + str(image_far))

    transform = sitk.VersorRigid3DTransform()
    transform.SetIdentity()
    transform.SetCenter(image_center)  # set center of rotation at reference volume center
    sitk.WriteTransform(transform, transform_filepath)
    logging.info(f"\tIdentity transform saved as : {transform_filepath}")
    return transform_filepath


def select_input_transform(identityTransform_filepath, counter, reference_volume_flag):
    """Identify prior alignment transform file (.tfm) as input initialization transform in next registration call."""
    if reference_volume_flag == 0:
        logging.info(f"Still calibrating reference volume. Input transform for registration : {identityTransform_filepath}")
        return identityTransform_filepath

    prior_transform_filepath = f"/working/alignTransform_{counter - 1:04d}_*.tfm"   # ignore volcount or groupcount, just use index counter
    matches = glob.glob(prior_transform_filepath)
    if matches:
        chosen_transform = max(matches, key=os.path.getmtime)
    else:
        chosen_transform = identityTransform_filepath

    # chosen_transform = identityTransform_filepath
    logging.info(f"Input transform for registration : {chosen_transform}")
    return chosen_transform


def run_MIregistration(reference_volume_filepath, target_filepaths, inputTransform_filepath, outputTransformLabel):
    """
    Compile all inputs for the sub-process run command to execute registration between the assigned reference volume and
    the specified target image slice(s).
    """
    # Ensure target_volume_filepaths is a list
    if isinstance(target_filepaths, str):
        target_filepaths = [target_filepaths]

    # Check that all required input files exist
    all_files = [(reference_volume_filepath, "Reference volume"), (inputTransform_filepath, "Input transform")]
    all_files += [(f, f"Target volume {i + 1}") for i, f in enumerate(target_filepaths)]

    for fpath, label in all_files:
        if not os.path.exists(fpath):
            logging.info(f"Registration FAILED. {label} file not found : {fpath}")
            return None

    # Call sms-mi-reg executable with local filepath inputs
    optimizer = "LN_BOBYQA"
    maxiterations = "1000"
    run_command = [
                      "/usr/src/moco/sms-mi-reg",
                      reference_volume_filepath,
                      inputTransform_filepath,
                      outputTransformLabel
                  ] + target_filepaths + [
                      "--optimizer", optimizer,
                      "--maxiter", maxiterations
                  ]

    logging.info(f"Registration run command : {run_command}")
    start_time = time.time()
    try:
        registration_process = subprocess.run(run_command, shell=False, text=True, capture_output=True, check=True)
        end_time = time.time()
        run_time = end_time - start_time
        logging.info(registration_process.stdout.replace('\n', '\n\t'))
        logging.info(f"Registration call elapsed runtime (sec) : {run_time:.10f}")
        return None
    except subprocess.CalledProcessError as e:
        logging.error(f"Registration FAILED : {e.stderr}")
        return None


def compose_transform_pair(transform1, transform2):
    """Calculate the composed transform between two given input transforms."""
    A0 = np.asarray(transform2.GetMatrix()).reshape(3, 3)
    c0 = np.asarray(transform2.GetCenter())
    t0 = np.asarray(transform2.GetTranslation())

    A1 = np.asarray(transform1.GetInverse().GetMatrix()).reshape(3, 3)
    c1 = np.asarray(transform1.GetInverse().GetCenter())
    t1 = np.asarray(transform1.GetInverse().GetTranslation())

    combined_mat = np.dot(A0,A1)
    combined_center = c1
    combined_translation = np.dot(A0, t1+c1-c0) + t0+c0-c1

    euler3d = sitk.Euler3DTransform()
    euler3d.SetCenter(combined_center)
    euler3d.SetTranslation(combined_translation)
    euler3d.SetMatrix(combined_mat.flatten())
    return euler3d


def calculate_displacement(euler3d_transform, radius=50):
    """Calculate framewise displacement from a Euler 3D transform (following Tisdall et al. 2012)"""
    # assumes Euler3D transform as input
    params = np.asarray(euler3d_transform.GetParameters())
    # logging.info(f"\tParameters (Euler3D) : {params}")

    # displacement calculation
    theta = np.abs(np.arccos(0.5 * (-1 + np.cos(params[0]) * np.cos(params[1]) + \
                                    np.cos(params[0]) * np.cos(params[2]) + \
                                    np.cos(params[1]) * np.cos(params[2]) + \
                                    np.sin(params[0]) * np.sin(params[1]) * np.sin(params[2]))))
    drot = radius * np.sqrt((1 - np.cos(theta)) ** 2 + np.sin(theta) ** 2)
    dtrans = np.linalg.norm(params[3:])
    displacement = drot + dtrans
    return displacement


def motion_table_to_dataframe(motion_table):
    """
    Convert the internal motion_table (list of dicts)
    into a pandas DataFrame suitable for plotting/export.
    """
    if not motion_table:
        return pd.DataFrame()

    df = pd.DataFrame(motion_table)
    column_order = [
        "reg_index",
        "X_rotation(rad)", "Y_rotation(rad)", "Z_rotation(rad)",
        "X_translation(mm)", "Y_translation(mm)", "Z_translation(mm)",
        "Displacement(mm)",
        "Cumulative_displacement(mm)",
        "Motion_flag",
    ]
    return df[column_order]


def convert_versor_to_euler(transform):
    """Convert a Versor Rigid 3D transform into an equivalent Euler 3D transform."""
    center = transform.GetCenter()
    translation = transform.GetTranslation()
    rotation_matrix = transform.GetMatrix()

    # Extract Euler angles (in radians) from the Versor rotation matrix
    R = np.array(rotation_matrix).reshape(3, 3)  # Convert to 3x3 matrix

    # Convert to Euler angles (ZYX convention)
    sy = np.sqrt(R[0, 0] ** 2 + R[1, 0] ** 2)
    singular = sy < 1e-6

    if not singular:
        x = np.arctan2(R[2, 1], R[2, 2])
        y = np.arctan2(-R[2, 0], sy)
        z = np.arctan2(R[1, 0], R[0, 0])
    else:
        x = np.arctan2(-R[1, 2], R[1, 1])
        y = np.arctan2(-R[2, 0], sy)
        z = 0

    # Create the Euler3DTransform
    euler_transform = sitk.Euler3DTransform()
    euler_transform.SetRotation(x, y, z)  # Angles are in radians
    euler_transform.SetTranslation(translation)
    euler_transform.SetCenter(center)
    return euler_transform


def check_reference_volume(transform_filepath, threshold=0.60):
    """
    Check framewise displacement of input alignment transform (from registration of new volume and provisional
    reference volume).
    Default acceptable motion threshold = 0.60 mm (25% of 2.4 mm voxel size of functional resting state sequence).
    """
    # Convert versor 3D rigid transform to Euler3D transform
    euler3d_transform = convert_versor_to_euler(sitk.ReadTransform(transform_filepath))
    displacement = calculate_displacement(euler3d_transform)
    logging.info(f"\tDisplacement : {displacement:.3f} mm (Threshold : {threshold:.3f} mm)")
    return 1 if displacement < threshold else 0


def make_test_transform(identity_transform_filepath, volcount):
    """
    Generate dummy transforms to test moco feedback.
    Set the range of values for Euler rotation angles and translations, then randomly generate 3 rotations and 3
    translations for the dummy transform.
    """
    identityTransform = sitk.ReadTransform(identity_transform_filepath)
    transform_center = identityTransform.GetCenter()

    test_transform = sitk.Euler3DTransform()
    test_transform.SetIdentity()
    test_transform.SetCenter(transform_center)

    # Generate 3 random rotations + 3 random translations
    rotation_range = (-0.25,0.25)       # Selection range of Euler rotations (rads)
    translation_range = (-10.0,10.0)    # Selection range of translations (mm)
    rotations = np.round(np.random.uniform(rotation_range[0], rotation_range[1], size=3), 4)    # round to 4 decimal places
    translations = np.round(np.random.uniform(translation_range[0], translation_range[1], size=3), 4)
    parameters = np.concatenate([rotations, translations])

    logging.info(f"\tRandom Euler3D test parameters for volume {volcount:04d} : {parameters}")
    test_transform.SetParameters(parameters)

    # Write out transform
    sitk.WriteTransform(test_transform, f"/working/testTransform_{volcount:04d}.tfm")
    return


def read_slice_timings_from_json(json_filepath):
    """Read series metadata file (.json) and extract the list of slice timings"""
    if not os.path.isfile(json_filepath):
        logging.error(f"JSON file not found: {json_filepath}")
        return {}

    with open(json_filepath, "r") as f:
        metadata = json.load(f)

    if "SliceTiming" not in metadata:
        logging.warning(f"'SliceTiming' field not found in {json_filepath}")
        return {}

    slice_timing_dict = metadata["SliceTiming"]
    # Convert to list ordered by slice index
    slice_timing = [slice_timing_dict[str(i)] if str(i) in slice_timing_dict else slice_timing_dict[i]
            for i in sorted(map(int, slice_timing_dict.keys()))]

    # Determine slice order by ascending acquisition time
    sorted_indices = sorted(range(len(slice_timing)), key=lambda i: slice_timing[i])
    logging.info(f"Extracted slice timings (ordered by index): {slice_timing}")
    logging.info(f"Slice acquisition order (ascending time): {sorted_indices}")
    return slice_timing


def monitor_directory(input_dir, moco_flag):
    """Monitor directory for new image files without deleting any."""
    # Initialization
    # -------------------------------------------
    if not os.path.exists(input_dir):
        logging.error(f"Directory '{input_dir}' does not exist.")
        return

    log_dir = input_dir
    moco_enabled = (moco_flag.lower() == "on") if isinstance(moco_flag, str) else False
    logging.info(f"Motion correction {'ENABLED' if moco_enabled else 'DISABLED'}.")

    VALID_EXTENSIONS = {'.json', '.txt', '.closeQ'}  # JDA: Read incoming metadata (.json) and group pointer files (.txt) only!
    logging.info(f"Monitoring directory [{VALID_EXTENSIONS}] : {input_dir} ...")

    # Reset all monitoring state variables
    # -------------------------------------------
    def reset_variables():
        return {
            "metadata_filepath": None,
            "reference_volume_filepath": None,
            "reference_volume_flag": 0,
            "itemcount": 0,     # running count of files seen (all types)
            "regcount": 0,      # running count of registrations run
            "volcount": 0,      # current image volume number
            "groupcount": 0,    # current slice group number
            "seen_files": set(),
            "previous_filesize": 0,
            "begintime": None,
            "slice_timings": None,
            "motion_table": [],
            "cumulative_displacement": 0.0,
            "motion_flag_count": 0
        }
    state = reset_variables()

    # Helper functions
    # -------------------------------------------
    def list_new_files():
        """Return list of valid, non-seen files sorted by modification time,
        ignoring files that disappear during the scan."""
        valid_files = []
        for f in os.listdir(input_dir):
            # Skip wrong extensions
            if os.path.splitext(f)[1] not in VALID_EXTENSIONS:
                continue
            # Skip files that have already been processed
            if f in state["seen_files"]:
                continue

            full_path = os.path.join(input_dir, f)
            try:    # attempt to ping for modification time
                mtime = os.path.getmtime(full_path)
            except FileNotFoundError:
                # File disappeared between os.listdir() and os.path.getmtime() (i.e. end of sequence consolidation)
                logging.info(f"Unable to ping file for mod time : {full_path}")
                continue

            valid_files.append((f, mtime))

        # Sort by mtime: oldest → newest
        valid_files.sort(key=lambda x: x[1])
        # Return only new files not yet processed
        return [f for (f, _) in valid_files]

    def pick_files_moco(new_files):
        """Moco: keep newest file only, discard older ones."""
        newest = max(
            new_files,
            key=lambda f: os.path.getmtime(os.path.join(input_dir, f))
        )
        for f in new_files:
            if f != newest:
                logging.info(f"Old file found. Skipping {f}")
                state["seen_files"].add(f)
        return [newest]

    def wait_for_complete_write(filepath, max_checks=200, delay=0.005):
        """Poll file size until it is stable for one check. True if stable, False if file disappeared."""
        checks = 0
        last_size = -1
        while checks < max_checks:
            if not os.path.exists(filepath):
                logging.warning(f"File no longer exists : {filepath}")  # file moved/deleted (i.e. end of sequence consolidation)
                return False
            size_now = os.path.getsize(filepath)
            if size_now == last_size:
                return True
            last_size = size_now
            checks += 1
            time.sleep(delay)
        # timed out waiting for stability, still proceed but warn
        logging.warning(f"File write did not stabilize after {max_checks} checks for {filepath}; proceeding anyway.")
        return True

    def process_metadata_file(filepath):
        """Handle incoming series metadata JSON."""
        if state["metadata_filepath"] is None:
            logging.info(f"Found series metadata file : {os.path.basename(filepath)}")
            state["metadata_filepath"] = filepath
        else:
            logging.info(f"Skipping {filepath} because metadata already loaded.")
        state["seen_files"].add(os.path.basename(filepath))

    def read_pointer_file(filepath):
        """Read listed filenames from pointer file."""
        with open(filepath, "r") as f:
            listed = [line.strip() for line in f if line.strip()]
        return [os.path.join(input_dir, name) for name in listed]

    def update_counters_from_pointer_file(pointer_filepath):
        """Extract volume count string and slice group count string from pointer filename."""
        basefilename = os.path.basename(pointer_filepath)
        match = re.search(r"volume_(\d{4})_group_(\d{4})", basefilename)
        if not match:
            raise ValueError(f"Filename does not match expected pattern: {basefilename}")
        state["volcount"] = int(match.group(1))
        state["groupcount"] = int(match.group(2))
        return

    def initialize_reference_volume(target_paths):
        """Set first batch as reference volume."""
        state["reference_volume_filepath"] = target_paths[0]
        logging.info(f"Provisional reference volume set to : {state['reference_volume_filepath']}")

        identity_transform_path = create_identityTransformFile(state["reference_volume_filepath"], state["regcount"], state["volcount"], state["groupcount"])

        slice_timings = read_slice_timings_from_json(state["metadata_filepath"])
        state["slice_timings"] = slice_timings

        if not slice_timings:
            logging.warning("No slice timings found in metadata.")
        else:
            logging.info(f"Number of slices per volume : {len(slice_timings)}")
            groups = np.unique(slice_timings)
            logging.info(f"Number of slice groups : {len(groups)}")
            logging.info(f"Number of slices per group : {len(slice_timings) / len(groups)}")

        return identity_transform_path

    def run_registration_batch(target_paths, identity_transform_path):
        """Run MI registration for one batch of image files within a single group pointer file."""
        state["regcount"] += 1
        logging.info(f"Running registration call {state['regcount']}")
        logging.info(f"Reference volume : {state['reference_volume_filepath']}")
        logging.info("Target items :")
        for idx, t in enumerate(target_paths, 1):
            logging.info(f"\t[{idx:02d}] {t}")

        input_transform = select_input_transform(identity_transform_path, state["regcount"], state["reference_volume_flag"])
        # input_transform = identity_transform_path

        output_string = f"{state['regcount']:04d}_{state['volcount']:04d}-{state['groupcount']:04d}"
        run_MIregistration(state["reference_volume_filepath"], target_paths, input_transform, output_string)

        # track_framewise_displacement(f"/working/alignTransform_{output_string}.tfm", input_transform, 0.6)

    def track_framewise_displacement(current_transform_filepath, prior_transform_filepath, motion_threshold):
        """Maintain cumulative ledger of calculated framewise displacements between transform pairs."""
        prior_transform = convert_versor_to_euler(sitk.ReadTransform(prior_transform_filepath))
        current_transform = convert_versor_to_euler(sitk.ReadTransform(current_transform_filepath))
        combined_transform = compose_transform_pair(prior_transform, current_transform)
        prior_params = prior_transform.GetParameters()
        current_params = current_transform.GetParameters()
        framewise_displacement = calculate_displacement(combined_transform)
        motion_flag = 1 if framewise_displacement > motion_threshold else 0

        state["motion_flag_count"] += motion_flag
        state["cumulative_displacement"] += framewise_displacement

        row = {
            "reg_index": state["regcount"],
            "X_rotation(rad)": current_params[0],
            "Y_rotation(rad)": current_params[1],
            "Z_rotation(rad)": current_params[2],
            "X_translation(mm)": current_params[3],
            "Y_translation(mm)": current_params[4],
            "Z_translation(mm)": current_params[5],
            "Displacement(mm)": framewise_displacement,
            "Cumulative_displacement(mm)": state["cumulative_displacement"],
            "Volume_index": state["volcount"],
            "Slice_group_index": state["groupcount"],
            "Motion_flag": motion_flag
        }
        state["motion_table"].append(row)

        def format_params(params, precision=4):
            return "(" + ", ".join(f"{p:.{precision}g}" for p in params) + ")"

        logging.info(f"=================================")
        logging.info(f"===== MOTION SUMMARY : {state['regcount']:04d} =====")
        logging.info(f"\tPrior parameters (Euler) : {format_params(prior_params, 4)}")
        logging.info(f"\tCurrent parameters (Euler) : {format_params(current_params, 4)}")
        logging.info(f"\tFramewise displacement (mm) : {framewise_displacement:04f}")
        logging.info(f"\tCumulative displacement (mm) : {state['cumulative_displacement']:04f}")
        logging.info(f"\tCumulative motion flags : {state['motion_flag_count']}")
        logging.info(f"\tCurrent volume count : {state['volcount']} (slice group {state['groupcount']})")
        logging.info(f"=================================")

    def maybe_update_reference(pointer_filepath):
        """Check reference volume transform and update if needed."""
        if state["reference_volume_flag"] == 1:
            return

        logging.info("Calibrating reference volume...")
        tfm = f"/working/alignTransform_{state['regcount']:04d}.tfm"
        if not os.path.exists(tfm):
            logging.warning(f"\tTransform file not found: {tfm}")
            logging.warning("\tSkipping reference calibration. Volume 0000 as default reference.")
            state["reference_volume_flag"] = 1
            return

        flag = check_reference_volume(tfm)
        if flag == 0:
            state["reference_volume_filepath"] = pointer_filepath
            logging.info("\tReference calibration failed.")
            logging.info(f"\tAssigning new provisional reference volume : {pointer_filepath}")
        else:
            logging.info("\tReference volume calibration successful!")
            state["reference_volume_flag"] = 1

    def handle_reset_trigger(filepath):
        """Handle close trigger file."""
        logging.info(f"Reset trigger detected : {os.path.basename(filepath)}")
        try:
            os.remove(filepath)
        except Exception as e:
            logging.error(f"Failed to delete reset trigger file {filepath}: {e}")

        # # Export motion table BEFORE wiping state
        # export_motion_table_csv(input_dir)

        # Reset all state
        nonlocal_state = reset_variables()
        state.update(nonlocal_state)
        logging.info("\n\n---- Local-queue-processor reset ----")
        return


    # Main monitoring loop
    # =====================================================================
    while True:
        new_files = list_new_files()
        if not new_files:
            time.sleep(0.005)
            continue

        if new_files:
            # Start timer of new session
            if state["begintime"] is None:
                reset_logging(log_dir)
                state["begintime"] = time.time()
                logging.info(f"Started monitoring at : {datetime.now()}")

            logging.info(f"Found {len(new_files)} new file(s) to process")

            # Determine if moco is "on" and need to "pick newest only" logic
            apply_moco_filter = moco_enabled and state["metadata_filepath"] is not None and state["reference_volume_filepath"] is not None
            if apply_moco_filter:
                new_files = pick_files_moco(new_files)  # Keep only newest file

            # Process each new file
            for fname in new_files:
                state["itemcount"] += 1
                new_filepath = os.path.join(input_dir, fname)
                ext = os.path.splitext(fname)[1]

                # Handle CLOSE trigger file(s)
                if ext == ".closeQ":
                    handle_reset_trigger(new_filepath)
                    continue

                # Handle series metadata file
                if ext == ".json":
                    process_metadata_file(new_filepath)
                    continue

                # Handle pointer file
                if ext == ".txt":
                    start_time = time.time()
                    if not wait_for_complete_write(new_filepath):
                        state["seen_files"].add(fname)
                        continue

                    # Populate pending target paths with the listed image files in the pointer file
                    update_counters_from_pointer_file(new_filepath)
                    target_paths = read_pointer_file(new_filepath)
                    logging.info(f"Added {len(target_paths)} file(s) from {fname} to pending targets.")

                    # First pointer file --> reference volume
                    if state["reference_volume_filepath"] is None:
                        identity_transform_path = initialize_reference_volume(target_paths)
                    else:
                        run_registration_batch(target_paths, identity_transform_path)
                        # Check reference volume status
                        maybe_update_reference(new_filepath)

                    # Mark file as processed
                    state["seen_files"].add(fname)

                    item_process_time = time.time() - start_time
                    total_process_time = time.time() - state["begintime"]
                    avg_item_process_time = total_process_time / state['itemcount']
                    # Logging block
                    logging.info(f"===================================")
                    logging.info(f"===== PROGRESS SUMMARY : {state['itemcount']:04d} =====")
                    logging.info(f"Current volume number : {state['volcount']} (slice group {state['groupcount']})")
                    logging.info(f"Processed item {state['itemcount']:04d} in {item_process_time:.3f} (s)")
                    logging.info(f"Total processed files : {state['itemcount']}")
                    logging.info(f"Total registration calls : {state['regcount']}")
                    logging.info(f"Total elapsed time (sec) : {total_process_time:.3f}")
                    logging.info(f"Average item processing time (sec) : {avg_item_process_time:.3f}")
                    logging.info(f"===================================")


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("input_directory", nargs="?", help="Path to image files")
    parser.add_argument("--moco", type=str, choices=["on", "off"])
    args = parser.parse_args()

    env_moco = os.environ.get("MOCO_FLAG", "off")
    moco_flag = args.moco if args.moco is not None else env_moco

    setup_logging(args.input_directory)
    monitor_directory(args.input_directory, moco_flag)

if __name__ == "__main__":
    main()
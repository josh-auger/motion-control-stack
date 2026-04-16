#!/usr/bin/env python3
"""
Title: run_nipype_docker_tsnr.py

Description:
Run a NiPype-based tSNR workflow inside the NiPype Docker container.
Outputs tSNR, mean, and std maps, plus a CSV table of values.

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Date of creation: October 17, 2025
"""

import os
import sys
import argparse
import subprocess
import datetime

def run_nipype_tsnr(infile, outdir=None, container="nipype/nipype:latest", workflow_script="nipype_tsnr_workflow.py"):
    """Run the NiPype tSNR workflow inside a Docker container."""

    # Verify input file
    if not os.path.exists(infile):
        sys.exit(f"Error : Input file not found: {infile}")
    infile = os.path.abspath(infile)

    # Use same directory as input if no outdir is given
    if outdir is None:
        outdir = os.path.dirname(os.path.abspath(infile))
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    # Current working directory (PWD)
    pwd = os.getcwd()
    workflow_script_host = os.path.join(pwd, workflow_script)
    if not os.path.exists(workflow_script_host):
        sys.exit(f"[ERROR] Workflow script not found in current directory: {workflow_script_host}")

    # Determine the parent directory to mount
    parent_mount = os.path.commonpath([infile, outdir])
    if not parent_mount:
        parent_mount = os.path.dirname(infile)

    # Prepare container paths
    workflow_script_container = f"/workflow/{workflow_script}"
    infile_container = f"/data{infile.replace(parent_mount, '')}"
    outdir_container = f"/data{outdir.replace(parent_mount, '')}"

    timestamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    log_file = os.path.join(outdir, f"nipype_tsnr_log_{timestamp}.txt")

    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{pwd}:/workflow",  # Mount PWD for workflow script
        "-v", f"{parent_mount}:/data",  # Mount input/output file directory
        container,
        "python", workflow_script_container,
        "--infile", infile_container,
        "--outdir", outdir_container
    ]

    print(f"Running NiPype tSNR workflow in Docker container:")
    print(" ".join(docker_cmd))
    print(f"Logging output to: {log_file}")

    with open(log_file, "w") as log:
        result = subprocess.run(docker_cmd, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode == 0:
        print("[COMPLETED] NiPype tSNR workflow completed successfully.")
        print(f"Results saved to : {outdir}")
    else:
        print("[FAILED] NiPype tSNR workflow failed. Check log file for details:")
        print(log_file)

def main():
    parser = argparse.ArgumentParser(description="Run NiPype tSNR workflow inside nipype/nipype Docker container.")
    parser.add_argument("--infile", required=True, help="Input 4D NIfTI (.nii.gz)")
    parser.add_argument("--outdir", required=False, help="Output directory for results")
    parser.add_argument("--container", default="nipype/nipype:latest", help="Docker image name (default: nipype/nipype:latest)")
    parser.add_argument("--script_in_container", default="nipype_tsnr_workflow.py", help="Path to NiPype workflow script inside container")
    args = parser.parse_args()

    run_nipype_tsnr(args.infile, args.outdir, args.container, args.script_in_container)

if __name__ == "__main__":
    main()


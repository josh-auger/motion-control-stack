#!/usr/bin/env python3
"""
Title: nipype_tsnr_workflow.py

Description:
Run NiPype TSNR workflow on a 4D fMRI dataset.
For use with "run_nipype_tsnr.py" script that calls nipype docker container. This python script must be present in the
same directory as the input 4D image data file.

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Date of creation: October 17, 2025
"""

import os
import argparse
import nibabel as nib
import numpy as np
import pandas as pd
from nipype.algorithms.confounds import TSNR
from nipype import Node, Workflow

def run_tsnr(infile, outdir):
    """Run NiPype TSNR workflow on a 4D fMRI dataset."""

    os.makedirs(outdir, exist_ok=True)
    print(f"Input file : {infile}")
    print(f"Output directory : {outdir}")

    # Nipype TSNR Node
    tsnr_node = Node(TSNR(in_file=infile), name="tsnr")
    tsnr_node.base_dir = outdir  # Working directory

    # Define workflow
    wf = Workflow(name="tsnr_workflow", base_dir=outdir)
    wf.add_nodes([tsnr_node])
    print("Running NiPype TSNR workflow...")
    wf.run()

    # Locate generated files
    mean_file = os.path.join(outdir, "tsnr_workflow", "tsnr", "tsnr_mean.nii.gz")
    std_file = os.path.join(outdir, "tsnr_workflow", "tsnr", "tsnr_stddev.nii.gz")
    tsnr_file = os.path.join(outdir, "tsnr_workflow", "tsnr", "tsnr.nii.gz")

    # Verify outputs exist
    for f in [mean_file, std_file, tsnr_file]:
        if not os.path.exists(f):
            print(f"[Warning] Expected output not found: {f}")

    # Load data for summary CSV
    print("Computing summary statistics...")
    tsnr_img = nib.load(tsnr_file)
    tsnr_data = tsnr_img.get_fdata()

    summary = {
        "mean_tsnr": float(np.mean(tsnr_data)),
        "median_tsnr": float(np.median(tsnr_data)),
        "std_tsnr": float(np.std(tsnr_data)),
        "max_tsnr": float(np.max(tsnr_data)),
        "min_tsnr": float(np.min(tsnr_data)),
    }

    summary_csv = os.path.join(outdir, "tsnr_summary.csv")
    pd.DataFrame([summary]).to_csv(summary_csv, index=False)

    print("TSNR workflow complete.")
    print(f"Outputs saved to : {outdir}")
    print(f"Summary CSV : {summary_csv}")

def main():
    parser = argparse.ArgumentParser(description="Run NiPype TSNR workflow (container-side).")
    parser.add_argument("--infile", required=True, help="Path to input 4D .nii or .nii.gz file")
    parser.add_argument("--outdir", required=True, help="Path to output directory")
    args = parser.parse_args()

    run_tsnr(args.infile, args.outdir)

if __name__ == "__main__":
    main()

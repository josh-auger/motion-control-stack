#!/usr/bin/env python3

import os
import argparse
import SimpleITK as sitk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns


# -------------------------------
# HELPER FUNCTIONS
# -------------------------------
def parse_transform_str(s):
    """Convert comma-separated string of transform params to a float list."""
    s_clean = str(s).strip("[]")
    return [float(x) for x in s_clean.split(" ")]

def weighted_cov(X, w):
    """
    Compute weighted covariance matrix.
    X : (n_samples, n_features)
    w : weights (n_samples,)
    """
    X_centered = X - np.average(X, axis=0, weights=w)
    cov = np.dot((X_centered * w[:, None]).T, X_centered) / w.sum()
    return cov

def versor_to_euler_deg(params):
    """
    Convert VersorRigid3D parameters to Euler3D angles (degrees) + translation.
    Input
        params = [vx, vy, vz, tx, ty, tz]
    Output
        [rx_deg, ry_deg, rz_deg, tx, ty, tz]
    """
    vx, vy, vz, tx, ty, tz = params
    versor = sitk.VersorRigid3DTransform()
    versor.SetParameters((vx, vy, vz, tx, ty, tz))

    R = np.array(versor.GetMatrix()).reshape(3, 3)

    euler = sitk.Euler3DTransform()
    euler.SetMatrix(R.flatten())

    rx = np.degrees(euler.GetAngleX())
    ry = np.degrees(euler.GetAngleY())
    rz = np.degrees(euler.GetAngleZ())
    return [rx, ry, rz, tx, ty, tz]

def convert_column_to_euler(df, colname):
    """Parse Versor parameters column and convert to Euler degrees."""
    param_list = df[colname].apply(parse_transform_str).to_list()
    param_array = np.array(param_list)
    euler = np.array([versor_to_euler_deg(p) for p in param_array])
    return euler

def build_versor_transform(params, fixed_params):
    """Build a VersorRigid3DTransform from an array of Versor parameters [vx, vy, vz, tx, ty, tz] and center of
    rotation [cx, cy, cz].
    """
    t = sitk.VersorRigid3DTransform()
    t.SetCenter(fixed_params)
    t.SetParameters(params)
    return t

def plot_parameter_distributions(data, labels, title):
    """Create 6 subplot histogram figure."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()
    for i in range(6):
        vals = data[:, i]
        mean = np.mean(vals)
        std = np.std(vals)
        axes[i].hist(vals, bins=50)
        axes[i].set_title(f"{labels[i]}\nmean={mean:.3e}, std={std:.3e}")
        axes[i].set_xlabel("Value")
        axes[i].set_ylabel("Count")

    fig.suptitle(title, fontsize=16)
    plt.tight_layout()
    plt.show()

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


# -------------------------------
# MAIN
# -------------------------------
def main(master_csv_path):
    initial_params_col = "VersorTransformParams"
    fixed_params_col = "FixedParams"
    initial_mi_col = "MutualInfo_initial"

    stage1_params_col = "VersorTransformParams_stage1"
    stage1_mi_col = "MutualInfo_stage1"

    final_params_col = "VersorTransformParams_final"
    final_mi_col = "MutualInfo_final"

    param_labels = [
        "rot_x_deg",
        "rot_y_deg",
        "rot_z_deg",
        "trans_x_mm",
        "trans_y_mm",
        "trans_z_mm"
    ]

    # -------------------------------
    # LOAD DATA
    # -------------------------------
    df = pd.read_csv(master_csv_path)
    df = df.dropna(subset=[initial_params_col, fixed_params_col, final_params_col, initial_mi_col, final_mi_col])
    print(f"\nLoaded {len(df)} valid search runs.")
    initial_euler = convert_column_to_euler(df, initial_params_col)
    print("Converted initial transforms to Euler angles.")
    stage1_euler = convert_column_to_euler(df, stage1_params_col)
    print("Converted stage1 transforms to Euler angles.")
    final_euler = convert_column_to_euler(df, final_params_col)
    print("Converted final transforms to Euler angles.")
    initial_mi_values = df[initial_mi_col].to_numpy()
    stage1_mi_values = df[stage1_mi_col].to_numpy()
    final_mi_values = df[final_mi_col].to_numpy()

    # -------------------------------
    # LOAD RUNTIME DATA
    # -------------------------------
    n_evals1 = df["n_evals1"].to_numpy()
    init_time1 = df["total_init_time1"].to_numpy()
    parallel_time1 = df["total_parallel_time1"].to_numpy()
    serial_time1 = df["total_serial_time1"].to_numpy()
    total_reg_time1 = df["total_reg_time1"].to_numpy()

    n_evals2 = df["n_evals2"].to_numpy()
    init_time2 = df["total_init_time2"].to_numpy()
    parallel_time2 = df["total_parallel_time2"].to_numpy()
    serial_time2 = df["total_serial_time2"].to_numpy()
    total_reg_time2 = df["total_reg_time2"].to_numpy()
    total_reg_time_all = df["total_reg_time_all"].to_numpy()

    n_evals = n_evals1 + n_evals2
    case_index = np.arange(len(df))

    # Number of evaluations vs total runtime
    plt.figure(figsize=(8, 6))
    plt.scatter(n_evals, total_reg_time_all, s=10, alpha=0.5)
    plt.xlabel("Number of function evaluations")
    plt.ylabel("Total registration time (seconds)")
    plt.title("Registration runtime vs number of function evaluations")
    plt.xlim(100,550)
    plt.ylim(0.5,5.5)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Normalize runtimes by number of function evaluations
    init_time_per_eval = (init_time1 + init_time2) / n_evals
    parallel_time_per_eval = (parallel_time1 + parallel_time2) / n_evals
    serial_time_per_eval = (serial_time1 + serial_time2) / n_evals
    total_time_per_eval = total_reg_time_all / n_evals
    # Scatter plot of runtimes per evaluation
    plt.figure(figsize=(10, 6))
    plt.scatter(case_index, init_time_per_eval, s=10, alpha=0.5, label="Init time / eval")
    plt.scatter(case_index, parallel_time_per_eval, s=10, alpha=0.5, label="Parallel time / eval")
    plt.scatter(case_index, serial_time_per_eval, s=10, alpha=0.5, label="Serial time / eval")
    plt.scatter(case_index, total_time_per_eval, s=10, alpha=0.5, label="Total reg time / eval")
    plt.xlabel("Search case index")
    plt.ylabel("Runtime per evaluation (seconds)")
    plt.title("Registration runtime breakdown normalized by number of function evaluations")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Initial vs final parameters
    # norm = mpl.colors.Normalize(vmin=1.680, vmax=1.705)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()
    for i in range(6):
        sc = axes[i].scatter(
            initial_euler[:, i],
            final_euler[:, i],
            c=final_mi_values,
            cmap="viridis",
            # norm=norm,
            s=10,
            alpha=0.6
        )
        axes[i].set_xlabel(f"Initial {param_labels[i]}")
        axes[i].set_ylabel(f"Final {param_labels[i]}")
        axes[i].set_title(f"Initial vs Final {param_labels[i]}")
        axes[i].set_xlim(-15, 15)
        axes[i].set_ylim(-0.5, 0.5)
    # Create ONE shared colorbar for all subplots
    cbar = fig.colorbar(sc, ax=axes, location="right", shrink=0.9, pad=0.05, aspect=20)
    cbar.set_label("Final Mutual Information")
    fig.suptitle("Final Alignment Parameter vs Initial Starting Point", fontsize=16)
    # plt.tight_layout()
    plt.show()

    # -------------------------------
    # MUTUAL INFORMATION LANDSCAPE
    # -------------------------------
    print("Plotting mutual information landscape by parameter.")
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()
    for i in range(6):
        axes[i].scatter(
            initial_euler[:, i],
            initial_mi_values,
            s=10,
            alpha=0.4
        )
        axes[i].set_xlabel(param_labels[i])
        axes[i].set_ylabel("Mutual Information")
        axes[i].set_title(f"MI vs {param_labels[i]}")
        axes[i].set_xlim(-15, 15)
        axes[i].set_ylim(0.6, 1.90)
    fig.suptitle("Initial Mutual Information Landscape Projection", fontsize=16)
    plt.tight_layout()
    plt.show()

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()
    for i in range(6):
        axes[i].scatter(
            stage1_euler[:, i],
            stage1_mi_values,
            s=10,
            alpha=0.4
        )
        axes[i].set_xlabel(param_labels[i])
        axes[i].set_ylabel("Mutual Information")
        axes[i].set_title(f"MI vs {param_labels[i]}")
        axes[i].set_ylim(1.0, 1.90)
    fig.suptitle("Stage 1 Mutual Information Landscape Projection", fontsize=16)
    plt.tight_layout()
    plt.show()

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()
    for i in range(6):
        axes[i].scatter(
            final_euler[:, i],
            final_mi_values,
            s=10,
            alpha=0.4
        )
        axes[i].set_xlabel(param_labels[i])
        axes[i].set_ylabel("Mutual Information")
        axes[i].set_title(f"MI vs {param_labels[i]}")
        axes[i].set_xlim(-0.5, 0.5)
        axes[i].set_ylim(1.80, 1.90)
    fig.suptitle("Final Mutual Information Landscape Projection", fontsize=16)
    plt.tight_layout()
    plt.show()

    # -------------------------------
    # DISPLACEMENTS
    # -------------------------------
    print("Computing initial and final framewise displacements from the known solution.")
    success_threshold = 1.87   # MI value to define a successful registration convergence
    df["success"] = df[final_mi_col] > success_threshold

    displacements_initial = []
    displacements_stage1 = []
    displacements_final = []
    success_flags = []
    identity_transform = sitk.VersorRigid3DTransform()
    identity_transform.SetIdentity()

    for idx, row in df.iterrows():
        try:
            init_params = parse_transform_str(row[initial_params_col])
            stage1_params = parse_transform_str(row[stage1_params_col])
            final_params = parse_transform_str(row[final_params_col])
            fixed_params = parse_transform_str(row[fixed_params_col])
            start_transform = build_versor_transform(init_params, fixed_params)
            stage1_transform = build_versor_transform(stage1_params, fixed_params)
            end_transform = build_versor_transform(final_params, fixed_params)

            identity_transform.SetCenter(fixed_params)      # Identity must share same rotation center

            # Compute displacement (Tisdall et al. metric)
            composed_euler_initial = compose_transform_pair(start_transform, identity_transform)
            displacement_initial = calculate_displacement(composed_euler_initial, radius=50)
            displacements_initial.append(displacement_initial)

            composed_euler_stage1 = compose_transform_pair(stage1_transform, identity_transform)
            displacement_stage1 = calculate_displacement(composed_euler_stage1, radius=50)
            displacements_stage1.append(displacement_stage1)

            composed_euler_final = compose_transform_pair(end_transform, identity_transform)
            displacement_final = calculate_displacement(composed_euler_final, radius=50)
            displacements_final.append(displacement_final)

            success_flags.append(row["success"])
        except Exception as e:
            print(f"Skipping row {idx}: {e}")

    displacements_initial = np.asarray(displacements_initial)
    displacements_stage1 = np.asarray(displacements_stage1)
    displacements_final = np.asarray(displacements_final)
    success_flags = np.asarray(success_flags)
    print(f"Computed initial and final displacements for {len(displacements_initial)} runs.")

    # Initial displacement vs Initial Mutual Information
    plt.figure(figsize=(8, 6))
    plt.scatter(
        displacements_initial,
        initial_mi_values,
        s=12,
        alpha=0.5
    )
    plt.xlabel("Initial transform displacement from identity (mm)")
    plt.ylabel("Initial mutual information")
    plt.title("Initial mutual information as a function of starting distance from identity")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Stage 1 displacement vs Stage 1 Mutual Information
    plt.figure(figsize=(8, 6))
    plt.scatter(
        displacements_stage1,
        stage1_mi_values,
        s=12,
        alpha=0.5
    )
    plt.xlabel("Stage 1 transform displacement from identity (mm)")
    plt.ylabel("Stage 1 mutual information")
    plt.title("Stage 1 mutual information as a function of distance from identity")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Final displacement vs Final Mutual Information
    plt.figure(figsize=(8, 6))
    plt.scatter(
        displacements_final,
        final_mi_values,
        s=12,
        alpha=0.5
    )
    plt.xlabel("Final transform displacement from identity (mm)")
    plt.ylabel("Final mutual information")
    plt.title("Final mutual information and final distance from identity")
    # plt.xlim(0.0, 1.0)
    # plt.ylim(1.865, 1.88)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Initial displacement vs Final Mutual Information
    plt.figure(figsize=(8, 6))
    plt.scatter(
        displacements_initial,
        final_mi_values,
        s=12,
        alpha=0.5
    )
    plt.xlabel("Initial transform displacement from identity (mm)")
    plt.ylabel("Final mutual information")
    plt.title("Final mutual information after registration vs starting distance from identity")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Overlay success threshold
    success_mask = success_flags
    failure_mask = ~success_flags
    n_success = np.sum(success_mask)
    n_failure = np.sum(failure_mask)
    plt.figure(figsize=(8, 6))
    # Successful registrations
    plt.scatter(
        displacements_initial[success_mask],
        final_mi_values[success_mask],
        color="green",
        alpha=0.6,
        s=14,
        label=f"Success (n={n_success})"
    )
    # Failed registrations
    plt.scatter(
        displacements_initial[failure_mask],
        final_mi_values[failure_mask],
        color="red",
        alpha=0.6,
        s=14,
        label=f"Failure (n={n_failure})"
    )
    # Success threshold line
    plt.axhline(
        success_threshold,
        linestyle="--",
        color="black",
        label=f"Success Threshold (MI = {success_threshold})"
    )
    plt.xlabel("Initial transform displacement from identity (mm)")
    plt.ylabel("Final mutual information")
    plt.title("Final mutual information after registration vs starting distance from identity")
    plt.xlim(-1, 25)
    plt.ylim(1.865, 1.88)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # ------------------------------------------------
    # INITIAL vs FINAL DISPLACEMENT
    # ------------------------------------------------
    # norm = mpl.colors.Normalize(vmin=1.680, vmax=1.705)
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(
        displacements_initial,
        displacements_final,
        c=final_mi_values,
        cmap="viridis",
        # norm=norm,
        s=12,
        alpha=0.6
    )
    plt.xlabel("Initial transform displacement from identity (mm)")
    plt.ylabel("Final transform displacement from identity (mm)")
    plt.title("Registration outcome vs starting point (coded by final MI)")
    plt.xlim(-1, 25)
    plt.ylim(0.0, 1.0)
    plt.grid(True)
    cbar = plt.colorbar(scatter)
    cbar.set_label("Final Mutual Information")
    plt.tight_layout()
    plt.show()

    # ------------------------------------------------
    # ADD DERIVED COLUMNS TO MASTER DATAFRAME
    # ------------------------------------------------
    df["case_index"] = np.arange(len(df))
    df["initial_displacement"] = displacements_initial
    df["final_displacement"] = displacements_final
    # Add final transform parameters
    for i, label in enumerate(param_labels):
        df[f"final_{label}"] = final_euler[:, i]

    # ------------------------------------------------
    # REPORT INTERESTING CASES
    # ------------------------------------------------
    cases_of_interest = set()
    print("\nCases with lowest final MI value")
    print("------------------------------------------------")
    lowest_mi_cases = df.nsmallest(5, "MutualInfo_final")
    cases_of_interest.update(lowest_mi_cases.index)
    print(lowest_mi_cases[[
        "case_index",
        "VersorTransformParams",
        "MutualInfo_final",
        "initial_displacement",
        "final_displacement"
    ]])

    print("\nCases with largest final displacement")
    print("------------------------------------------------")
    largest_final_disp = df.nlargest(5, "final_displacement")
    cases_of_interest.update(largest_final_disp.index)
    print(largest_final_disp[[
        "case_index",
        "VersorTransformParams",
        "MutualInfo_final",
        "initial_displacement",
        "final_displacement"
    ]])

    print("\nCases with largest magnitude final parameter values")
    print("------------------------------------------------")
    for label in param_labels:
        col = f"final_{label}"
        top_cases = df.iloc[np.argsort(np.abs(df[col]))[::-1][:2]]
        cases_of_interest.update(top_cases.index)
        print(f"\nTop 2 |{label}| cases:")
        print(top_cases[[
            col,
            "case_index",
            "VersorTransformParams",
            "MutualInfo_final",
            "initial_displacement",
            "final_displacement"
        ]])

    print("\nCases that started close but failed (final MI < success threshold (3.0))")
    print("------------------------------------------------")
    failed_near_identity = df[df["MutualInfo_final"] < success_threshold]
    lowest_initial_disp_failures = failed_near_identity.nsmallest(5, "initial_displacement")
    cases_of_interest.update(lowest_initial_disp_failures.index)
    print(lowest_initial_disp_failures[[
        "case_index",
        "VersorTransformParams",
        "MutualInfo_final",
        "initial_displacement",
        "final_displacement",
    ]])

    cases_of_interest = sorted(list(cases_of_interest))
    print(f"\n\nNumber of interesting cases : {len(cases_of_interest)}")
    print(f"Case indices : {cases_of_interest}")

# -------------------------------
# ENTRY POINT
# -------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze random search transform experiment.")
    parser.add_argument("master_csv_path",help="Path to the master CSV table")
    args = parser.parse_args()

    main(args.master_csv_path)
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
    final_params_col = "VersorTransformParams_final"
    initial_mi_col = "MutualInfo_initial"
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
    final_euler = convert_column_to_euler(df, final_params_col)
    print("Converted final transforms to Euler angles.")
    initial_mi_values = df[initial_mi_col].to_numpy()
    final_mi_values = df[final_mi_col].to_numpy()

    # -------------------------------
    # LOAD RUNTIME DATA
    # -------------------------------
    n_evals = df["n_evals"].to_numpy()
    init_time = df["total_init_time"].to_numpy()
    parallel_time = df["total_parallel_time"].to_numpy()
    serial_time = df["total_serial_time"].to_numpy()
    total_reg_time = df["total_reg_time"].to_numpy()
    case_index = np.arange(len(df))

    # -------------------------------
    # FUNCTION EVALUATIONS PER CASE
    # -------------------------------
    # # Histogram of number of evaluations of each search case
    # plt.figure(figsize=(8, 6))
    # plt.hist(n_evals, bins=50)
    # plt.xlabel("Number of function evaluations")
    # plt.ylabel("Count of search cases")
    # plt.title("Distribution of optimizer evaluations")
    # plt.tight_layout()
    # plt.show()

    # -------------------------------
    # RUNTIME BREAKDOWN
    # -------------------------------
    # # Histogram of runtime distribution
    # plt.figure(figsize=(8, 6))
    # plt.hist(total_reg_time, bins=50)
    # plt.xlabel("Total registration time (seconds)")
    # plt.ylabel("Count of search cases")
    # plt.title("Distribution of registration runtimes")
    # plt.tight_layout()
    # plt.show()

    # Number of evaluations vs total runtime
    plt.figure(figsize=(8, 6))
    plt.scatter(n_evals, total_reg_time, s=10, alpha=0.5)
    plt.xlabel("Number of function evaluations")
    plt.ylabel("Total registration time (seconds)")
    plt.title("Registration runtime vs number of function evaluations")
    # plt.xlim(100, 550)
    # plt.ylim(0.5, 5.5)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # # Scatter plot of runtimes per case
    # plt.figure(figsize=(10, 6))
    # plt.scatter(case_index, init_time, s=10, alpha=0.5, label="Initialization time")
    # plt.scatter(case_index, parallel_time, s=10, alpha=0.5, label="Parallel compute time")
    # plt.scatter(case_index, serial_time, s=10, alpha=0.5, label="Serial compute time")
    # plt.scatter(case_index, total_reg_time, s=10, alpha=0.5, label="Total registration time")
    # plt.xlabel("Search case index")
    # plt.ylabel("Runtime (seconds)")
    # plt.title("Registration runtime breakdown for each search case")
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()
    # plt.show()

    # Normalize runtimes by number of function evaluations
    init_time_per_eval = init_time / n_evals
    parallel_time_per_eval = parallel_time / n_evals
    serial_time_per_eval = serial_time / n_evals
    total_time_per_eval = total_reg_time / n_evals
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

    # -------------------------------
    # PLOT PARAMETER DISTRIBUTIONS
    # -------------------------------
    # plot_parameter_distributions(
    #     initial_euler,
    #     param_labels,
    #     "Distribution of Initial Search Parameters"
    # )

    # plot_parameter_distributions(
    #     final_euler,
    #     param_labels,
    #     "Distribution of Final Alignment Parameters"
    # )

    # # Difference in final - initial parameters
    # delta_params = final_euler - initial_euler
    # plot_parameter_distributions(
    #     delta_params,
    #     param_labels,
    #     "Optimizer Parameter Movement (Final - Initial)"
    # )

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
        # axes[i].set_xlim(-15, 15)
        # axes[i].set_ylim(-0.5, 0.5)
    # Create ONE shared colorbar for all subplots
    cbar = fig.colorbar(sc, ax=axes, location="right", shrink=0.9, pad=0.05, aspect=20)
    cbar.set_label("Final Mutual Information")
    fig.suptitle("Final Alignment Parameter vs Initial Starting Point", fontsize=16)
    # plt.tight_layout()
    plt.show()

    # -------------------------------
    # COVARIANCE MATRICES
    # -------------------------------
    # cov_matrix_initial = np.cov(initial_euler, rowvar=False)
    # print("\nUnweighted covariance matrix (6x6) of initial params :")
    # print(cov_matrix_initial)
    #
    # cov_matrix_initial_weighted = weighted_cov(initial_euler, final_mi_values)
    # print("\nWeighted covariance matrix (6x6) by MutualInfo_final of initial params:")
    # print(cov_matrix_initial_weighted)
    #
    # cov_matrix_final = np.cov(final_euler, rowvar=False)
    # print("\nUnweighted covariance matrix (6x6) of final params :")
    # print(cov_matrix_final)
    #
    # cov_matrix_final_weighted = weighted_cov(final_euler, final_mi_values)
    # print("\nWeighted covariance matrix (6x6) by MutualInfo_final of final params:")
    # print(cov_matrix_final_weighted)


    # -------------------------------
    # HEATMAP VISUALIZATION
    # -------------------------------
    # plt.figure(figsize=(8,6))
    # sns.heatmap(
    #     cov_matrix_initial,
    #     annot=True,
    #     fmt=".2e",
    #     xticklabels=param_labels,
    #     yticklabels=param_labels,
    #     cmap="viridis"
    # )
    # plt.title("Unweighted Covariance Matrix of Initial Parameters")
    # plt.tight_layout()
    # plt.show()
    #
    # plt.figure(figsize=(8,6))
    # sns.heatmap(
    #     cov_matrix_initial_weighted,
    #     annot=True,
    #     fmt=".2e",
    #     xticklabels=param_labels,
    #     yticklabels=param_labels,
    #     cmap="viridis"
    # )
    # plt.title("Weighted Covariance Matrix of Initial Parameters (by MutualInfo_final)")
    # plt.tight_layout()
    # plt.show()

    # plt.figure(figsize=(8, 6))
    # sns.heatmap(
    #     cov_matrix_final,
    #     annot=True,
    #     fmt=".2e",
    #     xticklabels=param_labels,
    #     yticklabels=param_labels,
    #     cmap="viridis"
    # )
    # plt.title("Unweighted Covariance Matrix of Final Alignment Parameters")
    # plt.tight_layout()
    # plt.show()
    #
    # plt.figure(figsize=(8, 6))
    # sns.heatmap(
    #     cov_matrix_final_weighted,
    #     annot=True,
    #     fmt=".2e",
    #     xticklabels=param_labels,
    #     yticklabels=param_labels,
    #     cmap="viridis"
    # )
    # plt.title("Weighted Covariance Matrix of Final Alignment Parameters (by MutualInfo_final)")
    # plt.tight_layout()
    # plt.show()

    # -------------------------------
    # PAIRWISE PARAMETER CORRELATION
    # -------------------------------
    # print("Compiling pairwise parameter correlation of final parameters.")
    # corr_df = pd.DataFrame(final_euler, columns=param_labels)
    # corr_df["MutualInfo"] = final_mi_values
    # # Set up the color map for MutualInfo
    # cmap = mpl.cm.viridis
    # norm = mpl.colors.Normalize(vmin=corr_df["MutualInfo"].min(), vmax=corr_df["MutualInfo"].max())
    # sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])  # for the colorbar
    # # Create PairGrid
    # g = sns.PairGrid(corr_df, vars=param_labels, corner=True)
    # # Diagonal: histograms
    # g.map_diag(sns.histplot, kde=False, color="skyblue")
    # # Lower triangle: scatter plots colored by MutualInfo
    # def scatter_with_color(x, y, **kwargs):
    #     plt.scatter(x, y, c=corr_df["MutualInfo"], cmap=cmap, norm=norm, alpha=0.7, edgecolor=None)
    # g.map_lower(scatter_with_color)
    # # Add colorbar to the figure (not axes)
    # g.fig.colorbar(sm, label="MutualInfo_final", orientation="vertical")
    # # Titles and layout
    # g.fig.suptitle("Pairwise Correlation of Final Parameters\n(color-coded by final MI value)",fontsize=16,y=1.0)
    # plt.tight_layout()
    # plt.show()

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
        # axes[i].set_xlim(-15, 15)
        # axes[i].set_ylim(0.6, 1.8)
    fig.suptitle("Initial Mutual Information Landscape Projection", fontsize=16)
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
        # axes[i].set_xlim(-0.5, 0.5)
        # axes[i].set_ylim(1.678, 1.705)
    fig.suptitle("Final Mutual Information Landscape Projection", fontsize=16)
    plt.tight_layout()
    plt.show()

    # -------------------------------
    # DISPLACEMENTS
    # -------------------------------
    print("Computing initial and final framewise displacements from the known solution.")
    success_threshold = 1.70   # MI value to define a successful registration convergence
    df["success"] = df[final_mi_col] > success_threshold

    displacements_initial = []
    displacements_final = []
    success_flags = []
    identity_transform = sitk.VersorRigid3DTransform()
    identity_transform.SetIdentity()

    for idx, row in df.iterrows():
        try:
            init_params = parse_transform_str(row[initial_params_col])
            final_params = parse_transform_str(row[final_params_col])
            fixed_params = parse_transform_str(row[fixed_params_col])
            start_transform = build_versor_transform(init_params, fixed_params)
            end_transform = build_versor_transform(final_params, fixed_params)

            identity_transform.SetCenter(fixed_params)      # Identity must share same rotation center

            # Compute displacement (Tisdall et al. metric)
            composed_euler_initial = compose_transform_pair(start_transform, identity_transform)
            displacement_initial = calculate_displacement(composed_euler_initial, radius=50)
            displacements_initial.append(displacement_initial)

            composed_euler_final = compose_transform_pair(end_transform, identity_transform)
            displacement_final = calculate_displacement(composed_euler_final, radius=50)
            displacements_final.append(displacement_final)

            success_flags.append(row["success"])
        except Exception as e:
            print(f"Skipping row {idx}: {e}")

    displacements_initial = np.asarray(displacements_initial)
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
    # plt.ylim(1.678, 1.705)
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
    # plt.xlim(-1, 25)
    # plt.ylim(1.678, 1.705)
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
    # plt.xlim(-1, 25)
    # plt.ylim(0.0, 1.0)
    plt.grid(True)
    cbar = plt.colorbar(scatter)
    cbar.set_label("Final Mutual Information")
    plt.tight_layout()
    plt.show()

    # ------------------------------------------------
    # EVALUATION OF REGISTRATION AGAINST RUNTIMES
    # ------------------------------------------------
    # # Initial displacement vs total registration runtime
    # plt.figure(figsize=(8, 6))
    # scatter = plt.scatter(
    #     displacements_initial,
    #     total_reg_time,
    #     c=final_mi_values,
    #     cmap="viridis",
    #     s=12,
    #     alpha=0.6
    # )
    # plt.xlabel("Initial transform displacement from identity (mm)")
    # plt.ylabel("Total registration runtime (sec)")
    # plt.title("Time to register vs starting point (coded by final MI)")
    # plt.grid(True)
    # cbar = plt.colorbar(scatter)
    # cbar.set_label("Final Mutual Information")
    # plt.tight_layout()
    # plt.show()

    # # Initial displacement vs number of function evaluations
    # plt.figure(figsize=(8, 6))
    # scatter = plt.scatter(
    #     displacements_initial,
    #     n_evals,
    #     c=final_mi_values,
    #     cmap="viridis",
    #     s=12,
    #     alpha=0.6
    # )
    # plt.xlabel("Initial transform displacement from identity (mm)")
    # plt.ylabel("Number of function evaluations")
    # plt.title("Optimizer evaluations vs starting point (coded by final MI)")
    # plt.grid(True)
    # cbar = plt.colorbar(scatter)
    # cbar.set_label("Final Mutual Information")
    # plt.tight_layout()
    # plt.show()

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

    # # Mark interesting cases on final vs initial displacement plot
    # highlight_initial = displacements_initial[cases_of_interest]
    # highlight_final = displacements_final[cases_of_interest]
    # plt.figure(figsize=(8, 6))
    # scatter = plt.scatter(
    #     displacements_initial,
    #     displacements_final,
    #     c=final_mi_values,
    #     cmap="viridis",
    #     s=12,
    #     alpha=0.6
    # )
    # # Overlay highlighted cases
    # plt.scatter(
    #     highlight_initial,
    #     highlight_final,
    #     facecolors="none",
    #     edgecolors="red",
    #     s=120,
    #     linewidths=1.5,
    #     zorder=3,
    #     label="Extreme failed cases"
    # )
    # # Overlay highlighted cases
    # plt.scatter(
    #     displacements_initial[lowest_initial_disp_failures.index],
    #     displacements_final[lowest_initial_disp_failures.index],
    #     facecolors="none",
    #     edgecolors="teal",
    #     s=120,
    #     linewidths=1.5,
    #     zorder=3,
    #     label="Almost there cases"
    # )
    # plt.xlabel("Initial transform displacement from identity (mm)")
    # plt.ylabel("Final transform displacement from identity (mm)")
    # plt.title("Registration outcome vs starting point (coded by final MI)")
    # plt.grid(True)
    # cbar = plt.colorbar(scatter)
    # cbar.set_label("Final Mutual Information")
    # plt.legend()
    # plt.tight_layout()
    # plt.show()

    # ------------------------------------------------
    # SAVE TRANSFORMS OF INTERESTING CASES
    # ------------------------------------------------
    # output_dir = "interesting_case_transforms"
    # os.makedirs(output_dir, exist_ok=True)
    # for idx in cases_of_interest:
    #     try:
    #         row = df.loc[idx]
    #         # Parse parameters
    #         initial_versor_params = parse_transform_str(row["VersorTransformParams"])
    #         final_versor_params = parse_transform_str(row["VersorTransformParams_final"])
    #         fixed_params = parse_transform_str(row["FixedParams"])
    #         # Build initial transform
    #         t_init = sitk.VersorRigid3DTransform()
    #         t_init.SetParameters(initial_versor_params)
    #         t_init.SetFixedParameters(fixed_params)
    #         # Build final transform
    #         t_final = sitk.VersorRigid3DTransform()
    #         t_final.SetParameters(final_versor_params)
    #         t_final.SetFixedParameters(fixed_params)
    #         # Output filenames
    #         filename_init = os.path.join(output_dir, f"inputTransform_{idx:04d}.tfm")
    #         filename_final = os.path.join(output_dir, f"alignTransform_{idx:04d}.tfm")
    #         # Save transform
    #         sitk.WriteTransform(t_init, filename_init)
    #         sitk.WriteTransform(t_final, filename_final)
    #         print(f"Saved: {filename_init}")
    #         print(f"Saved: {filename_final}")
    #     except Exception as e:
    #         print(f"Skipping case {idx}: {e}")

# -------------------------------
# ENTRY POINT
# -------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze random search transform experiment.")
    parser.add_argument("master_csv_path",help="Path to the master CSV table")
    args = parser.parse_args()

    main(args.master_csv_path)
# Title: analyze_regtrace.py

# Description:
# Post-processing of the trace of evaluated transform parameters and resulting mutual information metrics during a
# single registration optimization.
#

# Created on: July 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.interpolate import griddata

def load_csv_to_numpy(file_path):

    df = pd.read_csv(file_path)

    # Split 'EvaluatedTransformParams' into separate columns
    transform_params = df['EvaluatedTransformParams'].str.split(' ', expand=True).astype(float)
    transform_params.columns = [f'TransformParam_{i + 1}' for i in range(transform_params.shape[1])]

    # Split 'FixedParams' into separate columns
    fixed_params = df['FixedParams'].str.split(' ', expand=True).astype(float)
    fixed_params.columns = [f'FixedParam_{i + 1}' for i in range(fixed_params.shape[1])]

    combined_df = pd.concat([transform_params, fixed_params, df[['MutualInfo']]], axis=1)
    data_array = combined_df.to_numpy()

    return data_array

def find_max_value(column_vector):
    max_index = np.argmax(column_vector)
    max_value = column_vector[max_index]

    print(f"Maximum MI value : {max_value}")
    print(f"Maximum MI iteration : {max_index}")

    return max_index, max_value


def plot_iteration_series(data_array):
    transform_params = data_array[:, :6]
    mutual_info = data_array[:, -1]

    # Create a figure with 2 subplots (2 rows, 1 column)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Plot transform parameters
    x = np.arange(len(transform_params))
    for i in range(3):
        ax1.plot(x, transform_params[:, i], label=f'Rotation {["X", "Y", "Z"][i]}')
    for i in range(3, 6):
        ax1.plot(x, transform_params[:, i], linestyle='--', label=f'Translation {["X", "Y", "Z"][i - 3]}')

    ax1.set_title('Transform Parameters')
    ax1.set_ylabel('rads/mm')
    ax1.legend()

    # Plot mutual information
    ax2.plot(x, mutual_info, marker='o', markersize=4, color='r', label='Mutual Information')

    # Mark the maximum MI value
    max_mi_index, max_mi_value = find_max_value(mutual_info)
    ax2.plot(max_mi_index, max_mi_value, 'o', markersize=8, markerfacecolor='none', markeredgewidth=2, color='blue', label=f'Max MI [{max_mi_index}]={max_mi_value:.4f}')

    ax2.set_title('Mutual Information')
    ax2.set_xlabel('Optimization Iteration')
    ax2.set_ylabel('MI Value')
    ax2.legend()

    plt.tight_layout()
    plt.show()


def plot_3d_scatter(ax, x, y, z, x_label, y_label, z_label, color_map):
    scatter = ax.scatter(x, y, z, c=z, cmap=color_map)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    return scatter

def plot_transform_params_and_mi_3d(data_array):
    transform_params = data_array[:, :6]
    mutual_info = data_array[:, -1]
    color_map = 'viridis'

    fig = plt.figure(figsize=(18, 6))
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.05])

    # X Translation vs X Rotation
    ax1 = fig.add_subplot(gs[0], projection='3d')
    scatter1 = plot_3d_scatter(ax1, transform_params[:, 0], transform_params[:, 3], mutual_info,
                               'Rotation X (rads)', 'Translation X (mm)', 'Mutual Information', color_map)
    ax1.set_title('Rotation X vs Translation X')

    # Y Translation vs Y Rotation
    ax2 = fig.add_subplot(gs[1], projection='3d')
    scatter2 = plot_3d_scatter(ax2, transform_params[:, 1], transform_params[:, 4], mutual_info,
                               'Rotation Y (rads)', 'Translation Y (mm)', 'Mutual Information', color_map)
    ax2.set_title('Rotation Y vs Translation Y')

    # Z Translation vs Z Rotation
    ax3 = fig.add_subplot(gs[2], projection='3d')
    scatter3 = plot_3d_scatter(ax3, transform_params[:, 2], transform_params[:, 5], mutual_info,
                               'Rotation Z (rads)', 'Translation Z (mm)', 'Mutual Information', color_map)
    ax3.set_title('Rotation Z vs Translation Z')

    # # X Translation vs Y Translation
    # ax3 = fig.add_subplot(gs[2], projection='3d')
    # scatter3 = plot_3d_scatter(ax3, transform_params[:, 3], transform_params[:, 4], mutual_info,
    #                            'Translation X (mm)', 'Translation Y (mm)', 'Mutual Information', color_map)
    # ax3.set_title('Translation X vs Translation Y')

    # Add a single color bar for all subplots
    cbar_ax = fig.add_subplot(gs[3])
    cbar = fig.colorbar(scatter3, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Mutual Information')

    plt.tight_layout()
    plt.show()


def plot_2d_optimization_path(data_array, param1_idx=3, param2_idx=4, param1_label='X', param2_label='Y'):
    x_values = data_array[:, param1_idx]
    y_values = data_array[:, param2_idx]
    mutual_info = data_array[:, -1]

    # Create a DataFrame to handle duplicate (x, y) pairs
    data_df = pd.DataFrame({'x': x_values, 'y': y_values, 'z': mutual_info})
    # Group by (x, y) and keep only the max z for each pair
    max_z_df = data_df.groupby(['x', 'y'], as_index=False).max()

    # Extract filtered values
    x_values = max_z_df['x'].values
    y_values = max_z_df['y'].values
    mutual_info = max_z_df['z'].values

    # Create a grid for contour plot
    xi = np.linspace(x_values.min(), x_values.max(), 100)
    yi = np.linspace(y_values.min(), y_values.max(), 100)
    xi, yi = np.meshgrid(xi, yi)
    zi = griddata((x_values, y_values), mutual_info, (xi, yi), method='nearest')
    zi = np.clip(zi, mutual_info.min(), mutual_info.max())

    fig, ax = plt.subplots(figsize=(8, 8))

    # Contour plot for mutual information
    contour = ax.contourf(xi, yi, zi, levels=20, cmap='viridis')
    plt.colorbar(contour, label='Mutual Information')

    # Plot the optimization path
    ax.plot(x_values, y_values, marker='o', linestyle='-', color='red', markersize=3, linewidth=0)

    # Mark the first optimization iteration with a blue X marker
    ax.plot(x_values[0], y_values[0], marker='o', markersize=10, color='blue', markerfacecolor='none',
            markeredgewidth=2, label='Initial eval')

    # Mark the maximum mutual information value with an open blue circle
    max_mi_index = np.argmax(mutual_info)
    max_mi_value = mutual_info[max_mi_index]
    ax.plot(x_values[max_mi_index], y_values[max_mi_index], marker='x', markersize=10, color='blue',
            markerfacecolor='none', markeredgewidth=2, label=f'Max MI: {max_mi_value:.4f}')

    ax.set_title(f'Optimization Path for {param1_label} vs {param2_label}')
    ax.set_xlabel(f'{param1_label} (mm/rads)')
    ax.set_ylabel(f'{param2_label} (mm/rads)')
    ax.legend()

    plt.show()

def plot_mi_vs_transform_params(data_array):
    transform_params = data_array[:, :6]
    mutual_info = data_array[:, -1]

    param_labels = [
        'Rotation X (rads)',
        'Rotation Y (rads)',
        'Rotation Z (rads)',
        'Translation X (mm)',
        'Translation Y (mm)',
        'Translation Z (mm)'
    ]

    fig, axes = plt.subplots(2, 3, figsize=(24, 8), sharey=False)
    axes = axes.flatten()

    for i, ax in enumerate(axes):
        ax.scatter(transform_params[:, i], mutual_info, color='blue', s=15, alpha=0.6, edgecolors='none')

        # Compute and plot center line
        center_x = (np.min(transform_params[:, i]) + np.max(transform_params[:, i])) / 2
        ax.axvline(center_x, color='red', linestyle='--', linewidth=1.5)

        ax.set_xlabel(param_labels[i])
        ax.set_ylabel('Mutual Information')
        ax.set_title(f'Mutual Information vs {param_labels[i]}')
        ax.grid(True)

    plt.tight_layout(pad=3.0)
    plt.show()



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 plot_regtrace.py <csv_file_path>")
        sys.exit(1)

    csv_file_path = sys.argv[1]
    data_array = load_csv_to_numpy(csv_file_path)
    # print(data_array)

    # plot_iteration_series(data_array)
    plot_mi_vs_transform_params(data_array)
    # plot_transform_params_and_mi_3d(data_array)
    # plot_2d_optimization_path(data_array, param1_idx=3, param2_idx=4, param1_label='Trans X', param2_label='Trans Y')
    # plot_2d_optimization_path(data_array, param1_idx=3, param2_idx=5, param1_label='Trans X', param2_label='Trans Z')
    # plot_2d_optimization_path(data_array, param1_idx=4, param2_idx=5, param1_label='Trans Y', param2_label='Trans Z')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D



def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so spheres appear as spheres, cubes as cubes, etc.
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max(x_range, y_range, z_range)
    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)
    ax.set_xlim3d([x_middle - max_range/2, x_middle + max_range/2])
    ax.set_ylim3d([y_middle - max_range/2, y_middle + max_range/2])
    ax.set_zlim3d([z_middle - max_range/2, z_middle + max_range/2])

def visualize_coordinate_frame(direction_matrix, origin=[0, 0, 0]):
    """
    Visualizes a 3D coordinate frame and compares it with the global device frame.

    Parameters:
    - direction_matrix: 3x3 numpy array where each column is a direction vector (dir_x, dir_y, dir_z)
    - origin: 3-element list or array specifying the origin (default: [0, 0, 0])
    """
    origin = np.array(origin)
    direction_matrix = np.array(direction_matrix)
    assert direction_matrix.shape == (3, 3), "direction_matrix must be 3x3"

    # Normalize input direction matrix (column-wise)
    direction_matrix = direction_matrix / np.linalg.norm(direction_matrix, axis=0)

    print(f"Image origin : {origin}")
    print(f"Image direction_matrix : {direction_matrix}")

    # Global/device coordinate frame (identity)
    identity_matrix = np.eye(3)
    global_origin = np.array([-126.511627, -126.511627, -88.5])

    print(f"Device origin : {global_origin}")
    print(f"Device direction_matrix : {identity_matrix}")

    # Create figure and 3D axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot image frame origin marker
    ax.scatter(*origin, color='black', s=50)

    # Plot global/device frame origin marker
    ax.scatter(*global_origin, color='red', s=50, alpha=0.4)

    # Axis labels and colors
    axis_labels = ['X', 'Y', 'Z']
    colors = ['r', 'g', 'b']

    # --- Plot global/device frame (dotted, faded, different origin) ---
    vec_length = 15  # or adjust to fit your scene scale
    for i in range(3):
        vec = identity_matrix[:, i] * vec_length
        ax.quiver(*global_origin, *vec,
                  color=colors[i], linestyle='--', linewidth=1.5, alpha=0.4)

    # --- Plot input/image frame (solid, bold) ---
    for i in range(3):
        vec = direction_matrix[:, i] * vec_length
        ax.quiver(*origin, *vec,
                  color=colors[i], linewidth=1.5)

    legend_elements = [
        Line2D([0], [0], color='k', marker='o', linestyle='', label=f'Image Origin\n({origin[0]:.2f}, {origin[1]:.2f}, {origin[2]:.2f})'),
        Line2D([0], [0], color='r', marker='o', linestyle='', label=f'Device Origin\n({global_origin[0]:.2f}, {global_origin[1]:.2f}, {global_origin[2]:.2f})'),
        Line2D([0], [0], color='grey', linestyle='--', label='Device Frame'),
        Line2D([0], [0], color='grey', linestyle='-', label='Image Frame'),
    ]
    ax.legend(handles=legend_elements,loc='upper left', fontsize='small')

    # Labels and legend
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Image Coordinate Frame vs Device Frame')
    set_axes_equal(ax)
    plt.tight_layout()
    plt.show()


# Example usage
if __name__ == "__main__":
    # Define direction vectors
    dir_x = np.array([2.931520,-0.516906,-0.0])
    dir_y = np.array([0.516906,2.931520,0.0])
    dir_z = np.array([-0.0,-0.0,3.0])

    direction_matrix = np.column_stack((dir_x, dir_y, dir_z))
    origin = [-126.558140, -102.621125, -88.5]  # Default origin if not specified

    # Visualize coordinate frame
    visualize_coordinate_frame(direction_matrix, origin)

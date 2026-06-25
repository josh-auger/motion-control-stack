"""
Title: convert_transform_for_scanner.py

Description:
    Helper script to solve a hand-eye calibration problem, a brute force solution to find the unknown coordinate frame
    conversion between sms-mi-reg VVR results frame (which may differ from the specified image frame in NRRD header)
    and the global device frame of the scanner.

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Created on: 13 August 2025

"""

import os
import argparse
import numpy as np
from scipy.spatial.transform import Rotation as R
import SimpleITK as sitk


def read_transforms_from_dir(dir_path):
    """
    Reads all .tfm transform files from a directory and converts them into
    4x4 numpy transformation matrices.
    """
    if not os.path.isdir(dir_path):
        raise FileNotFoundError(f"Directory not found: {dir_path}")

    tfm_files = [f for f in os.listdir(dir_path) if f.lower().endswith('.tfm')]
    tfm_files.sort()  # Ensure consistent ordering

    mats = []
    count = 0
    for tfm_file in tfm_files:
        count += 1
        print(f"Reading transform {count:03d}")
        tfm_path = os.path.join(dir_path, tfm_file)
        transform = sitk.ReadTransform(tfm_path)
        params = transform.GetParameters()  # rotation (and possibly translation)
        print(f"Transformation parameters : {params}")
        fixed_params = transform.GetFixedParameters()  # center of rotation
        print(f"Fixed parameters : {fixed_params}")
        rot_mat = np.array(transform.GetMatrix()).reshape(3, 3)
        trans_vec = np.array(transform.GetTranslation())

        # Construct 4x4 matrix
        mat_4x4 = np.eye(4)
        mat_4x4[:3, :3] = rot_mat
        mat_4x4[:3, 3] = trans_vec

        mats.append(mat_4x4)

    return mats


def read_and_negate_transforms_from_dir(dir_path):
    """
    Reads all .tfm transform files from a directory and converts them into
    4x4 numpy transformation matrices.
    """
    if not os.path.isdir(dir_path):
        raise FileNotFoundError(f"Directory not found: {dir_path}")

    tfm_files = [f for f in os.listdir(dir_path) if f.lower().endswith('.tfm')]
    tfm_files.sort()  # Ensure consistent ordering

    mats = []
    count = 0
    for tfm_file in tfm_files:
        count += 1
        print(f"Reading transform {count:03d}")
        tfm_path = os.path.join(dir_path, tfm_file)
        transform = sitk.ReadTransform(tfm_path)
        # print(f"Read transform : ")
        # print(transform)
        params = transform.GetParameters()  # rotation (and possibly translation)
        print(f"Read parameters : {params}")
        fixed_params = transform.GetFixedParameters()  # center of rotation
        print(f"Read fixed parameters : {fixed_params}")

        params = [params[x] * -1 for x in range(6)]
        print(f"Negated parameters : {params}")
        transform.SetParameters(params)
        # print(f"Negate transform : ")
        # print(transform)

        rot_mat = np.array(transform.GetMatrix()).reshape(3, 3)
        trans_vec = np.array(transform.GetTranslation())

        # Construct 4x4 matrix
        mat_4x4 = np.eye(4)
        mat_4x4[:3, :3] = rot_mat
        mat_4x4[:3, 3] = trans_vec

        mats.append(mat_4x4)

    return mats


def skew(v):
    """Return the skew-symmetric matrix of a vector."""
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])


def solve_hand_eye_tsai_lenz(A_mats, B_mats):
    """
    Solve AX = XB for X using the Tsai–Lenz separable solution method. Decompose the equation into 3x3 rotation matrix
    and 3x1 translation vector.
    R_A*R_X = R_X*R_B
    R_A*t_X + t_A = R_X*t_B + t_X
    """
    assert len(A_mats) == len(B_mats), "A and B must have the same number of motions"
    n = len(A_mats)

    # --- Step 1: Solve rotation RX ---
    L = []
    r = []
    for i in range(n):
        R_a = A_mats[i][:3, :3]
        R_b = B_mats[i][:3, :3]
        alpha = R.from_matrix(R_a).as_rotvec()
        beta = R.from_matrix(R_b).as_rotvec()

        # skew(alpha + beta) * p = beta - alpha
        L.append(skew(alpha + beta))
        r.append((beta - alpha).reshape(3, 1))

    L = np.vstack(L)
    r = np.vstack(r)

    # Solve L * p = r for p (rotation vector of RX)
    p, _, _, _ = np.linalg.lstsq(L, r, rcond=None)
    theta = np.linalg.norm(p)
    R_x = R.from_rotvec((p.flatten())).as_matrix()

    # --- Step 2: Solve translation tX using linear least squares---
    C = []
    d = []
    for i in range(n):
        R_a = A_mats[i][:3, :3]
        t_a = A_mats[i][:3, 3]
        R_b = B_mats[i][:3, :3]
        t_b = B_mats[i][:3, 3]

        C.append(R_a - np.eye(3))
        d.append((R_x @ t_b - t_a).reshape(3, 1))  # ensure column vector

    C = np.vstack(C)    # size 3*n x 3
    d = np.vstack(d)    # size 3*n x 1

    t_x, _, _, _ = np.linalg.lstsq(C, d, rcond=None)
    # t_x = np.zeros(3)   # Ideally, there is no translation b/w coordinate frames (e.g. VVR is accurate)

    # --- Step 3: Assemble X ---
    X = np.eye(4)
    X[:3, :3] = R_x
    X[:3, 3] = t_x.flatten()

    return X


def rigid_transform_3D(A, B):
    """
    Kabsch algorithm: find the rigid transform components (rotation, R  and translation, t) that maps the
    point coordinates A to the point coordinates B, such that R @ Ai + t = Bi is satisfied using least-squares fit.
    A : source point set, estimated locations (N x 3)
    B : target point set, known locations (N x 3)
    """
    assert A.shape == B.shape
    N = A.shape[0]

    # compute centroids, the mean of each point set
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)

    # subtract centroids to get point clouds centered at the origin (first solve for rotation without translation)
    AA = A - centroid_A
    BB = B - centroid_B

    # Covariance: how points in A correlate to points in B. Each entry of H tells how much each axis of A aligns with each axis of B
    H = AA.T @ BB   # (3xN)(Nx3) = (3x3)

    # # Singular value decomposition. Find rotation that minimizes squared error b/w R @ A and B
    # U, S, Vt = np.linalg.svd(H)
    # R = Vt.T @ U.T
    #
    # # reflection check to make sure determinant = +1, otherwise flip last column of V
    # if np.linalg.det(R) < 0:
    #     Vt[-1, :] *= -1
    #     R = Vt.T @ U.T

    # Eigen-decomposition of H.T*H
    eigvals, eigvecs = np.linalg.eigh(H.T @ H)

    # Take eigenvector corresponding to largest eigenvalue
    v = eigvecs[:, np.argmax(eigvals)]
    u = H @ v
    u /= np.linalg.norm(u)

    # Construct rotation using outer product (closest rotation matrix)
    R = np.outer(u, v)

    # Ensure R is a proper rotation (determinant = +1)
    if np.linalg.det(R) < 0:
        u = -u
        R = np.outer(u, v)

    # Translation: difference in centroids after applying the found rotation, R, to points of B
    t = centroid_B - R @ centroid_A
    return R, t


def solve_hand_eye_point_method(A_mats, B_mats, p=np.array([1, 1, 1])):
    """
    Solve AX = XB using point correspondences generated by a single reference point p.
    """
    assert len(A_mats) == len(B_mats), "A and B must have the same number of motions"
    n = len(A_mats)

    p_A_list = []
    p_B_list = []

    # For each transform pair, A and B, apply each transform to the point, p
    # for A, B in zip(A_mats, B_mats):
    #     # Transform reference point by A (from VVR) and by B (expected moco response) to get to the moved target location
    #     p_A = A[:3, :3] @ p + A[:3,3]
    #     p_B = B[:3, :3] @ p + B[:3,3]
    #     p_A_list.append(p_A)
    #     p_B_list.append(p_B)

    # For each point, p, apply transform A and B to each point
    print('A:', A_mats[0])
    print('B:', B_mats[0])

    for point in p:
        p_A = A_mats[0][:3, :3] @ point + A_mats[0][:3, 3]
        p_B = B_mats[0][:3, :3] @ point + B_mats[0][:3, 3]
        p_A_list.append(p_A)
        p_B_list.append(p_B)

    p_A_array = np.stack(p_A_list, axis=0)  # reference point mapped via A to moved target, shape (n, 3)
    p_B_array = np.stack(p_B_list, axis=0)  # reference point mapped via B to moved target, shape (n, 3)
    print(f"p_A_array: {p_A_array.shape}")
    print(f"{p_A_array}")
    print(f"p_B_array: {p_B_array.shape}")
    print(f"{p_B_array}")

    # Solve for rigid transform components that maps point coordinates B to the point coordinates A
    R_x, t_x = rigid_transform_3D(p_A_array, p_B_array)
    # t_x = np.zeros(3)   # Ideally, there is no translation b/w coordinate frames (e.g. VVR is accurate)

    # Assemble transformation matrix, X, from solved rotation and translation
    X = np.eye(4)
    X[:3, :3] = R_x
    X[:3, 3] = t_x
    return X


# ===== Example usage =====
def main():
    # parser = argparse.ArgumentParser(description="Solve AX = XB using Tsai–Lenz hand–eye calibration")
    # parser.add_argument("--dirA", required=True, help="Path to directory of A Versor3DRigid transforms (*.tfm files)")
    # parser.add_argument("--dirB", required=True, help="Path to directory of B Versor3DRigid transforms (*.tfm files)")
    # args = parser.parse_args()
    # A_mats = read_transforms_from_dir(args.dirA)
    # B_mats = read_and_negate_transforms_from_dir(args.dirB)

    # dir_A = transforms found by sms-mi-reg VVR to align the reference position to the moved target
    # dir_B = moco feedback Euler3D transforms sent to the scanner to move the FoV
    dir_A = "/home/jauger/Radiology_Research/Scan_data/20250813_Phantom_scan/fromPythonFireServer/hand_eye_calibration_analysis/smsmireg_retroVVR_alignTransforms_Euler3D_20250813"
    dir_B = "/home/jauger/Radiology_Research/Scan_data/20250813_Phantom_scan/fromPythonFireServer/hand_eye_calibration_analysis/mocoTransforms_20250813"
    A_mats = read_transforms_from_dir(dir_A)
    B_mats = read_and_negate_transforms_from_dir(dir_B) # negate the input moco feedback to get the expected moco response to re-center the FoV on the phantom

    print(f"Read {len(A_mats)} transforms for source A (to be mapped to B)")
    print(f"Read {len(B_mats)} transforms for target B (expected results)")

    # Set numpy print formatting
    np.set_printoptions(precision=3, suppress=True)

    # solve hand-eye calibration using Tsai-Lenz approach
    X = solve_hand_eye_tsai_lenz(A_mats, B_mats)
    print("Estimated X (Tsai–Lenz) :\n", X)

    # solve hand-eye calibration using point cloud approach
    points = [
        np.array([-126.5116, -126.5116, -88.5]),
        np.array([-126.5116, 126.5116, -88.5]),
        np.array([126.5116, -126.5116, -88.5]),
        np.array([126.5116, 126.5116, -88.5]),
        np.array([-126.5116, -126.5116, 88.5]),
        np.array([-126.5116, 126.5116, 88.5]),
        np.array([126.5116, -126.5116, 88.5]),
        np.array([126.5116, 126.5116, 88.5])
    ]
    X_point = solve_hand_eye_point_method(A_mats, B_mats, points)
    print("Estimated X (point-based method):\n", X_point)


if __name__ == "__main__":
    main()

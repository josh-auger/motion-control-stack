"""
Title: convert_transform_for_scanner.py

Description:
    Convert the final alignment transform from sms-mi-reg into the device coordinate system for sending back to the
    scanner.

Author: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital
Created on: 3 April 2025

"""

import SimpleITK as sitk
import numpy as np
import argparse
import ctypes
import time
import os


def recenter_transform(transform, new_center):
    # Calculate the equivalent transform centered at the origin
    A0 = np.asarray(transform.GetMatrix()).reshape(3, 3)
    c0 = np.asarray(transform.GetCenter())
    t0 = np.asarray(transform.GetTranslation())
    combined_mat = A0
    c1 = new_center
    combined_center = c1
    center_difference = c1 - c0
    combined_translation = t0 + np.dot(A0, center_difference) - center_difference

    new_transform = construct_transform_object(combined_mat, combined_center, combined_translation, "Affine")
    return new_transform


def convert_patientLPS_to_deviceFrame(M_LPS):
    """Apply the coordinate system transformation M' = T * M to convert a transform from LPS+ frame to the global device
    coordinate frame.
    image LPS+ X == device X => patient left
    image LPS+ Y => patient posterior (back); device Y => patient anterior (front)
    image LPS+ Z => patient superior (head); device Z => patient inferior (feet)
    """
    # Construct 4x4 transformation matrix of rotation and translation
    T_LPStoDev = np.array([
        [1, 0, 0, 0],
        [0, -1, 0, 0],
        [0, 0, -1, 0],
        [0, 0, 0, 1]
    ], dtype=float)
    print(f"\nCoordinate transformation matrix T_LPStoDev (LPS+ to device frame) :")
    print(f"{T_LPStoDev}")

    # Apply transformation
    M_dev = T_LPStoDev @ M_LPS
    return M_dev


def convert_imageFrame_to_deviceFrame(transform, R_img):
    """Apply a coordinate system chance of basis, M' = T_inverse * M, to the alignment transform in order to convert
    from the image coordinate frame to the patient LPS+ frame and then to the global device coordinate frame."""
    # Image direction matrix (columns are direction vectors) from image header direction vectors
    print(f"\n********** Debugging coordinate frame conversion **********")
    print(f"\nImage header direction cosine matrix :")
    print(R_img)

    # Extract affine components of the registration alignment transform
    A_img = np.asarray(transform.GetMatrix()).reshape(3, 3)
    t_img = np.asarray(transform.GetTranslation())
    center = np.asarray(transform.GetCenter())
    print(f"VVR transform rotation matrix (in image frame?) :")
    print(A_img)
    print(f"VVR transform translation (in image frame?) :")
    print(t_img)
    # Construct a 4x4 matrix of the alignment transform components
    M_img = np.eye(4)
    M_img[:3, :3] = A_img
    M_img[:3, 3] = t_img
    print(f"VVR transform matrix M (in image frame?) :")
    print(M_img)

    # # Construct a 4x4 matrix to define the coordinate frame transformation, T, from image frame to patient LPS+ (back to identity directions)
    # T = np.eye(4)
    # T[:3, :3] = R_img
    # # T[:3, 3] = offset_vector # where offset_vector is the delta b/w image center and device center
    # T_inv = np.linalg.inv(T)
    # print(f"\nCoordinate transformation matrix T (image frame to LPS+) :")
    # print(T_inv)

    # Apply change of basis to alignment transform matrix, M, to convert to patient LPS+ frame
    # M_lps = T_inv @ M_img
    # print(f"Transform matrix M in LPS+ frame :")
    # print(M_lps)

    # JDA: Keep transform in sms-mi-reg output frame (it is not the expected image frame or LPS+, but is close to being in the device frame)
    M_lps = M_img.copy()

    # # Apply change of basis to convert from LPS+ frame to device coordinate frame
    # M_dev = convert_patientLPS_to_deviceFrame(M_lps)
    # A_dev = M_dev[:3, :3]
    # t_dev = M_dev[:3, 3]
    # print(f"\nTransformation matrix M_new in device frame :")
    # print(M_dev)

    # JDA: Brute force an x-axis parameters flip to put sms-mi-reg results in the global device frame for moco feedback
    euler_transform_lps = construct_transform_object(M_lps[:3, :3], center, M_lps[:3, 3], "Euler3D")
    params = np.asarray(euler_transform_lps.GetParameters())
    params_flipX = params.copy()
    params_flipX[0] *= -1
    params_flipX[3] *= -1
    euler_transform_dev = sitk.Euler3DTransform()
    euler_transform_dev.SetParameters(params_flipX)
    euler_transform_dev.SetCenter(center)
    A_dev = np.asarray(euler_transform_dev.GetMatrix()).reshape(3, 3)
    t_dev = np.asarray(euler_transform_dev.GetTranslation())

    fmt = lambda arr: np.array(arr, dtype=float).round(6).tolist()
    print(f"\nVVR transform euler params : {fmt(params)}")
    print(f"Device frame (X-flipped) transform euler params : {fmt(params_flipX)}")
    print(f"Rotation matrix in device frame :")
    print(A_dev)
    print(f"Translation in device frame :")
    print(t_dev)

    # Construct new affine transform from the components in the device frame, centered at [0,0,0]
    new_transform_affine = construct_transform_object(A_dev, center, t_dev, "Affine")
    new_transform_euler = construct_transform_object(A_dev, center, t_dev, "Euler3D")
    print(f"\nConverted euler transform (in device frame) :")
    print(new_transform_euler)

    return new_transform_affine


def construct_transform_object(rotation_matrix, center, translation, transformtype):
    if transformtype == "Affine":
        new_transform = sitk.AffineTransform(3)
        new_transform.SetCenter(center)
        new_transform.SetMatrix(rotation_matrix.ravel())
        new_transform.SetTranslation(translation)

    elif transformtype == "VersorRigid3D":
        new_transform = sitk.VersorRigid3DTransform()
        new_transform.SetCenter(center)
        new_transform.SetMatrix(rotation_matrix.ravel())
        new_transform.SetTranslation(translation)

    elif transformtype == "Euler3D":
        new_transform = sitk.Euler3DTransform()
        new_transform.SetCenter(center)
        new_transform.SetMatrix(rotation_matrix.ravel())
        new_transform.SetTranslation(translation)

    else:
        raise ValueError("Invalid transform type!")

    return new_transform


def convert_transform_for_moco(transform, direction_matrix):
    # Convert to an Affine transform
    A0 = np.asarray(transform.GetMatrix()).reshape(3, 3)
    c0 = np.asarray(transform.GetCenter())
    t0 = np.asarray(transform.GetTranslation())
    new_transform = construct_transform_object(A0, c0, t0, "Affine")
    print(f"Input affine transform (in image frame?) :")
    print(new_transform)
    print(f"Input transform as Euler3D (in image frame?) :")
    print(construct_transform_object(A0, c0, t0, "Euler3D"))

    # Calculate equivalent transform centered at the origin [0,0,0], which is the device isocenter
    new_transform = recenter_transform(new_transform, np.zeros(3))
    print(f"\nRe-center affine transform about origin :")
    print(new_transform)

    # Convert transform written w.r.t. reference image coordinate frame to be w.r.t. device coordinate frame
    new_transform = convert_imageFrame_to_deviceFrame(new_transform, direction_matrix)
    print(f"\nConverted affine transform (in device frame) :")
    print(new_transform)

    return new_transform


class MyMocoFeedbackData(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ('mB_mocoInfoValid', ctypes.c_bool),             # 1 byte
        ('mD_vecTrans', ctypes.c_double * 3),            # 24 bytes
        ('mD_rotMat', ctypes.c_double * 9),              # 72 bytes
        ('mI_acquisitionTime', ctypes.c_int * 1),        # 4 bytes
        ('mI_frameNo', ctypes.c_int * 1)                 # 4 bytes
    ]  # Total = 105 bytes


def package_transform_as_datastruct(affine_transform, acquisition_time, frame_number):
    """ Packages MOCO affine transform into a MyMocoFeedbackData struct. """
    feedback = MyMocoFeedbackData()
    feedback.mB_mocoInfoValid = True

    # Extract rotation matrix and translation vector
    matrix = np.asarray(affine_transform.GetMatrix()).reshape(3, 3)
    translation = np.asarray(affine_transform.GetTranslation())
    # Assign translation vector (x, y, z)
    for i in range(3):
        feedback.mD_vecTrans[i] = translation[i]
    # Assign rotation matrix (row-major order)
    flat_matrix = matrix.flatten()
    for i in range(9):
        feedback.mD_rotMat[i] = flat_matrix[i]

    # Assign acquisition time --> pull from first slice image[0] header, head[0].acquisition_time_stamp
    feedback.mI_acquisitionTime[0] = acquisition_time
    # Assign frame number --> set here (stays paired to this feedback) and then assigned to subsequent image acquisitions by Siemens MOCO
    feedback.mI_frameNo[0] = frame_number

    return feedback


def format_moco_struct(feedback_struct):
    """Converts the struct into a nicely formatted string for saving as a .txt file for reference.
    Float values are truncated to 6 decimal places."""
    def fmt(val):
        if isinstance(val, float):
            return f"{val:.6f}"
        return val

    lines = []
    lines.append(f"mB_mocoInfoValid: {feedback_struct.mB_mocoInfoValid}")
    lines.append(f"mD_vecTrans: {[fmt(v) for v in feedback_struct.mD_vecTrans]}")
    lines.append(f"mD_rotMat: {[fmt(v) for v in feedback_struct.mD_rotMat]}")
    lines.append(f"mI_acquisitionTime: {feedback_struct.mI_acquisitionTime[0]}")
    lines.append(f"mI_frameNo: {feedback_struct.mI_frameNo[0]}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Convert a 3D rigid transform to an affine transform.")
    parser.add_argument("--transform", required=True, help="Input VersorRigid3D transform from sms-mi-reg.")
    parser.add_argument("--refimage", required=True, help="Input reference image used to define image frame directions.")
    parser.add_argument("--outputtransform", required=True, help="Output filename for the converted transform.")
    args = parser.parse_args()

    # Read transform and image
    read_transform = sitk.ReadTransform(args.transform)
    refimage = sitk.ReadImage(args.refimage)
    directions = np.array(refimage.GetDirection()).reshape(3, 3)
    dir_x = directions[:, 0] / np.linalg.norm(directions[:, 0])
    dir_y = directions[:, 1] / np.linalg.norm(directions[:, 1])
    dir_z = directions[:, 2] / np.linalg.norm(directions[:, 2])
    print("\nImage Direction Matrix (image frame axes relative to device frame):")
    print(f"  X direction  : {dir_x}")
    print(f"  Y direction : {dir_y}")
    print(f"  Z direction : {dir_z}\n")

    # Convert to Device Coordinate System
    direction_matrix = np.column_stack((dir_x, dir_y, dir_z))
    moco_transform = convert_transform_for_moco(read_transform, direction_matrix)

    # Write converted transform
    sitk.WriteTransform(moco_transform, args.outputtransform)

    # Derive struct file path using os.path
    base, _ = os.path.splitext(args.outputtransform)
    output_struct_path = base + ".txt"

    # Package and write moco struct
    moco_struct = package_transform_as_datastruct(moco_transform, int(time.time()), 5)
    with open(output_struct_path, 'w') as f:
        f.write(format_moco_struct(moco_struct))

if __name__ == "__main__":
    main()
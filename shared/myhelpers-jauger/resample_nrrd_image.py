
import argparse
import SimpleITK as sitk
import os


def resample_nrrd_volume(input_nhdr_path, upsample_factor=0.9, output_nhdr_path=None):
    """
    Upsample a 3D NRRD volume using SimpleITK and write out a detached NRRD (.nhdr + .raw).

    Parameters
    ----------
    input_nhdr_path : str
        Path to the input .nhdr file.
    upsample_factor : float
        Factor to scale voxel spacing (e.g., 0.9 reduces spacing by 10% → higher resolution).
    output_nhdr_path : str, optional
        Path to output .nhdr file. If None, appends '_upsampled' to input filename.

    Returns
    -------
    output_nhdr_path : str
        Path to the written .nhdr file.
    """

    # Resolve output path
    if output_nhdr_path is None:
        base, ext = os.path.splitext(input_nhdr_path)
        output_nhdr_path = base + "_upsampled.nhdr"

    # Read image
    image = sitk.ReadImage(input_nhdr_path)

    # Ensure 3D
    if image.GetDimension() != 3:
        raise ValueError("Input image is not 3D.")

    # Original properties
    original_size = image.GetSize()
    original_spacing = image.GetSpacing()
    origin = image.GetOrigin()
    direction = image.GetDirection()

    # Compute new spacing and size
    new_spacing = [s * upsample_factor for s in original_spacing]
    new_size = [int(round(osz * ospc / nspc)) for osz, ospc, nspc in zip(original_size, original_spacing, new_spacing)]

    # Resample
    resampled = sitk.Resample(
        image,
        new_size,
        sitk.Transform(),
        sitk.sitkLinear,
        origin,
        new_spacing,
        direction,
        0,  # default pixel value
        image.GetPixelID()
    )

    # Write as detached NRRD (NHDR header file + RAW image data file)
    writer = sitk.ImageFileWriter()
    writer.SetFileName(output_nhdr_path)
    writer.SetUseCompression(False)

    # Force detached header/raw format
    writer.SetImageIO("NrrdImageIO")
    writer.Execute(resampled)
    return output_nhdr_path


def main():
    parser = argparse.ArgumentParser(description="Downsample an ITK image.")
    parser.add_argument("--inputnrrd", required=True, help="Path to the input image file.")
    parser.add_argument("--outputnrrd", required=False, help="Filename for the output image.")
    parser.add_argument("--factor", type=float, default=2.0, help="Downsampling factor applied to all axes (default: 2.0).")
    args = parser.parse_args()

    resampled_image_path = resample_nrrd_volume(args.inputnrrd, args.factor)
    print(f"Downsampled image saved as: {resampled_image_path}")


if __name__ == "__main__":
    main()
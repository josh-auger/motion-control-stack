


import SimpleITK as sitk
import pydicom
import argparse

def read_nrrd_metadata(nrrd_file_path):
    # Read the NRRD file using SimpleITK
    image = sitk.ReadImage(nrrd_file_path)

    # Get the metadata dictionary
    metadata_keys = image.GetMetaDataKeys()

    # Print out the metadata header information
    print("NRRD Metadata Header Information:")
    for key in metadata_keys:
        value = image.GetMetaData(key)
        print(f"{key}: {value}")

def read_dicom_metadata(dicom_file_path):
    # Read the DICOM file using pydicom
    dataset = pydicom.dcmread(dicom_file_path)

    # Print out the metadata header information, excluding large/binary data elements
    print("DICOM Metadata Header Information:")
    for elem in dataset:
        # Check if the element is text-based and not too large
        if elem.VR != 'OB' and elem.VR != 'OW' and len(str(elem.value)) < 1024:
            print(f"{elem.tag} {elem.name}: {elem.value}")
        else:
            print(f"{elem.tag} {elem.name}: [Non-text or large data element]")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Read and print metadata header information from NRRD or DICOM files.")
    parser.add_argument('--nrrdfile', type=str, help="Path to the NRRD file")
    parser.add_argument('--dicomfile', type=str, help="Path to the DICOM file")

    # Parse the command line arguments
    args = parser.parse_args()

    # Check if an NRRD file was provided
    if args.nrrdfile:
        read_nrrd_metadata(args.nrrdfile)

    # Check if a DICOM file was provided
    if args.dicomfile:
        read_dicom_metadata(args.dicomfile)


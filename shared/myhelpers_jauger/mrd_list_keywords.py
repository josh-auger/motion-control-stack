
# Title: mrd_list_keywords.py

# Description:
# List the keyword hierarchy and number of individual images in an MRD file

# Created on: March 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

import h5py
import ismrmrd

def print_mrd_contents(group, indent=0):
    for key, item in group.items():
        if isinstance(item, h5py.Group):
            print("  " * indent + f"Group: {key}")
            print_mrd_contents(item, indent + 1)
        else:
            try:
                # Attempt to access the value attribute
                value = item.value
                print("  " * indent + f"Dataset: {key}: {value}")
            except AttributeError:
                # If value attribute does not exist, print the dataset name only
                print("  " * indent + f"Dataset: {key}")

def print_mrd_keywords(file_path):
    try:
        with h5py.File(file_path, 'r') as f:
            print("List of keys (groups) in the MRD file:")
            print_mrd_contents(f)
    except IOError:
        print("Error: Could not open file or file does not exist.")
    except Exception as e:
        print("Error:", e)



if __name__ == "__main__":
    # file_path = input("Enter the path to the MRD file (.h5): ")
    file_path = "/home/jauger/Radiology_Research/SLIMM_data/20231010_restingstate/KelvinMRD_saved_data/func-bold_task-rest960_run-01slimmon_2024-03-14-183816_29.h5"
    print_mrd_keywords(file_path)

    data = ismrmrd.Dataset(file_path, '/dataset', True)

    num_images = data.number_of_images("images_0")
    print("\nNumber of images: ", num_images)

    img400_header = data.read_image("images_0",400).meta
    img400_attr = data.read_image("images_0", 400).attribute_string
    img400_data = data.read_image("images_0", 400).data

    data.close()


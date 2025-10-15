#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "numpy",
#     "simpleitk",
#     "pydicom",
# ]
# ///

import SimpleITK as sitk
import argparse
import sys
import time
import os
import shutil
import numpy as np
import pydicom

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert a nifti or 3d-tiff to dicom given a donor dicom image",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="increase verbosity")
    parser.add_argument("-s", "--seriesdescription")
    parser.add_argument("-n", "--seriesnumber")
    parser.add_argument("nifti", help="nifti or 3d-tiff image")
    parser.add_argument("donor", help="dicom donor image")
    parser.add_argument("dicomdir", help="dicom output directory")
    return parser.parse_args(argv)


def writeSlices(series_tag_values,
                new_img: sitk.Image,
                out_dir,
                depth: int,
                writer: sitk.ImageFileWriter):
    image_slice = new_img[:, :, depth]

    # Tags shared by the series
    for tag, value in series_tag_values:
        try:
            image_slice.SetMetaData(tag, value)
        except Exception:
            print(f"tag: {tag}, value: {value}.")

    # Slice-specific tags
    image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))  # Creation Date
    image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S"))  # Creation Time
    image_slice.SetMetaData(
        "0020|0032",
        "\\".join(map(str, new_img.TransformIndexToPhysicalPoint((0, 0, depth)))),
    )  # Image Position
    image_slice.SetMetaData("0020|0013", str(depth))  # Instance Number
    image_slice.SetMetaData("0018|0088", str(new_img.GetSpacing()[2])) # Slice spacing, this is apparantly not used and not written by the writer...
    
    # Write slice
    writer.SetFileName(os.path.join(out_dir, str(depth).rjust(6, "0") + ".dcm"))
    writer.Execute(image_slice)


def prepare_image(nifti_input, is_tiff):
    nii_img = sitk.ReadImage(nifti_input)
    
    if not is_tiff:
        print("Converting the nifti to 16bit")
        arr = sitk.GetArrayFromImage(nii_img)
        max_val = np.amax(arr)
        img_int16 = arr * (np.iinfo(np.int16).max / max_val)
        img_int16b = img_int16.astype(np.int16)
        new_img = sitk.GetImageFromArray(img_int16b)
        new_img.CopyInformation(nii_img)
        
        new_img = sitk.DICOMOrient(new_img, "LPS")
    else:
        new_img = nii_img

    return new_img


def copy_tags(reader, new_img, seriesdesc, seriesnumber):
    tags_to_copy = [
        "0002|0002", "0010|0010", "0010|0020", "0010|0030", "0010|0040",
        "0020|000D", "0020|000d", "0020|0010", "0008|0016", "0008|0020",
        "0008|0022", "0008|0023", "0008|0030", "0008|0032", "0008|0033",
        "0008|0050", "0008|0060", "0008|0080", "0018|0050",
    ]

    modification_time = time.strftime("%H%M%S")
    modification_date = time.strftime("%Y%m%d")
    direction = new_img.GetDirection()


    
    series_tag_values_a = [
        (
            k,
            reader.GetMetaData(k)
            .encode('utf-8', 'surrogateescape')  # convert str -> bytes (preserving invalids)
            .decode('latin-1')                   # decode bytes as latin-1
        )
        for k in tags_to_copy if reader.HasMetaDataKey(k)
    ]
    series_tag_values_b = [
        ("0008|0031", modification_time),
        ("0008|0021", modification_date),
        ("0008|0008", "DERIVED\\SECONDARY"),
        ("0020|000e",
         "1.2.826.0.1.3680043.2.1125." + modification_date + ".1" + modification_time),
        ("0020|0037", "\\".join(map(str, (
            direction[0], direction[3], direction[6],
            direction[1], direction[4], direction[7],
        )))),
        ("0008|103e", seriesdesc),
        ("0020|0011", str(seriesnumber)),
    ]
    return series_tag_values_a + series_tag_values_b


def check_inputs(nifti_input, donor_dcm):
    """Validate that the input files exist."""
    if not os.path.exists(donor_dcm):
        return False, f"{donor_dcm} does not exist"
    if not os.path.exists(nifti_input):
        return False, f"{nifti_input} does not exist"
    return True, ""



def detect_input_type(nifti_input):
    """Return whether the input is a TIFF or NIfTI."""
    _, img_ext = os.path.splitext(nifti_input)
    is_tiff = (img_ext == ".tiff")
    return is_tiff


def read_dicom_metadata(donor_dcm, verbose=False):
    """Read DICOM metadata from donor image."""
    reader = sitk.ImageFileReader()
    reader.SetFileName(donor_dcm)
    reader.LoadPrivateTagsOn()
    reader.ReadImageInformation()

    if verbose:
        for k in reader.GetMetaDataKeys():
            v = reader.GetMetaData(k)
            try:
                print(f"({k}) = \"{v}\"")
            except Exception:
                print("An exception occurred")

    return reader


def prepare_output_dir(dcm_output):
    """Clean and recreate the output directory."""
    if os.path.exists(dcm_output):
        shutil.rmtree(dcm_output)
    os.makedirs(dcm_output, exist_ok=True)


def write_slices(new_img, series_tag_values, dcm_output, writer: sitk.ImageFileWriter):
    """Write all slices of the new image as DICOMs."""
    for i in range(new_img.GetDepth()):
        writeSlices(series_tag_values, new_img, dcm_output, i, writer)


def get_dicom_files_with_slice_location(directory):
    """
    Recursively searches a directory for DICOM files that contain
    the SliceLocation (0020,1041) tag.

    Parameters:
        directory (str): Path to the root directory containing DICOM files.

    Returns:
        list[str]: Alphabetically sorted list of file paths.
    """
    dicom_files_with_slice_location = []

    for root, _, files in os.walk(directory):
        for filename in files:
            filepath = os.path.join(root, filename)

            try:
                ds = pydicom.dcmread(filepath, stop_before_pixels=True, force=True)
                if 'SliceLocation' in ds:
                    dicom_files_with_slice_location.append(filepath)
            except Exception:
                # Skip files that are not valid DICOM
                continue

    return sorted(dicom_files_with_slice_location)


def get_all_dicom_files(directory):
    """
    Recursively searches a directory for all valid DICOM files.

    Parameters:
        directory (str): Path to the root directory containing DICOM files.

    Returns:
        list[str]: Alphabetically sorted list of all valid DICOM file paths.
    """
    all_dicom_files = []

    for root, _, files in os.walk(directory):
        for filename in files:
            filepath = os.path.join(root, filename)

            try:
                # Try reading minimal metadata to confirm it's a DICOM file
                pydicom.dcmread(filepath, stop_before_pixels=True, force=True)
                all_dicom_files.append(filepath)
            except Exception:
                # Skip files that are not valid DICOM
                continue

    return sorted(all_dicom_files)


def copy_location_thickness(from_list, to_list):
    print("Copying slice location, thickness and spacing between slices.")
    for from_file, to_file in zip(from_list, to_list):
        from_dcm = pydicom.dcmread(from_file, stop_before_pixels=True, force=True)
        to_file_data = pydicom.dcmread(to_file, stop_before_pixels=False, force=True)
        to_file_data.SliceLocation = from_dcm.SliceLocation
        to_file_data.SliceThickness = from_dcm.SliceThickness
        to_file_data.SpacingBetweenSlices = from_dcm.SpacingBetweenSlices
        to_file_data.save_as(to_file)


def convert_and_write2(nifti_input, donor_dcm, dcm_output, verbose=False,
                      seriesdesc="ITKsimple - KUL_NIS", seriesnumber=""):
    donor_dcm_dir = os.path.dirname(donor_dcm)
    donor_dcm_list = get_dicom_files_with_slice_location(donor_dcm_dir)
    donor_dcm = donor_dcm_list[0]
    print(donor_dcm)
    valid, msg = check_inputs(nifti_input, donor_dcm)
    if not valid:
        print(msg)
        return 1
    
    is_tiff = detect_input_type(nifti_input)
    print("Assuming input is a 3d-tiff" if is_tiff else "Assuming input is nifti")
    
    reader = read_dicom_metadata(donor_dcm, verbose)
    new_img = prepare_image(nifti_input, is_tiff)
    
    series_tag_values = copy_tags(reader, new_img, seriesdesc, seriesnumber)
    print("Incorporating the following dicom tags:")
    print(series_tag_values)

    prepare_output_dir(dcm_output)

    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()
    write_slices(new_img, series_tag_values, dcm_output, writer)

    # Copy SliceLocation and SliceThickness tags to be able to use the MPR in pacs enterprise.
    copy_location_thickness(donor_dcm_list, get_all_dicom_files(dcm_output))

    return 0        


def convert_and_write(nifti_input, donor_dcm, dcm_output, verbose=False,
                      seriesdesc="ITKsimple - KUL_NIS", seriesnumber=""):
    valid, msg = check_inputs(nifti_input, donor_dcm)
    if not valid:
        print(msg)
        return 1
    
    is_tiff = detect_input_type(nifti_input)
    print("Assuming input is a 3d-tiff" if is_tiff else "Assuming input is nifti")
    
    reader = read_dicom_metadata(donor_dcm, verbose)

    
    new_img = prepare_image(nifti_input, is_tiff)
    

    series_tag_values = copy_tags(reader, new_img, seriesdesc, seriesnumber)
    print("Incorporating the following dicom tags:")
    print(series_tag_values)

    prepare_output_dir(dcm_output)

    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()
    write_slices(new_img, series_tag_values, dcm_output, writer)

    return 0


def main(argv=None):
    args = parse_args(argv)
    return convert_and_write2(
        nifti_input=args.nifti,
        donor_dcm=args.donor,
        dcm_output=args.dicomdir,
        verbose=args.verbose,
        seriesdesc=args.seriesdescription or "IKTsimple - KUL_NIS",
        seriesnumber=args.seriesnumber or "",
    )


if __name__ == "__main__":
    sys.exit(main())

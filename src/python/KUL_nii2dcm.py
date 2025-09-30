#!/usr/bin/env -S uv run --script

import SimpleITK as sitk
import argparse
import sys
import time
import os
import shutil
import numpy as np


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


def writeSlices(series_tag_values, new_img, out_dir, i, writer):
    image_slice = new_img[:, :, i]

    # Tags shared by the series
    for tag, value in series_tag_values:
        image_slice.SetMetaData(tag, value)

    # Slice-specific tags
    image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))  # Creation Date
    image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S"))  # Creation Time
    image_slice.SetMetaData(
        "0020|0032",
        "\\".join(map(str, new_img.TransformIndexToPhysicalPoint((0, 0, i)))),
    )  # Image Position
    image_slice.SetMetaData("0020|0013", str(i))  # Instance Number

    # Write slice
    writer.SetFileName(os.path.join(out_dir, str(i).rjust(6, "0") + ".dcm"))
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
        "0008|0050", "0008|0060", "0008|0080",
    ]

    modification_time = time.strftime("%H%M%S")
    modification_date = time.strftime("%Y%m%d")
    direction = new_img.GetDirection()

    series_tag_values_a = [
        (k, reader.GetMetaData(k))
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
        ("0020|0011", seriesnumber),
    ]
    return series_tag_values_a + series_tag_values_b


def convert_and_write(nifti_input, donor_dcm, dcm_output, verbose=False,
                      seriesdesc="IKTsimple - KUL_NIS", seriesnumber=""):
    if not os.path.exists(donor_dcm):
        print(donor_dcm + " does not exist")
        return 1
    if not os.path.exists(nifti_input):
        print(nifti_input + " does not exist")
        return 1

    _, img_ext = os.path.splitext(nifti_input)
    is_tiff = (img_ext == ".tiff")
    print("Assuming input is a 3d-tiff" if is_tiff else "Assuming input is nifti")

    # Read donor dicom metadata
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

    # Prepare the image
    new_img = prepare_image(nifti_input, is_tiff)

    # Collect dicom tags
    series_tag_values = copy_tags(reader, new_img, seriesdesc, seriesnumber)

    print("Incorporating the following dicom tags:")
    print(series_tag_values)

    # Clean and recreate output dir
    if os.path.exists(dcm_output):
        shutil.rmtree(dcm_output)
    os.makedirs(dcm_output, exist_ok=True)

    # Writer
    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()

    # Write slices
    for i in range(new_img.GetDepth()):
        writeSlices(series_tag_values, new_img, dcm_output, i, writer)

    return 0


def main(argv=None):
    args = parse_args(argv)
    return convert_and_write(
        nifti_input=args.nifti,
        donor_dcm=args.donor,
        dcm_output=args.dicomdir,
        verbose=args.verbose,
        seriesdesc=args.seriesdescription or "IKTsimple - KUL_NIS",
        seriesnumber=args.seriesnumber or "",
    )


if __name__ == "__main__":
    sys.exit(main())


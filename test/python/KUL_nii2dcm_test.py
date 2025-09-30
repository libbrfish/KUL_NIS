#!/usr/bin/env -S uv run --script

import os
import tempfile
import shutil
import numpy as np
import uuid
import pytest
import time
import sys
from unittest.mock import MagicMock, patch

import SimpleITK as sitk
import KUL_nii2dcm as convert_to_dicom


@pytest.fixture
def fake_image():
    """Create a small fake 3D image using SimpleITK."""
    arr = np.zeros((5, 5, 3), dtype=np.int16)
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing([1.0, 1.0, 1.0])
    return img



def test_writeSlices_calls_writer(fake_image, tmp_path):
    tags = [("0008|0060", "MR")]
    mock_writer = MagicMock()

    convert_to_dicom.writeSlices(tags, fake_image, str(tmp_path), 0, mock_writer)

    mock_writer.SetFileName.assert_called_once()
    mock_writer.Execute.assert_called_once()


def test_prepare_image_nifti_conversion(monkeypatch):
    """Ensure NIfTI branch rescales to int16."""
    fake_arr = np.ones((5, 5, 5), dtype=np.float32)
    fake_img = sitk.GetImageFromArray(fake_arr)

    monkeypatch.setattr(sitk, "ReadImage", lambda x: fake_img)

    new_img = convert_to_dicom.prepare_image("dummy.nii.gz", is_tiff=False)

    assert isinstance(new_img, sitk.Image)
    assert new_img.GetPixelIDTypeAsString() == "16-bit signed integer"


def test_prepare_image_tiff(monkeypatch):
    """Ensure TIFF branch does not rescale."""
    fake_arr = np.ones((5, 5, 5), dtype=np.int16)
    fake_img = sitk.GetImageFromArray(fake_arr)

    monkeypatch.setattr(sitk, "ReadImage", lambda x: fake_img)

    new_img = convert_to_dicom.prepare_image("dummy.tiff", is_tiff=True)

    assert isinstance(new_img, sitk.Image)
    # Same dtype as input
    assert new_img.GetPixelIDTypeAsString() == "16-bit signed integer"


def test_copy_tags_returns_expected_keys(fake_image, monkeypatch):
    """Check if copy_tags collects metadata and adds new tags."""
    mock_reader = MagicMock()
    mock_reader.HasMetaDataKey.return_value = True
    mock_reader.GetMetaData.side_effect = lambda k: "value"

    tags = convert_to_dicom.copy_tags(mock_reader, fake_image, "desc", "123")

    keys = [k for k, _ in tags]
    assert "0008|103e" in keys  # Series Description
    assert "0020|0011" in keys  # Series Number
    assert "0020|000e" in keys  # Series Instance UID


def test_convert_and_write_missing_files(tmp_path):
    """Should return 1 when input files are missing."""
    nifti = tmp_path / "input.nii.gz"
    donor = tmp_path / "donor.dcm"
    outdir = tmp_path / "out"

    # Neither exists
    rc = convert_to_dicom.convert_and_write(str(nifti), str(donor), str(outdir))
    assert rc == 1


def test_convert_and_write_tiff(monkeypatch, tmp_path, fake_image):
    """Run through TIFF branch without actual IO."""
    donor = tmp_path / "donor.dcm"
    nifti = tmp_path / "input.tiff"
    outdir = tmp_path / "out"
    donor.write_text("dummy")
    nifti.write_text("dummy")

    # Mock SimpleITK pieces
    monkeypatch.setattr(sitk, "ReadImage", lambda x: fake_image)
    mock_reader = MagicMock()
    mock_reader.GetMetaDataKeys.return_value = []
    mock_reader.HasMetaDataKey.return_value = False
    mock_reader.GetMetaData.return_value = "val"
    monkeypatch.setattr(sitk, "ImageFileReader", lambda: mock_reader)
    monkeypatch.setattr(convert_to_dicom, "writeSlices", lambda *a, **k: None)

    rc = convert_to_dicom.convert_and_write(str(nifti), str(donor), str(outdir))
    assert rc == 0
    assert outdir.exists()


def test_main_runs(monkeypatch, tmp_path, fake_image):
    """Integration test for main()."""
    donor = tmp_path / "donor.dcm"
    nifti = tmp_path / "input.tiff"
    outdir = tmp_path / "out"
    donor.write_text("dummy")
    nifti.write_text("dummy")

    monkeypatch.setattr(sitk, "ReadImage", lambda x: fake_image)
    mock_reader = MagicMock()
    mock_reader.GetMetaDataKeys.return_value = []
    mock_reader.HasMetaDataKey.return_value = False
    mock_reader.GetMetaData.return_value = "val"
    monkeypatch.setattr(sitk, "ImageFileReader", lambda: mock_reader)
    monkeypatch.setattr(convert_to_dicom, "writeSlices", lambda *a, **k: None)

    rc = convert_to_dicom.main([str(nifti), str(donor), str(outdir)])
    assert rc == 0


@pytest.fixture
def fake_nifti_image():
    arr = np.arange(27).reshape((3, 3, 3)).astype(np.int16)
    img = sitk.GetImageFromArray(arr)
    return img


def create_fake_donor_dicom(tmp_path, patient_name="John^Doe"):
    import numpy as np
    arr = np.zeros((2, 2), dtype=np.int16)
    img = sitk.GetImageFromArray(arr)
    img.SetMetaData("0010|0010", patient_name)
    img.SetMetaData("0008|0060", "MR")
    img.SetMetaData("0008|0016", "1.2.840.10008.5.1.4.1.1.4")
    
    dicom_file = tmp_path / "donor.dcm"
    writer = sitk.ImageFileWriter()
    writer.SetFileName(str(dicom_file))
    writer.Execute(img)
    return dicom_file


# def test_patient_name_integrated(tmp_path, fake_nifti_image):
#     # Prepare donor and input images
#     donor_dcm = create_fake_donor_dicom(tmp_path, patient_name="John^Doe")
#     nifti_file = tmp_path / "input.nii.gz"
#     sitk.WriteImage(fake_nifti_image, str(nifti_file))
#     out_dir = tmp_path / "dicoms"

#     # Patch writeSlices to record metadata instead of writing real files
#     recorded_metadata = []

#     def fake_writeSlices(series_tag_values, new_img, out_dir_arg, i, writer):
#         tag_dict = dict(series_tag_values)
#         recorded_metadata.append(tag_dict)

#     # Run conversion
#     rc = convert_to_dicom.convert_and_write(
#         nifti_input=str(nifti_file),
#         donor_dcm=str(donor_dcm),
#         dcm_output=str(out_dir),
#         verbose=False,
#     )
#     assert rc == 0

#     # Build series_tag_values manually to test actual metadata
#     reader = sitk.ImageFileReader()
#     reader.SetFileName(str(donor_dcm))
#     reader.LoadPrivateTagsOn()
#     reader.ReadImageInformation()
#     new_img = fake_nifti_image
#     assert reader.HasMetaDataKey("0010|0010")
#     series_tag_values = convert_to_dicom.copy_tags(reader, new_img, "Test Series", "1")

#     # Check that Patient Name is present in the tags
#     tag_dict = dict(series_tag_values)
#     assert "0010|0010" in tag_dict
#     assert tag_dict["0010|0010"] == "John^Doe"


@pytest.fixture
def fake_donor_dicom(tmp_path):
    """Create a minimal valid donor DICOM with Patient Name."""
    dicom_file = tmp_path / "donor.dcm"
    
    arr = np.zeros((2, 2), dtype=np.int16)  # minimal size
    img = sitk.GetImageFromArray(arr)
    patient_name = "John^Doe"
    img.SetMetaData("0010|0010", patient_name)           # Patient Name
    img.SetMetaData("0008|0060", "MR")                  # Modality
    img.SetMetaData("0008|0016", "1.2.840.10008.5.1.4.1.1.4")  # SOP Class UID

    writer = sitk.ImageFileWriter()
    writer.SetFileName(str(dicom_file))
    writer.KeepOriginalImageUIDOn()
    writer.Execute(img)

    return dicom_file, patient_name


def test_donor_contains_patient_name(fake_donor_dicom):
    dicom_file, expected_name = fake_donor_dicom

    # Read donor DICOM
    reader = sitk.ImageFileReader()
    reader.SetFileName(str(dicom_file))
    reader.LoadPrivateTagsOn()
    reader.ReadImageInformation()

    tags = {k : reader.GetMetaData(k) for k in reader.GetMetaDataKeys()}

    # Assert that Patient Name tag exists
    assert reader.HasMetaDataKey("0010|0010"), "Patient Name tag missing in donor DICOM"

    # Assert that the value is correct
    patient_name = reader.GetMetaData("0010|0010")
    assert patient_name == expected_name, f"Patient Name is {patient_name}, expected {expected_name}"

    # Optional: print all keys and values for debugging
    for k in reader.GetMetaDataKeys():
        print(f"{k}: {reader.GetMetaData(k)}")


# @pytest.fixture
# def fake_image():
#     """Small 3D image to pass to copy_tags."""
#     arr = np.zeros((2, 2, 2), dtype=np.int16)
#     img = sitk.GetImageFromArray(arr)
#     return img


# def test_copy_tags_copies_patient_name(tmp_path, fake_donor_dicom, fake_image):
#     dicom_path, expected_name = fake_donor_dicom

#     # Read donor DICOM metadata
#     reader = sitk.ImageFileReader()
#     reader.SetFileName(dicom_path)
#     reader.LoadPrivateTagsOn()
#     reader.ReadImageInformation()

#     for k in reader.GetMetaDataKeys():
#         try:
#             v = reader.GetMetaData(k)
#             print(f"{k}: {v}")
#         except Exception as e:
#             print(f"{k}: <Could not read value: {e}>")

#     assert reader.HasMetaDataKey("0010|0010"), "Donor DICOM does not have Patient Name tag"
#     patient_name_from_reader = reader.GetMetaData("0010|0010")
#     assert patient_name_from_reader == expected_name, "Patient Name in donor DICOM does not match expected value"


#     # Call copy_tags
#     series_tag_values = convert_to_dicom.copy_tags(
#         reader=reader,
#         new_img=fake_image,
#         seriesdesc="Test Series",
#         seriesnumber="1"
#     )

#     # Convert to dictionary for easy checking
#     tag_dict = dict(series_tag_values)

#     # Assert Patient Name is copied correctly
#     assert "0010|0010" in tag_dict
#     assert tag_dict["0010|0010"] == expected_name
    

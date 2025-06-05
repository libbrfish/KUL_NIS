#!/usr/bin/env .venv/bin/python

import nibabel as nib
import numpy as np
import argparse
import json
import os
import sys

def compute_nifti_percentiles(nifti_path, mask_path=None):
    try:
        img = nib.nifti1.load(nifti_path)
        data = img.get_fdata()
    except Exception as e:
        print(f"Error loading NIfTI file: {e}", file=sys.stderr)
        sys.exit(1)

    if mask_path:
        try:
            mask_img = nib.nifti1.load(mask_path)
            mask_data = mask_img.get_fdata()

            if mask_data.shape != data.shape:
                print("Error: Mask and image shapes do not match.", file=sys.stderr)
                sys.exit(1)

            mask = mask_data > 0
            data = data[mask]
        except Exception as e:
            print(f"Error loading or applying mask: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        data = data[np.isfinite(data)]

    data_flat = data.flatten()
    percentiles = {f"p{(10 * p)}": float(np.percentile(data_flat, (10 * p))) for p in range(1, 10)}
    percentiles["p99"] = float(np.percentile(data_flat, 99))
    percentiles["p01"] = float(np.percentile(data_flat, 1))
    
    return percentiles

def main():
    parser = argparse.ArgumentParser(description="Compute percentiles (0–100) of a NIfTI image, optionally using a mask.")
    parser.add_argument("nifti_path", help="Path to the input NIfTI file (.nii or .nii.gz)")
    parser.add_argument("-m", "--mask", help="Optional mask NIfTI file to include only specific voxels (non-zero values).")
    parser.add_argument("-o", "--output", help="Output JSON file path. If not provided, prints to stdout.")

    args = parser.parse_args()

    if not os.path.exists(args.nifti_path):
        print(f"Error: File '{args.nifti_path}' not found.", file=sys.stderr)
        sys.exit(1)

    if args.mask and not os.path.exists(args.mask):
        print(f"Error: Mask file '{args.mask}' not found.", file=sys.stderr)
        sys.exit(1)

    percentiles = compute_nifti_percentiles(args.nifti_path, args.mask)

    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(percentiles, f, indent=2)
        except Exception as e:
            print(f"Error writing output file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print(json.dumps(percentiles, indent=2))

if __name__ == "__main__":
    main()

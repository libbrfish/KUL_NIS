#!/bin/sh
'''exec' uv run --script --project "$(dirname -- "$(realpath -- "$0")")" "$0" "$@"
' '''



"""
intensity_match.py

Perform histogram-based intensity matching between a source and reference NIfTI image
using Bernstein polynomial histogram fitting and multi-segment linear remapping.

Usage:
  python intensity_match.py \
    --source subj01_T1.nii.gz \
    --reference template_T1.nii.gz \
    --output subj01_T1_matched.nii.gz \
    [--mask mask.nii.gz]
"""

import os
from typing import cast
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.special import comb
from scipy.signal import savgol_filter


# ======================================================
# --- Step 1: Multi-segment linear intensity mapping ---
# ======================================================
def multi_linear_match(img, src_vals, ref_vals, mask=None):
    """
    Multi-segment linear intensity remapping (R version translated to Python).
    """
    src_vals = np.asarray(src_vals, dtype=float)
    ref_vals = np.asarray(ref_vals, dtype=float)

    if len(src_vals) != len(ref_vals):
        raise ValueError("Source and reference vectors must have the same length.")
    if len(src_vals) < 2:
        raise ValueError("At least two points are required for interpolation.")
    if np.any(np.diff(src_vals) <= 0):
        raise ValueError("Source values must be strictly increasing.")

    img = img.astype(np.float32)
    out_img = np.copy(img)

    if mask is not None:
        if mask.shape != img.shape:
            raise ValueError("Mask dimensions must match image dimensions.")
        mask = mask != 0
    else:
        mask = np.ones_like(img, dtype=bool)

    # Below first breakpoint
    below_mask = (img < src_vals[0]) & mask
    out_img[below_mask] = (
        ((img[below_mask] - src_vals[0]) / (src_vals[1] - src_vals[0]))
        * (ref_vals[1] - ref_vals[0]) + ref_vals[0]
    )

    # Between breakpoints
    for i in range(len(src_vals) - 1):
        interval_mask = (img >= src_vals[i]) & (img < src_vals[i + 1]) & mask
        out_img[interval_mask] = (
            ((img[interval_mask] - src_vals[i]) / (src_vals[i + 1] - src_vals[i]))
            * (ref_vals[i + 1] - ref_vals[i]) + ref_vals[i]
        )

    # Above last breakpoint
    above_mask = (img >= src_vals[-1]) & mask
    out_img[above_mask] = (
        ((img[above_mask] - src_vals[-2]) / (src_vals[-1] - src_vals[-2]))
        * (ref_vals[-1] - ref_vals[-2]) + ref_vals[-2]
    )

    return out_img


# ======================================================
# --- Step 2: Bernstein polynomial histogram analysis ---
# ======================================================

def bernstein_basis(x, n, k):
    return comb(n, k) * (x ** k) * ((1 - x) ** (n - k))


def bernstein_basis_second_derivative(x, n, k):
    def safe_basis(x, n, k):
        if k < 0 or k > n:
            return 0.0
        return bernstein_basis(x, n, k)
    n2 = n * (n - 1)
    return n2 * (
        safe_basis(x, n - 2, k - 2)
        - 2 * safe_basis(x, n - 2, k - 1)
        + safe_basis(x, n - 2, k)
    )


def bernstein_fit(x, y, degree=8, reg_lambda=0.0):
    x_min, x_max = x.min(), x.max()
    x_norm = (x - x_min) / (x_max - x_min)

    B = np.array([bernstein_basis(x_norm, degree, k) for k in range(degree + 1)]).T
    A = B.T @ B + reg_lambda * np.eye(degree + 1)
    b = B.T @ y
    coeffs = np.linalg.solve(A, b)

    def model(x_query):
        xq = np.clip((x_query - x_min) / (x_max - x_min), 0, 1)
        Bq = np.array([bernstein_basis(xq, degree, k) for k in range(degree + 1)]).T
        return Bq @ coeffs

    def second_derivative(x_query):
        xq = np.clip((x_query - x_min) / (x_max - x_min), 0, 1)
        B2 = np.array([bernstein_basis_second_derivative(xq, degree, k)
                       for k in range(degree + 1)]).T
        scale = 1 / (x_max - x_min)
        return B2 @ coeffs * (scale ** 2)

    return model, second_derivative, (x_min, x_max)


def extract_histogram(data, bins=512, percentile=99.8):
    flattened = data[np.isfinite(data)]
    flattened = flattened[flattened > 0]
    p = np.percentile(flattened, percentile)
    filtered = flattened[flattened <= p]
    hist, edges = np.histogram(filtered, bins=bins, density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, hist, p


def analyze_histogram(nifti_path, degree=8, reg_lambda=0.0, save_plot=False):
    data = cast(nib.Nifti1Image, nib.load(nifti_path)).get_fdata()
    x, y, p99 = extract_histogram(data)

    model, second_deriv_func, (x_min, x_max) = bernstein_fit(x, y, degree, reg_lambda)
    y_fit = model(x)
    y_smooth = savgol_filter(y_fit, window_length=15, polyorder=3)

    # Find mode
    mode = x[np.argmax(y_smooth)]

    # Find max 2nd derivative to the right of mode
    x_vals = np.linspace(x_min, x_max, 2000)
    ddy_vals = second_deriv_func(x_vals)

    # Only look to the right of the mode
    mask = x_vals > mode
    if not np.any(mask):
        return None, None

    x_right = x_vals[mask]
    ddy_right = ddy_vals[mask]

    # Find indices where the sign of the second derivative changes (crosses zero)
    sign_change = np.where(np.diff(np.sign(ddy_right)) != 0)[0]
    if len(sign_change) == 0:
        return None, None

    # First zero crossing (closest to the mode)
    i0 = sign_change[0]
    # Linear interpolation between the points around zero crossing
    x0 = x_right[i0]
    x1 = x_right[i0 + 1]
    y0 = ddy_right[i0]
    y1 = ddy_right[i0 + 1]
    x_zero = x0 - y0 * (x1 - x0) / (y1 - y0)

    max_sd_x = x_zero

    
    # x_vals = np.linspace(x_min, x_max, 1000)
    # ddy_vals = second_deriv_func(x_vals)
    # mask_right = x_vals > mode
    # x_right = x_vals[mask_right]
    # ddy_right = ddy_vals[mask_right]
    # max_sd_x = x_right[np.argmax(ddy_right)] if len(x_right) > 0 else None

    if save_plot:
        plt.figure(figsize=(10, 5))
        plt.plot(x, y, label="Histogram", alpha=0.5)
        plt.plot(x, y_fit, label="Bernstein Fit", color="orange")
        plt.axvline(mode, color="red", ls="--", label=f"Mode: {mode:.2f}")
        if max_sd_x is not None:
            plt.axvline(max_sd_x, color="blue", ls="--", label=f"Shoulder: {max_sd_x:.2f}")
        plt.title(os.path.basename(nifti_path))
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.splitext(nifti_path)[0] + "_histfit.png")
        plt.close()

    return {
        "mode": float(mode),
        "max_second_deriv": float(max_sd_x) if max_sd_x else float(p99),
        "p99": float(p99)
    }


# ======================================================
# --- Step 3: Full pipeline ---
# ======================================================

def match_intensities(source_path, reference_path, output_path, mask_path=None):
    print(f"Analyzing source: {source_path}")
    src_landmarks = analyze_histogram(source_path, save_plot=True)

    print(f"Analyzing reference: {reference_path}")
    ref_landmarks = analyze_histogram(reference_path, save_plot=True)

    # Align corresponding intensity landmarks
    src_vals = [src_landmarks["mode"], src_landmarks["max_second_deriv"], src_landmarks["p99"]]
    ref_vals = [ref_landmarks["mode"], ref_landmarks["max_second_deriv"], ref_landmarks["p99"]]

    print("Source landmarks:", src_vals)
    print("Reference landmarks:", ref_vals)

    # Load NIfTI data
    src_img = cast(nib.Nifti1Image, nib.load(source_path))
    src_data = src_img.get_fdata()
    mask_data = cast(nib.Nifti1Image, nib.load(mask_path)).get_fdata() if mask_path else None

    print("Performing intensity remapping...")
    matched_data = multi_linear_match(src_data, src_vals, ref_vals, mask_data)

    print(f"Saving output image → {output_path}")
    nib.save(
        nib.Nifti1Image(matched_data.astype(np.float32), src_img.affine, src_img.header),
        output_path
    )

    print("Intensity matching complete.")


# ======================================================
# --- CLI Entrypoint ---
# ======================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Intensity match source NIfTI to reference using Bernstein histogram landmarks.")
    parser.add_argument("--source", required=True, help="Source NIfTI file")
    parser.add_argument("--reference", required=True, help="Reference NIfTI file")
    parser.add_argument("--output", required=True, help="Output NIfTI file")
    parser.add_argument("--mask", help="Optional binary mask NIfTI file")

    args = parser.parse_args()

    match_intensities(args.source, args.reference, args.output, args.mask)

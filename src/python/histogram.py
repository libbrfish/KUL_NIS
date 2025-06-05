#!/usr/bin/env .venv/bin/python

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import os
from scipy.special import comb
from scipy.signal import savgol_filter


def load_nifti(path):
    img = nib.nifti1.load(path)
    return img.get_fdata()


def extract_histogram(data, bins=512):
    flattened = data.flatten()
    flattened = flattened[~np.isnan(flattened)]
    flattened = flattened[flattened > 0]
    p99 = np.percentile(flattened, 99)
    filtered = flattened[flattened <= p99]
    hist, edges = np.histogram(filtered, bins=bins, density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, hist, p99


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


def bernstein_fit(x, y, degree, reg_lambda=0.0):
    x_min, x_max = x.min(), x.max()
    x_norm = (x - x_min) / (x_max - x_min)

    B = np.array([bernstein_basis(x_norm, degree, k) for k in range(degree + 1)]).T

    n_coeffs = degree + 1
    A = B.T @ B + reg_lambda * np.eye(n_coeffs)
    b = B.T @ y
    coeffs = np.linalg.solve(A, b)

    def model(x_query):
        xq_norm = (x_query - x_min) / (x_max - x_min)
        xq_norm = np.clip(xq_norm, 0, 1)
        Bq = np.array([bernstein_basis(xq_norm, degree, k) for k in range(degree + 1)]).T
        return Bq @ coeffs

    def second_derivative(x_query):
        xq_norm = (x_query - x_min) / (x_max - x_min)
        xq_norm = np.clip(xq_norm, 0, 1)
        B2 = np.array([bernstein_basis_second_derivative(xq_norm, degree, k) for k in range(degree + 1)]).T
        scale = 1 / (x_max - x_min)
        return B2 @ coeffs * (scale ** 2)

    return model, second_derivative, (x_min, x_max)


def analyze_extrema(x, y_smooth):
    dy = np.gradient(y_smooth, x)
    ddy = np.gradient(dy, x)

    mode_index = np.argmax(y_smooth)
    mode = x[mode_index]

    analysis = {
        "mode": float(mode),
    }
    return analysis, ddy


def find_max_second_deriv_right_of_mode(second_deriv_func, mode, x_min, x_max, num_points=1000):
    x_vals = np.linspace(x_min, x_max, num_points)
    ddy_vals = second_deriv_func(x_vals)
    mask = x_vals > mode
    if not np.any(mask):
        return None, None
    x_right = x_vals[mask]
    ddy_right = ddy_vals[mask]
    max_index = np.argmax(ddy_right)
    return x_right[max_index], ddy_right[max_index]


def main(nifti_path, degree=8, reg_lambda=0.0, output_path="bernstein_fit.png"):
    x, y, p99 = extract_histogram(load_nifti(nifti_path))
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    model, second_deriv_func, (x_min, x_max) = bernstein_fit(x, y, degree, reg_lambda)
    y_fit = model(x)

    # Smooth fit for mode and shoulder detection
    y_smooth = savgol_filter(y_fit, window_length=15, polyorder=3)
    analysis, _ = analyze_extrema(x, y_smooth)

    # Find max second derivative to the right of mode
    max_sd_x, max_sd_val = find_max_second_deriv_right_of_mode(second_deriv_func, analysis["mode"], x_min, x_max)
    analysis['max_sd'] = float(max_sd_x)
    analysis['p99'] = float(p99)
    
    # Plot histogram + Bernstein fit + mode + shoulders + max second derivative
    plt.figure(figsize=(12, 6))
    plt.plot(x, y, label="Original Histogram", alpha=0.4)
    plt.plot(x, y_fit, label=f"Bernstein Fit (deg={degree}, λ={reg_lambda})", color="darkorange", lw=2)
    plt.axvline(analysis["mode"], color="red", linestyle="--", label=f"Mode: {analysis['mode']:.2f}")
    if max_sd_x is not None:
        plt.scatter([max_sd_x], [model(max_sd_x)], color="blue", s=100, marker="x", label=f"Max 2nd Deriv (right of mode): {max_sd_x:.2f}")
    plt.title("Histogram Fit using Bernstein Polynomial")
    plt.xlabel("Intensity")
    plt.ylabel("Density")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

    # Plot analytic second derivative + max second derivative point
    ddy = second_deriv_func(x)
    plt.figure(figsize=(12, 4))
    plt.plot(x, ddy, label="Analytic Second Derivative of Bernstein Fit", color="blue")
    plt.axhline(0, color="black", lw=0.8, linestyle="--")
    if max_sd_x is not None:
        plt.scatter([max_sd_x], [max_sd_val], color="red", s=100, marker="x", label=f"Max 2nd Deriv (right of mode): {max_sd_x:.2f}")
    plt.title("Second Derivative of Bernstein Polynomial Fit")
    plt.xlabel("Intensity")
    plt.ylabel("Second Derivative")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    second_deriv_path = output_path.replace(".png", "_second_derivative.png")
    plt.savefig(second_deriv_path)
    plt.close()

    print("Bernstein fit analysis:")
    print(analysis)
    if max_sd_x is not None:
        print(f"Max second derivative right of mode at x={max_sd_x:.4f} with value={max_sd_val:.6f}")
    else:
        print("No max second derivative point found to the right of the mode.")
    print(f"Plots saved to: {output_path} and {second_deriv_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Bernstein polynomial fit to histogram of NIfTI image intensities.")
    parser.add_argument("nifti_image", help="Path to NIfTI image file (.nii or .nii.gz)")
    parser.add_argument("--degree", type=int, default=8, help="Degree of Bernstein polynomial")
    parser.add_argument("--reg_lambda", type=float, default=0.0, help="Regularization parameter lambda")
    parser.add_argument("--output", default="bernstein_fit.png", help="Output PNG file path")

    args = parser.parse_args()
    main(args.nifti_image, args.degree, args.reg_lambda, args.output)

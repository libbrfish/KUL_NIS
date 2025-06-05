#!/usr/bin/env Rscript

library(RNifti)
library(optparse)

multi_linear_match <- function(img, src_vals, ref_vals, mask = NULL) {
  if (length(src_vals) != length(ref_vals)) {
    stop("Source and reference vectors must have the same length.")
  }
  if (length(src_vals) < 2) {
    stop("At least two points are required for linear interpolation.")
  }

  out_img <- img

  # If a mask is provided, only modify masked voxels
  if (!is.null(mask)) {
    if (!all(dim(mask) == dim(img))) stop("Mask dimensions must match image dimensions.")
    mask <- mask != 0  # Convert to logical
  } else {
    mask <- array(TRUE, dim = dim(img))
  }

  # Below first breakpoint (extrapolate)
  mask_below <- (img < src_vals[1]) & mask
  out_img[mask_below] <- ((img[mask_below] - src_vals[1]) / (src_vals[2] - src_vals[1])) *
                         (ref_vals[2] - ref_vals[1]) + ref_vals[1]

  # Intervals
  for (i in seq_len(length(src_vals) - 1)) {
    m <- (img >= src_vals[i]) & (img < src_vals[i + 1]) & mask
    out_img[m] <- ((img[m] - src_vals[i]) / (src_vals[i + 1] - src_vals[i])) *
                  (ref_vals[i + 1] - ref_vals[i]) + ref_vals[i]
  }

  # Above last breakpoint (extrapolate)
  mask_above <- (img >= src_vals[length(src_vals)]) & mask
  out_img[mask_above] <- ((img[mask_above] - src_vals[length(src_vals) - 1]) /
                          (src_vals[length(src_vals)] - src_vals[length(src_vals) - 1])) *
                          (ref_vals[length(ref_vals)] - ref_vals[length(ref_vals) - 1]) +
                          ref_vals[length(ref_vals) - 1]

  return(out_img)
}

# CLI Setup
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input NIfTI file"),
  make_option(c("-o", "--output"), type="character", help="Output NIfTI file"),
  make_option(c("-s", "--src"), type="character", help="Source percentiles (comma-separated)"),
  make_option(c("-r", "--ref"), type="character", help="Reference percentiles (comma-separated)"),
  make_option(c("-m", "--mask"), type="character", default=NULL,
              help="Optional binary mask NIfTI file (masked voxels remain unchanged)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Safety check
src_clean <- gsub('"', '', opt$src)  # remove any embedded quotes
src_vals <- as.numeric(strsplit(src_clean, ",")[[1]])
ref_clean <- gsub('"', '', opt$ref)  # remove any embedded quotes
ref_vals <- as.numeric(strsplit(ref_clean, ",")[[1]])

if (any(is.na(src_vals)) || any(is.na(ref_vals))) {
  stop("Could not parse --src or --ref values. Ensure they're numeric and comma-separated, e.g., --src \"1,2,3\"")
}

if (length(src_vals) != length(ref_vals)) {
  stop("Source and reference percentiles must be of the same length.")
}

if (any(diff(src_vals) <= 0)) {
  stop("Source percentiles must be strictly increasing.")
}

img <- readNifti(opt$input)
mask <- if (!is.null(opt$mask)) readNifti(opt$mask) else NULL

matched_img <- multi_linear_match(img, src_vals, ref_vals, mask)
writeNifti(matched_img, opt$output)

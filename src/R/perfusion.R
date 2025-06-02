# install.packages("remotes")
remotes::install_github("jonclayden/RNifti")
install.packages("tidyverse")

library(RNifti)
library(tidyverse)

options(warn = 0, error = NULL)

setwd("/media/KUL/DATA2/devos/Projects/DSC/resources/nifti/DSC/")

image <- readNifti(file = "corr_DSC.nii.gz")
# Make sure there are no values smaller than 0
image <- (image - min(image) + 1)

dim(image)
pixunits(image)

TE <- 25

## Baseline image consisting of the first 20 timepoints of the 4D dataset, prior
## to influx of contrast
baseline_img <- image[, , , 1:20]
## Average of 4D time series prior to contrast influx
baseline_avg <- asNifti(apply(baseline_img, MARGIN = c(1, 2, 3), FUN = mean),
                        reference = image)

## R2star
R2 <- image
for (timepoint in seq(dim(image)[4])) {
  R2[, , , timepoint] <- (-1 / TE * (log(R2[, , , timepoint] / baseline_avg)))
}
writeNifti(R2, "R2star.nii.gz", template = "corr_DSC.nii.gz")

## Last 20 timepoints of the 4D dataset, after contrast influx
last_20_img <- image[, , , 60:80]
## Noise standard deviation of baseline
baseline_sd <- apply(baseline_img, MARGIN = c(1, 2, 3), FUN = sd)
## Average of last 20 timepoint after contrast influx
last_20_avg <- apply(last_20_img, MARGIN = c(1, 2, 3), FUN = mean)

## Write nifti files
writeNifti(baseline_avg, "R_avg.nii.gz", template = "corr_DSC.nii.gz")
writeNifti(baseline_sd, "R_sd.nii.gz", template = "corr_DSC.nii.gz")

## Create brain mask of average signal 10x above noise level
mask <- (baseline_avg - (10 * baseline_sd) > 0) |>
  apply(MARGIN = c(1, 2, 3),
        FUN = as.numeric)
writeNifti(mask, "R_mask.nii.gz", template = "corr_DSC.nii.gz")


## Better mask using hd-bet
# install.packages("processx")
library(processx)
test <- run("hd-bet", c("-i", "R_avg.nii.gz", "-s", "1"), echo = TRUE)
file.remove("R_avg_bet.nii.gz")
mask <- readNifti(file = "R_avg_bet_mask.nii.gz")

## Make a mask of non-enhancing voxels. Enhancement is defined as all voxels who
## in the final time steps of the 4D series have a higher signal than the
## baseline + 1 x sd
non_enhancing_mask <- ((last_20_avg - baseline_avg - baseline_sd) < 0) |>
  apply(
    MARGIN = c(1, 2, 3),
    FUN = as.numeric
  )
non_enhancing_mask <- (non_enhancing_mask * mask) |>
  asNifti(reference = image)

## Average all the non-enhancing voxels to remove the spatial dependency.
mean_R2_of_timepoint <- function(x) {
  mean(R2[, , , x][non_enhancing_mask > 0])
}

deltaR2_star_avg <- tibble(timepoint = seq(1:80),
                           intensity = sapply(seq(80), mean_R2_of_timepoint))

ggplot(deltaR2_star_avg, aes(timepoint, intensity)) +
  geom_path()

## Calculate int_0^t{deltaR2_star_avg (t)dt'} :

trapezoidal_int <- function(l) {
  ## Start by dividing the first and last element by 2
  l[[1]] <- (l[[1]] / 2)
  l[[length(l)]] <- (l[[length(l)]] / 2)
  ## The integral is equal to the sum
  sum(l)
}

int_deltaR2_helper<- function(t) {
  deltaR2_star_avg[0:t, "intensity"] |> trapezoidal_int()
}

int_deltaR2 <- tibble(
  timepoint = seq(1:80),
  value = sapply(seq(80), int_deltaR2_helper)
)
ggplot(int_deltaR2, aes(timepoint, value)) +
  geom_path()

## Create the design matrix:
K2_design_matrix <- cbind(deltaR2_star_avg$intensity, int_deltaR2$value)
lmfit(K2_design_matrix, R2[20, 20, 20, ])$be[[2]]
lm.fit(K2_design_matrix, R2[20,20,20,])$coefficients[[2]]

library(Rfast)
calculate_K2 <- function(v) {
  #use lmfit from the rfast library to make it faster!
  Rfast::lmfit(K2_design_matrix, v) |>
    (\(x) -x$be[[2]]) ()
}

library(microbenchmark)

microbenchmark(calculate_K2(R2[82, 33, 20, ]))

K2 <- apply(X = R2, MARGIN = c(1, 2, 3), FUN = calculate_K2)
K2_correction <- trapezoidal_int(int_deltaR2$value)

## rCBV berekenen:
## rCBV is estimated by integrating the relaxivity-time curve
rCBV <- apply(R2, MARGIN = c(1, 2, 3), FUN = trapezoidal_int)
(rCBV * mask) |> asNifti() |> writeNifti("rCBV_uncorrected.nii.gz")

rCBV_corrected <- rCBV + K2 * K2_correction
(rCBV_corrected * mask) |>
  asNifti(reference = image) |>
  writeNifti("rCBV_corrected.nii.gz")


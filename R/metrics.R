# This file contains functions that return structural metrics about point clouds.
# Since we are getting structural metrics on small clips we should just use
# the whole cloud. Therefore, each function has the signature
# f(las) -> numeric(1) so we can use it with lidR::cloud_metrics.

library(lidR)
library(Rdimtools)

# See discussion on fractal dimension here:
# https://github.com/r-lidar/lidR/discussions/428. This returns the slope d
# in Eq. 3 of https://spj.science.org/doi/10.34133/remotesensing.0001.
fractal_dimension <- function(las) {
  xyz <- cbind(las$X, las$Y, las$Z)
  
  est.boxcount(xyz)$estdim
}

# Shannon entropy has been proposed as a way to account for scaling issues
# with fractal dimension. See Eq. 4-6 in
# https://spj.science.org/doi/10.34133/remotesensing.0001.
fractal_dimension_shannon <- function(las, return_fit=FALSE) {
  # Determine set of length scales. This should depend on the length scale
  # of the las file in question.
  
  # Scaling factors
  s <- 2^-c(1:6)
  
  # Convert scaling factors to voxel sizes
  las_length_scale <- ext(las)[2] - ext(las)[1]
  stopifnot(las_length_scale > 0)
  vox_res <- las_length_scale * s
  
  # Calculate voxel counts weighted by Shannon entropy
  Ns <- vapply(
    vox_res, 
    function(v) {
      pt_counts <- voxel_metrics(las, ~length(X), res=v)$V1
      pi <- pt_counts / sum(pt_counts)
      Ns_ <- exp(-sum(pi * log(pi)))
      return(Ns_)
    },
    FUN.VALUE=0
  )
  
  # Now do the regression
  fit <- lm(log(Ns) ~ log(1/s))
  
  if (return_fit) {
    return(fit)
  } else {
  # Extract the slope
    return(unname(coef(fit)[2]))
  }
}

rumple_index_las <- function(las) {
  rumple_index(las$X, las$Y, las$Z)
}

# Below functions are inspired by
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14040

canopy_relief_ratio <- function(las) {
  zmean <- mean(las$Z)
  
  (zmean - min(las$Z)) / (max(las$Z) - zmean)
}

canopy_cover <- function(las, cutoff=2) {
  # Percentage of first returns > cutoff
  freturns <- las$Z[las$ReturnNumber == 1]
  
  sum(freturns > cutoff) / length(freturns)
}

foliar_height_diversity <- function(las) {
  entropy(las$Z) * log(max(las$Z))
}

rugosity <- function(las) {
  freturns <- las$Z[las$ReturnNumber == 1]
  
  sd(freturns)
}

outer_canopy_height <- function(las) {
  freturns <- las$Z[las$ReturnNumber == 1]
  
  median(freturns)
}




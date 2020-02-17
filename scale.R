scale_med <- function(mat){
  # mat is nxp; we want to scale the columns.
  n <- nrow(mat)
  p <- ncol(mat)
  
  # Center.
  mat <- sweep(mat, 2, miscTools::colMedians(mat, na.rm=TRUE), '-')
  
  # Compute MAD and check for zero-variance voxels.
  mad <- 1.4826 * miscTools::colMedians(abs(mat), na.rm=TRUE)
  zero_mad <- mad == 0
  if(any(zero_mad)){
    if(all(zero_mad)){
      stop("All voxels are zero-variance.\n")
    } else {
      warning(cat("Warning: ", sum(zero_mad),
                  " zero-variance voxels (out of ", length(zero_mad),
                  "). These will be set to zero for estimation of the covariance.\n", sep=""))
    }
  }
  
  # Scale.
  scale_col <- function(col, v){ return(ifelse(v != 0, col/v, 0)) }
  mat_scaled <- sweep(mat, 2, mad, scale_col)
  return(mat_scaled)
}
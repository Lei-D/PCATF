library(far) # for orthonormalize


genV = function(p, d){
  # generate V from Gaussian
  V = matrix(rnorm(p*d), nrow = p, ncol = d)
  if(d == 1){
    V = orthonormalization(V, basis = FALSE, norm = TRUE) # normalize
  } else {
    V = orthonormalization(V, basis = TRUE, norm = TRUE)[, 1:d] # orthonormalize
  }
  return(V)
}


sigma0 = function(n, p, n2, d, mu1, gamma1, SNR, 
                  V, group, mu0, gamma0){
  V.sum = sum(abs(apply(V, 1, mean)))
  V.square = sum(apply(V, 1, function(x)sqrt(mean(x^2))))
  sigma0 = (n*sqrt(2*pi)*group*abs(mu1-mu0)*V.sum/SNR - 
                (gamma1*n2 + gamma0*(n-n2))*V.square)/n/p
  if (sigma0 < 0) cat("sigma0 negative")
  return(sigma0)
}


SNR = function(n, p, n2, d, mu1, gamma1, sigma0, 
               V, group, mu0, gamma0){
  V.sum = sum(abs(apply(V, 1, mean)))
  V.square = sum(apply(V, 1, function(x)sqrt(mean(x^2))))
  num = n*sqrt(2*pi)*group*abs(mu1-mu0)*V.sum
  denum = sigma0*n*p + (gamma1*n2 + gamma0*(n-n2))*V.square
  SNR = num/denum
  return(SNR)
}


# data generation function
dataModel = function(n, p, d, mu1, gamma1, 
                     outlier.index = NULL, 
                     SNR = NULL, sigma0 = NULL,
                     group, V, mu0, gamma0, seed = NULL){
  m = length(outlier.index)
  
  # generate index for outliers
  index = rep(FALSE, n)
  index[outlier.index] = TRUE
  
  # set seed in simulation when we want to keep U the same
  if (is.null(seed) == FALSE) {
    set.seed(seed)
  }
  
  # generate rows in U from mixture Gaussian
  U = matrix(NA, nrow = n, ncol = d)
  U[index==TRUE, 1:d] = rnorm(n = d*m, mean = mu1, sd = gamma1)
  U[index==FALSE, 1:d] = rnorm(n = d*(n-m), mean = mu0, sd = gamma0)
  
  # extract V in right dimension
  if (ncol(V) != d) {
    V = V[, 1:d, drop = FALSE]
  }
  
  # calculate SNR or sigma0 when input is the other
  if(is.null(sigma0) == FALSE & is.null(SNR) == TRUE){
    SNR = SNR(n = n, p = p, n2 = m, d = d, 
              mu1 = mu1, gamma1 = gamma1, 
              sigma0 = sigma0, V = V, group = group, 
              mu0 = mu0, gamma0 = gamma0)
  } else if (is.null(sigma0) == TRUE & is.null(SNR) == FALSE) {
    sigma0 = sigma0(n = n, p = p, n2 = m, d = d, 
                    mu1 = mu1, gamma1 = gamma1, 
                    SNR = SNR, V = V, group = group, 
                    mu0 = mu0, gamma0 = gamma0)
  }
  
  # generate X
  X = U %*% t(V) + sigma0 * matrix(rnorm(n*p), n, p)
  
  # extract population mean
  mu = rep(mu0, n)
  mu[index] = mu1
  
  output = list(X = X, U = U, V = V, 
                SNR = SNR,
                sigma0 = sigma0, 
                index = index,
                mu = mu)
  return(output)
}


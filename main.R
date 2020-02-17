library(glmgen)
library(far) # for orthonormalize

# PCA leverage with trend filtering
TF = function(X, K, lambda, niter.max = 1000, tol = 1e-8){
  n = nrow(X)
  p = ncol(X)
  U = matrix(NA, nrow = n, ncol = K)
  V = matrix(NA, nrow = p, ncol = K)
  D = rep(NA, K)
  for(k in 1:K){
    decomp = svd(X)
    u = decomp$u[, k]
    v = decomp$v[, k]
    d = decomp$d[k]
    for(i in 1:niter.max){
      u.last = u
      u = glmgen::trendfilter(y = scale(X %*% v, center = FALSE, scale = d), 
                              x = 1:nrow(X), 
                              k = 0, lambda = lambda)$beta
      v = orthonormalization(crossprod(X, u) , basis = FALSE, norm = TRUE) # normalize
      diff = sqrt(mean((u - u.last)^2))
      # print(i)
      # print(diff)
      if(diff < tol){
        break
      }
    }
    # print(paste("niter", i))
    d = crossprod(u, X %*% v)[1, 1]
    X = X - d * tcrossprod(u, v)
    U[, k] = u
    V[, k] = v
    D[k] = d
  }
  output = list(U = U, V = V, D = D)
  return(output)
}

# compare estimated U for population mean, PCA leverage and our method
compare.methods = function(data, lambda, d){
  n = nrow(data$U)
  
  # population mean
  mu = matrix(rep(data$mu, d), nrow = n, ncol = d, byrow = FALSE)
  mu.colnorm = apply(mu, 2, function(x)sqrt(sum(x^2)))
  U = scale(mu, center = FALSE, scale = mu.colnorm)
  # U = orthonormalization(mu, basis = FALSE, norm = TRUE) # normalize
  pop = diag(tcrossprod(U))
  
  # PCA leverage
  svd.sample = svd(scale_med(data$X))
  PCAleverage = diag(tcrossprod(svd.sample$u[, 1:d, drop = FALSE]))
  
  df = data.frame(method = c(rep("Population Mean",n), 
                             rep("PCA Leverage", n)),
                  index = c(1:n, 1:n),
                  l2norm = c(pop, PCAleverage),
                  outlier = rep(data$index, 2), 
                  stringsAsFactors = FALSE)
  
  # trend filtering
  for(i in 1:length(lambda)) {
    our = TF(X = scale_med(data$X), K = d, lambda = lambda[i])
    our.row.norm = diag(tcrossprod(our$U))
    df.new = data.frame(method = rep(paste("lambda", round(lambda[i], 4)), n),
                        index = 1:n,
                        l2norm = our.row.norm,
                        outlier = data$index, 
                        stringsAsFactors = FALSE)
    df = bind_rows(df, df.new)
  }
  # df$method = factor(df$method, levels = unique(df$method))
  
  return(df)
}

id_out.leverage = function(leverage){
  cutoffs = 3:5 * median(leverage)
  names(cutoffs) = paste0(as.character(3:5), 'med')
  out.lev3 = (leverage > cutoffs[[1]])
  out.lev4 = (leverage > cutoffs[[2]])
  out.lev5 = (leverage > cutoffs[[3]])
  
  out = data.frame(out.lev3, out.lev4, out.lev5)
  names(out) = c('3 x median','4 x median','5 x median')
  result = list(outliers = out, cutoffs = cutoffs)
  return(result)
}

our = function(X, nPC, lambda, id_out = TRUE, scale = FALSE) {
  # Center and scale robustly.
  if(scale == TRUE){
    X = scale_med(X)
  }
  
  our = TF(X = X, K = nPC, lambda = lambda)
  U = our$U
  measure = diag(tcrossprod(U))
  
  params = list(choosePCs = "fixPC", method = "leverage")
  result = list(params = params, PCs = U, leverage = measure, robdist = NULL,
                inMCD = NULL, outliers = NULL, cutoffs = NULL)
  
  # Label outliers, if requested.
  if(id_out){
    id_out = id_out.leverage(measure)
    result$outliers = id_out$outliers
    result$cutoffs = id_out$cutoffs
  }
  
  class(result) = c('clever', class(result))
  return(result)
}
library(isotone)

Gbeta <- function(beta, ecdf_Y, X) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  sub_vals <- X %*% beta
  zhat <- as.vector(sapply(ecdf_Y, inverse_cdf_z, arg=sub_vals))
  
  res <- t(X) %*% (zhat - X %*% beta)
  return (res)
}
# FIXME values of matrix in Z are too large
Gbeta_grad <-function(beta, ecdf_Y, X, pert=1e-4) {
  n <- nrow(X)
  m <- ncol(X)
  
  X_bar <- X^(-1)
  
  # print(X_bar)
  
  Z_bar <- matrix(0, nrow=n, ncol=n)
  sub_vals <- X %*% beta
  zhat <- as.vector(sapply(ecdf_Y, inverse_cdf_z, arg=sub_vals))
  
  for(i in 1:n) {
    for(j in 1:n) {
      Z_bar[i, j] <- 1/(dnorm(zhat[j] - t(X[i, ]) %*% beta) + pert)
    }
  }
  print(Z_bar)
  grad_G_beta <- -n*t(X) %*% (Z_bar %*% X_bar - X)
  
  return (grad_G_beta)
}

compute_G_beta <- function(n) {
  data <- simulate_rank_regression_data(n, 6)
  
  X <- data$X
  Y <- data$Y
  
  ranks_Y <- rank(Y)
  ecdf_Y <- ranks_Y/(n+1)
  ecdf_Y
  
  beta <- data$beta
  g_beta <- Gbeta(ecdf_Y, X, beta)
  
  res <- list("gt_beta"=beta, "G_beta"=g_beta)
  return (res)
}

find_f1_coefs_fixed_point <- function(Y, X, tol=1e-9, max_iter=300) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  A <- t(X) %*% X
  P <- inv(A) %*% t(X)
  
  m <- ncol(X)
  n <- nrow(X)
  coefs <- runif(m, min=-10, max=10)
  
  ranks_Y <- rank(Y)
  empirical_cdf_Y <- ranks_Y/(n+1)
  
  err <- 1
  
  for(iter in 1:max_iter) {
    print(paste("iter  ", iter))
    
    sub_vals <- X %*% coefs
    yhat <- as.vector(sapply(empirical_cdf_Y, inverse_cdf_z, arg=sub_vals))
    
    a <- P %*% yhat
    err <- sum((a - coefs)**2)**0.5
    
    coefs <- a
    
    print(paste("error is ", err))
    print("coefs are ")
    print(coefs)
    if(err < tol) {
      print("converged")
      break
    }
  }
  
  return(coefs)
}

# https://arxiv.org/pdf/2103.13435.pdf
# Note computationally expensive, works 100 sample size
# FIXME some times produce variance zero
# FIXME check the sign of beta 
# FIXME change alg for Gausssian noise
rank.reg.prl.dim1 <-function(Y, X) {
  prl.F <- function(beta, Y, X) {
    n <- length(Y)
    I <- matrix(Y, n, n) - matrix(Y, n, n, byrow = T)
    I[I < 0] <- 0
    I[I > 0] <- 1
    I <- I[row(I) != col(I)]
    
    m_X <- X %*% beta
    v <- matrix(m_X, n, n) - matrix(m_X, n, n, byrow = T)
    v <- v[row(v) != col(v)]
    
    res_pava <- gpava(v, I)
    
    return(res_pava)
  }
  
  res_pava <- prl.F(1, Y, X)
  
  sort_x <- res_pava$x[order(res_pava$z)]
  sort_z <- res_pava$z[order(res_pava$z)]
  
  var_eps <- 0
  for(i in 1:(length(sort_x)-1)) {
    term <- (((sort_z[i] + sort_z[i+1])/2)^2)*(sort_x[i+1] - sort_x[i])
    var_eps <- var_eps + term
  }
  
  var_eps <- var_eps/sqrt(2)
  est_beta <- 1/var_eps
  
  return(est_beta)
}
cdf_z <- function(y, subtract_values) {
  res = 0
  for(i in subtract_values) {
    res = res + pnorm(y - i)
  }
  res = res / length(subtract_values)
  return(res)
}

inverse <- function(f, lower, upper){
  function(y, arg){
    if(f(lower, arg) - y > 0) {
      return(lower)
    }
    if (f(upper, arg) - y < 0) {
      return(upper)
    }
    uniroot(function(x){f(x, arg) - y}, lower = lower, upper = upper, tol=1e-3,)[1]$root
  }
}

inverse_cdf_z <- inverse(cdf_z, -10000000, 10000000)

find_f1_coefs_fixed_point_stochastic <- function(Y, X, batch_size=64, 
                                                 max_iter=100, tol=1e-9) {
  print("starting fixed point algorithm")
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  
  m <- ncol(X)
  n <- nrow(X)
  coefs <- matrix(runif(m, min=-10, max=10), m, 1)
  # coefs <- beta
  
  err <- 1
  for(iter in 1:max_iter) {
    random_inds <- sample(1:n, batch_size)
    batch_X <- X[random_inds, ]
    batch_Y <- Y[random_inds, ]
    
    A <- t(batch_X) %*% batch_X
    if(nrow(A) > 1) {
      P <- solve(A) %*% t(batch_X)
    } else {
      P <- 1/A[1, 1] * t(batch_X)
    }
    ranks_Y <- rank(batch_Y)
    empirical_cdf_Y <- ranks_Y/(batch_size+1)
    
    # print(paste("iter for fixed point alg  ", iter))
    sub_vals <- batch_X %*% coefs
    yhat <- as.vector(sapply(empirical_cdf_Y, inverse_cdf_z, arg=sub_vals))
    
    a <- P %*% yhat
    err <- sum((a - coefs)**2)**0.5
    
    coefs <- a
    
    # print(paste("error is ", err))
    # print("coefs are ")
    # print(coefs)
    if(err < tol) {
      print("converged")
      break
    }
  }
  
  return(coefs)
}

# possible penalties are ell1 and ell2
find_f1_coefs_expected_rank_algorithm <- function(Y, X, lamb=10, penalty="ell2") {
  print(paste("starting expected rank ", penalty))
  G_j_beta <- function(j, beta, X) {
    val <- 0
    for(i in 1:nrow(X)) {
      val <- val + pnorm(sum((X[j, ] - X[i, ])*beta)/2**0.5)
    }
    return (1/2 + val)
  }
  
  S_beta <- function(beta, X, ranks_Y) {
    val <- 0
    for(j in 1:nrow(X)) {
      val <- val + (ranks_Y[j] - G_j_beta(j, beta, X))**2
    }
    
    if (penalty=="ell1") {
      val <- val + lamb*sum(abs(beta))
    } else if(penalty=="ell2") {
      val <- val + lamb*sum((beta)**2)
    } else {
      print("no penalty")
    }
    
    return (val)
  }
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  
  m <- ncol(X)
  n <- nrow(X)
  coefs <- matrix(runif(m, min=-10, max=10), m, 1)
  ranks_Y <- rank(Y)
  
  est_beta <- optim(par=coefs, fn=S_beta, method = "BFGS", X=X, ranks_Y=ranks_Y)
  
  return(est_beta$par)
}
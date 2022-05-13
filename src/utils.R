set.seed(12)

comute.RB.Var.MSE <- function(estimates, gt_param) {
  # Relative Bias
  rb <- (colMeans(estimates) - gt_param)/gt_param
  
  # Variance
  mean_ests <- matrix(colMeans(estimates), nrow = nrow(estimates), ncol = ncol(estimates), byrow = T)
  vr <- colMeans((estimates - mean_ests)^2)
  
  # MSE
  gt_ests <- matrix(gt_param, nrow = nrow(estimates), ncol = ncol(estimates), byrow = T)
  ms <- colMeans((estimates - gt_ests)^2)
  
  res <- list("RB"=rb, "Var"=vr, "MSE"=ms)
  
  return(res)
}

simulate_rank_regression_data <- function(n, m, beta_max=10) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  
  noise <- rnorm(n)
  X <- matrix(rnorm(n*m), n, m)
  
  beta <- runif(n=m, min=-beta_max, max=beta_max)
  
  # make one third of the coefs zero
  k <- m %/% 3
  zero_inds = sample(1:m, k)
  beta[zero_inds] <- 0
  
  Z <- X %*% beta + noise
  Y <- exponent(Z, 1/3) + 4.7
  
  res <- list("X"=X, "Y"=Y, "beta"=beta, "noise"=noise, "Z"=Z)
  return(res)
}

simulate_rank_regression_data_fixed_beta <- function(n, m, beta) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  
  noise <- rnorm(n)
  X <- matrix(rnorm(n*m), n, m)
  
  Y <- X %*% beta + noise
  Y <- exponent(Y, 1/3) + 4.7
  
  res <- list("X"=X, "Y"=Y, "beta"=beta)
  return(res)
}

# identifiable model
simulate_bivariate_pnl_idf <- function(n) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  X <- matrix(rnorm(n), n, 1)
  beta <- runif(1, -100, 100)
  noise <- rnorm(n)
  Y <- (X^2) %*% beta + noise
  Y <- exponent(Y, 1/3) + 4.7
  data <- cbind(X, Y)
  res <- list("data"=data, "beta"=beta)
  
  return(res)
}

# non identifiable model
simulate_bivariate_pnl <- function(n) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  X <- matrix(rnorm(n), n, 1)
  beta <- runif(1, -100, 100)
  noise <- rnorm(n)
  Y <- X %*% beta + noise
  Y <- exponent(Y, 1/3) + 4.7
  data <- cbind(X, Y)
  res <- list("data"=data, "beta"=beta)
  
  return(res)
}

simulate_mult_pnl <- function(n, fj2_func="cube") {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  
  noise2 = rnorm(n)
  noise3 = rnorm(n)
  noise4 = rnorm(n)
  
  beta2 <- runif(1, -100, 100)
  beta3 <- runif(1, -100, 100)
  beta4 <- runif(1, -100, 100)
  
  x1 <- rnorm(n)
  x2 <- 6.7 + beta2*x1^2 + noise2
  x2 <- exponent(x2, 1/3) + 1
  x3 <- 1.2 + beta3*x1^2 + noise3
  x3 <- exponent(x3, 1/5) + 5
  if(fj2_func == "cube") {
    x4 <- (beta4*(x2 + x3) + noise4)^3
  } else {
    x4 <- beta4*(x2 + x3) + noise4
  }
  
  data <- cbind(x1, x2, x3, x4)
  data <- data.frame(data)
  
  return(data)
}
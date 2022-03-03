source("rank_regression.R")
library(dHSIC)

simulate_bivariate_pnl <- function(n=1000) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  X <- matrix(rnorm(n), n, 1)
  beta <- runif(1, -10, 10)
  noise <- rnorm(n)
  Y <- X %*% beta + noise
  Y <- exponent(Y, 1/3) + 4.7
  data <- cbind(X, Y)
  res <- list("data"=data, "beta"=beta)
  
  return(res)
}

find_estimates <- function(Y, X, rank_alg) {
  est_beta <- NaN
  if(rank_alg=="expected_rank") {
    est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, lamb = 10)
  } else if(rank_alg=="fixed_point") {
    est_beta <- find_f1_coefs_fixed_point_stochastic(Y, X)
  } else {
    print("no such rank regression algorithm. Options are expected_rank and fixed_point")
    return()
  }
  
  ranks_Y <- rank(Y)
  empirical_cdf_Y <- ranks_Y/(length(Y)+1)
  sub_vals <- X %*% est_beta
  print("estimate f_2_inv_y")
  
  # f_2_inv_y <- c()
  # for(e_y in empirical_cdf_Y) {
  #   print(e_y)
  #   f_2_inv_y <- c(f_2_inv_y, inverse_cdf_z(e_y, arg = sub_vals))
  # }
  f_2_inv_y <- as.vector(sapply(empirical_cdf_Y, inverse_cdf_z, arg=sub_vals))
  est_errors <- f_2_inv_y - sub_vals
  
  print("compute dhsic value")
  dhsic_val <- dhsic(X, est_errors, kernel = "gaussian.fixed")$dHSIC
  
  res <- list("est_beta"=est_beta, "est_errors"=est_errors, "dhsic_val"=dhsic_val)
  
  return(res)
}

find_bivariate_direction <- function(data, f_name, rank_alg) {
  X <- data[, 1]
  Y <- data[, 2]
  
  res1 <- find_estimates(Y, X, rank_alg)
  res2 <- find_estimates(X, Y, rank_alg)
  
  est_direction <- ""
  if(res1$dhsic_val < res2$dhsic_val) {
    est_direction <- "1 -> 2"
  } else {
    est_direction <- "2 -> 1"
  }
  
  return(list("res1"=res1, "res2"=res2, "est_direction"=est_direction, "f_name"= f_name))
}

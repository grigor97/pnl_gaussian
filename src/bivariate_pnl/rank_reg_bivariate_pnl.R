source("rank_regression/rank_regression.R")
library(dHSIC)

find_estimates <- function(Y, X, rank_alg, lamb) {
  est_beta <- NaN
  if(rank_alg=="expected_l2_rank") {
    est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, 
                                                      lamb = lamb, penalty = "ell2")
  } else if(rank_alg=="expected_l1_rank") {
    est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, 
                                                      lamb = lamb, penalty = "ell1")
  } else if(rank_alg=="fixed_point") {
    est_beta <- find_f1_coefs_fixed_point_stochastic(Y, X)
  } else {
    print("no such rank regression algorithm. Options are expected_l2_rank, 
          expected_l1_rank and fixed_point")
    return()
  }
  
  ranks_Y <- rank(Y)
  empirical_cdf_Y <- ranks_Y/(length(Y)+1)
  sub_vals <- X %*% est_beta
  # print("estimate f_2_inv_y")
  
  f_2_inv_y <- as.vector(sapply(empirical_cdf_Y, inverse_cdf_z, arg=sub_vals))
  est_errors <- f_2_inv_y - sub_vals
  
  # print("compute dhsic value")
  dhsic_val <- dhsic(X, est_errors, kernel = "gaussian.fixed")$dHSIC
  
  res <- list("est_beta"=est_beta, "est_errors"=est_errors, "dhsic_val"=dhsic_val)
  
  return(res)
}

find_bivariate_direction <- function(data, rank_alg, lamb, file_name="") {
  X <- data[, 1]
  Y <- data[, 2]
  
  XX <- cbind(X, X^2)
  YY <- cbind(Y, Y^2)
  
  res1 <- find_estimates(Y, XX, rank_alg, lamb)
  res2 <- find_estimates(X, YY, rank_alg, lamb)
  
  est_direction <- ""
  if(res1$dhsic_val < res2$dhsic_val) {
    est_direction <- "1 -> 2"
  } else {
    est_direction <- "2 -> 1"
  }
  
  print("dhsic vals ")
  print(res1$dhsic_val)
  print(res2$dhsic_val)
  print(est_direction)
  
  return(list("res1"=res1, "res2"=res2, "est_direction"=est_direction, 
              "file_name"= file_name))
}

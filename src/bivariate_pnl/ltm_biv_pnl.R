source("linear_transofrmation_models/noise_estimation_ltm.R")
library(dHSIC)

min.dhsic.val <- function(Y, X, beta_alg, f2_inv_alg) {
  res_noise_est <- noise.est(Y, X, beta_alg, f2_inv_alg)
  est_noise <- res_noise_est$est_noise
  est_f2_inv <- res_noise_est$est_f2_inv
  est_beta <- res_noise_est$est_beta
  
  # print("compute dhsic value")
  dhsic_val <- dhsic(X, est_noise, kernel = "gaussian.fixed")$dHSIC
  
  res <- list("est_beta"=est_beta, "est_noise"=est_noise, 
              "est_f2_inv"=est_f2_inv, "dhsic_val"=dhsic_val)
  
  return(res)
}

find_bivariate_direction <- function(data, beta_alg, f2_inv_alg) {
  X <- data[, 1]
  Y <- data[, 2]
  
  XX <- cbind(X, X^2)
  YY <- cbind(Y, Y^2)
  
  res1 <- min.dhsic.val(Y, XX, beta_alg, f2_inv_alg)
  res2 <- min.dhsic.val(X, YY, beta_alg, f2_inv_alg)
  
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
  
  return(list("res1"=res1, "res2"=res2, "est_direction"=est_direction))
}

source("linear_transofrmation_models/beta_estimation_ltm.R")
source("linear_transofrmation_models/f2_estimation_ltm.R")

beta.est <- function(Y, X, beta_alg) {
  est_beta <- NaN
  if(beta_alg=="prl") {
    est_beta <- beta.est.prl.gaussian(Y, X)$est_beta
  } else if (beta_alg=="oprl") {
    est_beta <- beta.est.oprl.gaussian(Y, X)$est_beta
  } else if(beta_alg=="expected_l2_rank") {
    est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, 
                                                      lamb = 10, penalty = "ell2")
  } else if(beta_alg=="expected_l1_rank") {
    est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, 
                                                      lamb = 10, penalty = "ell1")
  } else if(beta_alg=="fixed_point") {
    est_beta <- find_f1_coefs_fixed_point_stochastic(Y, X)
  } else {
    print("no such beta estimation algorithm. Options are expected_l2_rank, 
          expected_l1_rank, fixed_point, prl, oprl")
    return()
  }
  
  return(est_beta)
}

f2.inv.est <- function(Y, X, est_beta, f2_inv_alg) {
  est_f2_inv <- NaN
  if(f2_inv_alg=="rl") {
    est_f2_inv <- f2.inv.est.rl(Y, X, est_beta)
  } else if(f2_inv_alg=="rank_reg") {
    est_f2_inv <- f2.inv.est.rank.reg(Y, X, est_beta)
  } else {
    print("no such f2 inverse estimation algorithm. Options are rl and rank_reg")
    return()
  }
  return(est_f2_inv)
}

noise.est <- function(Y, X, beta_alg, f2_inv_alg) {
  est_beta <- beta.est(Y, X, beta_alg)
  print(paste("est beta is  ", est_beta))
  est_f2_inv <- f2.inv.est(Y, X, est_beta, f2_inv_alg)
  est_noise <- est_f2_inv - (X %*% est_beta)
  
  return(list("est_beta"=est_beta, "est_f2_inv"=est_f2_inv, "est_noise"=est_noise))
}

# source("utils.R")
# data <- simulate_rank_regression_data(1000, 1)
# X <- data$X
# Y <- data$Y
# gt_beta <- data$beta
# gt_Z <- data$Z
# gt_noise <- data$noise
# gt_beta
# 
# res <- noise.est(Y, X, "prl", "rank_reg")
# est_beta <- res$est_beta
# est_Z <- res$est_f2_inv
# est_noise <- res$est_noise
# est_beta
# 
# mean((est_noise - gt_noise)**2)
# mean((est_Z - gt_Z)**2)
# library(ggplot2)
# ggplot() + geom_point(aes(x=gt_noise, y=est_noise)) +
#   xlim(c(-4, 4)) + ylim(c(-4, 4))
#   labs(x="true values", y="predicted values")
# 
# ggplot() + geom_point(aes(x=gt_Z, y=est_Z)) +
#   labs(x="true values", y="predicted values")

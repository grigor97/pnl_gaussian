source("multivariate_pnl/full_hsic_pnl.R")
source("multivariate_pnl/rank_reg_pnl.R")
source("utils.R")
library(rjson)
library(parallel)

run_mult_full_hsic <- function(sample_size) {
  data <- simulate_mult_pnl(n=sample_size)
  
  df <- data.frame(data)
  
  res <- order_recovery_full_hsic(df)
  print(res)
  
  return(res)
}

run_mult_rank_reg <- function(sample_size, rank_alg, lamb) {
  data <- simulate_mult_pnl(n=sample_size)
  
  df <- data.frame(data)
  
  res <- order_recovery_rank_reg(df, rank_alg, lamb)
  print(res)
  
  return(res)
}

# res_hsic <- run_mult_full_hsic(100)
res_rank <- run_mult_rank_reg(100, rank_alg="expected_l2_rank", lamb=10)




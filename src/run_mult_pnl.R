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
# res_rank <- run_mult_rank_reg(500, rank_alg="expected_l1_rank", lamb=10)
# res_rank <- run_mult_rank_reg(500, rank_alg="expected_l2_rank", lamb=10)
# res_rank <- run_mult_rank_reg(500, rank_alg="fixed_point", lamb=10)

num_datasets <- 100
sample_size <- 1000
alg <- "fixed_point"

numCores <- detectCores() - 1
print(paste("num coress  --- ", numCores))


dummy_fun <- function(i) {
  if(alg == "full_hsic") {
    return(run_mult_full_hsic(sample_size))
  } else {
    return(run_mult_rank_reg(sample_size, rank_alg=alg, lamb=10))
  }
}

system.time(
  results <- mclapply(1:num_datasets, dummy_fun, mc.cores = numCores)
)
results <- unlist(results)
pred_order <- matrix(results, ncol = 4, nrow=length(results)/4, byrow = T)

res <- list("alg"=alg, "num_datasets"=num_datasets, "sample_size"=sample_size,
            "pred_order"=pred_order)

save_file <- "../res/mult_pnl/results4"
save_file <- paste(save_file, alg, sep="_")
save_file <- paste(save_file, num_datasets, sep="_")
save_file <- paste(save_file, ".json", sep="")

json_data <- toJSON(res)
write(json_data, save_file)





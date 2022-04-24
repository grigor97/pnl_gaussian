source("multivariate_pnl/full_hsic_pnl.R")
source("multivariate_pnl/ltm_mult_pnl.R")
source("utils.R")
library(rjson)
# library(parallel)
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores-1)
registerDoParallel(cl) 

run_mult_full_hsic <- function(sample_size) {
  data <- simulate_mult_pnl(n=sample_size)
  
  df <- data.frame(data)
  
  res <- order_recovery_full_hsic(df)
  print(res)
  
  return(res)
}

run_mult_ltm <- function(sample_size, beta_alg, f2_inv_alg) {
  data <- simulate_mult_pnl(n=sample_size)
  
  df <- data.frame(data)
  
  res <- order_recovery_ltm(df, beta_alg, f2_inv_alg)
  print(res)
  
  return(res)
}

# res_hsic <- run_mult_full_hsic(100)
# res_rank <- run_mult_ltm(100, "prl", "rank_reg")



num_datasets <- 4
sample_size <- 100
alg = 'ltm'
beta_alg <- "prl"
f2_inv_alg <- "rank_reg"

numCores <- detectCores() - 1
print(paste("num coress  --- ", numCores))


dummy_fun <- function(i) {
  if(alg == "full_hsic") {
    return(run_mult_full_hsic(sample_size))
  } else {
    return(run_mult_ltm(sample_size, beta_alg, f2_inv_alg))
  }
}

res <- foreach(i=1:num_datasets, .combine="rbind", .packages = c("dHSIC")) %dopar% {
  dummy_fun(i)
}

res
result <- list("alg"=alg, "num_datasets"=num_datasets, "sample_size"=sample_size,
            "pred_order"=res, "beta_alg"=beta_alg, "f2_inv_alg"=f2_inv_alg)

save_file <- "../res/mult_pnl/results4"
save_file <- paste(save_file, alg, sep="_")
save_file <- paste(save_file, beta_alg, sep="_")
save_file <- paste(save_file, f2_inv_alg, sep="_")
save_file <- paste(save_file, num_datasets, sep="_")
save_file <- paste(save_file, ".json", sep="")

json_data <- toJSON(result)
write(json_data, save_file)

stopCluster(cl)

# res1 <- fromJSON(file = save_file)
# matrix(res1$pred_order, length(res1$pred_order) %/% 4, 4) == res



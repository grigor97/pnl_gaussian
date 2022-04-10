# library(parallel)
library(doParallel)
library(rjson)
source("rank_regression/rank_regression.R")
source("utils.R")

no_cores <- detectCores()
cl <- makeCluster(no_cores-1)
registerDoParallel(cl) 

n <- 1000
m <- 1
betas <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 10, 30, 50, 70, 100, 1000)
num_datasets <- 100
save_file_name <- "../res/rank_regression/"
alg_name <- "normal_scores"
alg <- function(Y, X) {
  res = lin.tr.models.normal.scores(Y, X)
  return(res)
}

run_rank_reg_alg <- function(n, m, betas){
  est_betas <- matrix(0, nrow=m, ncol=length(betas))
  
  for (j in 1:length(betas)) {
    print(paste("beta ", j))
    beta <- betas[j]

    data <- simulate_rank_regression_data_fixed_beta(n, m, beta)
    X <- data$X
    Y <- data$Y
    
    res_alg <- alg(Y, X)
    est_beta <- res_alg$est_beta
    est_betas[ , j] <- est_beta
  }
  return(est_betas)
}

save_result <- function(est_betas, betas, n, m, alg_name, file_name) {
  res <- list("betas"=betas, "est_betas"=est_betas,
              "num_datasets"=nrow(est_betas), "num_betas"=length(betas), 
              "n"=n, "m"=m)
  json_data <- toJSON(res)
  
  file_name <- paste(file_name, "all_betas", sep='')
  file_name <- paste(file_name, alg_name, sep='_')
  file_name <- paste(file_name, n, sep='_')
  file_name <- paste(file_name, m, sep='_')
  file_name <- paste(file_name, nrow(est_betas), sep='_')
  file_name <- paste(file_name, length(betas), sep='_')
  file_name <- paste(file_name, '.json', sep='')
  
  write(json_data, file_name)
  
  return(res)
}

res <- foreach(i=1:num_datasets, .combine="rbind", .packages = c("EnvStats")) %dopar% {
  run_rank_reg_alg(n=n, m=m, betas)
}

res

save_result(est_betas=res, betas=betas, n=n, m=m, 
            alg_name=alg_name, file_name=save_file_name)

stopCluster(cl)

# data <- fromJSON(file = paste(save_file_name, "all_betas_prl_100_1_6_3.json", sep = ""))
# est_betas <- matrix(data$est_betas, data$num_datasets, data$num_betas)
# sum(est_betas - res)

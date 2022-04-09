# library(parallel)
library(doParallel)
library(rjson)
source("rank_regression/rank_regression.R")

no_cores <- detectCores()
cl <- makeCluster(no_cores-1)
registerDoParallel(cl) 

run_rank_reg_alg <- function(n, m, max_iter, batch_size, lamb){
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  betas <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 10, 30, 50, 70, 100, 1000)
  est_betas <- matrix(0, nrow=m, ncol=length(betas))
  
  # generate data
  X <- matrix(rnorm(n*m), n, m)
  for (j in 1:length(betas)) {
    print(paste("beta ", j))
    beta <- betas[j]
    # generate noise and corresponding y
    noise <- rnorm(n)
    Y <- X %*% beta + noise
    Y <- exponent(Y, 1/3) + 4.7
    
    # est_beta <- find_f1_coefs_fixed_point_stochastic(Y, X, batch_size=batch_size,
    #                                                         max_iter=max_iter, tol=1e-9)
    # est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, lamb=lamb)
    # est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, lamb=lamb, penalty = "ell1")
    est_beta <- find_f1_coefs_monte_carlo_algorithm(Y, X)
    # est_beta <- find_f1_coefs_cond_monte_carlo_algorithm(Y, X)

    est_betas[ , j] <- est_beta
  }
  return(est_betas)
}

save_result <- function(est_betas, betas, n, m, lamb, alg_name, file_name) {
  res <- list("betas"=betas, "est_betas"=est_betas,
              "num_datasets"=nrow(est_betas), "num_betas"=length(betas), 
              "n"=n, "m"=m, "lamb"=lamb)
  json_data <- toJSON(res)
  
  file_name <- paste(file_name, "all_betas", sep='')
  file_name <- paste(file_name, alg_name, sep='_')
  file_name <- paste(file_name, lamb, sep='_')
  file_name <- paste(file_name, n, sep='_')
  file_name <- paste(file_name, m, sep='_')
  file_name <- paste(file_name, nrow(est_betas), sep='_')
  file_name <- paste(file_name, length(betas), sep='_')
  file_name <- paste(file_name, '.json', sep='')
  
  write(json_data, file_name)
  
  return(res)
}


res <- foreach(i=1:100, .combine="rbind") %dopar% {
  run_rank_reg_alg(n=1000, m=1, max_iter=100, batch_size=64, lamb=10)
}

res

betas <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 10, 30, 50, 70, 100, 1000)
save_result(est_betas=res, betas=betas, n=1000, m=1, lamb=0, 
            alg_name="mc", file_name="../res/rank_regression/1000/")

# res1 <- fromJSON(file = "../res/rank_regression/1000/all_betas_cond_mc_0_1000_1_100_13.json")
# est_betas <- matrix(res1$est_betas, res1$num_datasets, res1$num_betas)
# sum(est_betas - res)

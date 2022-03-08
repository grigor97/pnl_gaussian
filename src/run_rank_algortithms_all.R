source("rank_regression.R")


run_rank_regression_algorithms <- function(n, m, max_iter, batch_size, lamb){
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  betas <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 10, 30, 50, 70, 100, 1000)
  
  fixed_est_betas <- matrix(0, nrow=m, ncol=length(betas))
  expected_est_betas <- matrix(0, nrow=m, ncol=length(betas))
  expected_l1_est_betas <- matrix(0, nrow=m, ncol=length(betas))
  
  # generate data
  X <- matrix(rnorm(n*m), n, m)
  
  for (j in 1:length(betas)) {
    print(paste("beta ", j))
    beta <- betas[j]
    
    # generate noise and corresponding y
    noise <- rnorm(n)
    Y <- X %*% beta + noise
    Y <- exponent(Y, 1/3) + 4.7
    # print("fixed alg started")
    # fixed_pred_beta <- find_f1_coefs_fixed_point_stochastic(Y, X, batch_size=batch_size,
    #                                                         max_iter=max_iter, tol=1e-9)
    # print("exp alg started")
    # expected_pred_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, lamb=lamb)
    
    print("exp_l1 alg started")
    expected_l1_est_beta <- find_f1_coefs_expected_rank_l1_algorithm(Y, X, lamb=lamb)
    
    # fixed_est_betas[ , j] <- fixed_pred_beta
    # expected_est_betas[ , j] <- expected_pred_beta
    expected_l1_est_betas[ , j] <- expected_l1_est_beta
  }
  return (list("fixed_est_betas"=fixed_est_betas, 
               "expected_est_betas"=expected_est_betas,
               "expected_l1_est_betas"=expected_l1_est_betas,
               "betas"=betas, "n"=n, "m"=m, "lamb"=lamb))
}

parse_result_and_save <- function(results) {
  betas <- results[[1]]$betas
  n <- results[[1]]$n
  m <- results[[1]]$m
  lamb <- results[[1]]$lamb
  # fixed_betas <- matrix(0, length(results), length(betas))
  # exp_betas <- matrix(0, length(results), length(betas))
  exp_l1_betas <- matrix(0, length(results), length(betas))
  
  for(i in 1:length(results)) {
    # fixed_betas[i, ] <- results[[i]]$fixed_est_betas
    # exp_betas[i, ] <- results[[i]]$expected_est_betas
    exp_l1_betas[i, ] <- results[[i]]$expected_l1_est_betas
  }
  
  res <- list("betas"=betas, "exp_l1_betas"=exp_l1_betas, #"fixed_betas"=fixed_betas, "exp_betas"=exp_betas,
              "num_datasets"=nrow(exp_l1_betas), "num_betas"=length(betas), 
              "n"=n, "m"=m, "lamb"=lamb)
  json_data <- toJSON(res)
  
  file_name <- "../res/all_betas_l1_lamb_"
  file_name <- paste(file_name, lamb, sep='')
  file_name <- paste(file_name, "_", sep='')
  file_name <- paste(file_name, n, sep='')
  file_name <- paste(file_name, "_", sep='')
  file_name <- paste(file_name, m, sep='')
  file_name <- paste(file_name, "_", sep='')
  file_name <- paste(file_name, nrow(exp_l1_betas), sep='')
  file_name <- paste(file_name, "_", sep='')
  file_name <- paste(file_name, length(betas), sep='')
  
  write(json_data, file_name)
  
  return(res)
}

dummy_fun <- function(i) {
  res <- run_rank_regression_algorithms(n=1000, m=1, max_iter=100, batch_size=64, lamb=10)
  return(res)
}

numCores <- detectCores() - 1
print(paste("num coress  --- ", numCores))

system.time(
  results <- mclapply(c(seq(1, 100)), dummy_fun, mc.cores = numCores)
)

results
print("finished algorithms")


parse_result_and_save(results)

source("rank_reg_bivariate_pnl.R")
source("utils.R")

run_bivariate_for_sim_data <- function(rank_alg="expected_l2_rank", lamb=10, 
                                       sample_size=100) {
  data <- simulate_bivariate_pnl(sample_size)
  res <- find_bivariate_direction(data$data, rank_alg, lamb=lamb, file_name = "")
  if(res$est_direction == "1 -> 2")
  return(T)
}

numCores <- detectCores() - 1
print(paste("num coress  --- ", numCores))

dummy_fun <- function(i) {
  return(run_bivariate_for_sim_data())
}

system.time(
  results <- mclapply(1:4, dummy_fun, mc.cores = numCores)
)

results <- unlist(results)
acc <- sum(results)/length(results)
print(paste("finished algorithms, accuracy is ", acc))

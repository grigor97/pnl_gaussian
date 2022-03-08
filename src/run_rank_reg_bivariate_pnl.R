source("rank_reg_bivariate_pnl.R")

run_bivariate_for_sim_data <- function(num_datasets=100, rank_alg="expected_rank") {
  results <- c()
  for(i in 1:num_datasets) {
    print(paste("iter ", i))
    data <- simulate_bivariate_pnl(n=100)
    res <- find_bivariate_direction(data$data, rank_alg)
    results <- c(results, res$est_direction)
  }
  
  return(table(results))
}

run_bivariate_for_sim_data()

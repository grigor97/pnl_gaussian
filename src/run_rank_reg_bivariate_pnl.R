source("rank_reg_bivariate_pnl.R")
source("utils.R")
library(rjson)

run_bivariate_for_sim_data <- function(rank_alg, lamb, 
                                       sample_size) {
  data <- simulate_bivariate_pnl(sample_size)
  res <- find_bivariate_direction(data$data, rank_alg, lamb=lamb, file_name = "")
  return(res$est_direction == "1 -> 2")
}

algs <- c("expected_l2_rank")#, "expected_l1_rank", "fixed_point")
lambdas <- c(0.1, 3, 10, 30, 100)
sample_sizes <- c(1000)

all_accs <- c()



num_datasets <- 100

numCores <- detectCores() - 1
print(paste("num coress  --- ", numCores))

for(n in sample_sizes) {
  for(alg in algs) {
    for(lamb in lambdas) {
      dummy_fun <- function(i) {
        return(run_bivariate_for_sim_data(rank_alg=alg, lamb=lamb, 
                                          sample_size=n))
      }
      
      system.time(
        results <- mclapply(1:num_datasets, dummy_fun, mc.cores = numCores)
      )
      
      results <- unlist(results)
      acc <- sum(results)/length(results)
      print(paste("for alg ", alg))
      print(paste("and lamb ", lamb))
      print(paste("finished, accuracy is ", acc))
      all_accs <- c(all_accs, acc)
    }
  }
}

res <- list("lambdas"=lambdas, "all_accs"=all_accs, "algs"=algs)

res

json_data <- toJSON(res)
write(json_data, "../res/bivariate_pnl/results_exp_l2.json")


# ress <- fromJSON(file = "../res/bivariate_pnl/results.json")
# ress

source("bivariate_pnl/ltm_bivariate_pnl.R")
source("utils.R")
library(rjson)
library(parallel)

run_bivariate_for_sim_data <- function(beta_alg, f2_inv_alg, sample_size) {
  data <- simulate_bivariate_pnl_idf(sample_size)
  res <- find_bivariate_direction(data$data, beta_alg, f2_inv_alg)
  return(res$est_direction == "1 -> 2")
}

# run_bivariate_for_sim_data("prl", "rank_reg", 100)

est_beta_alg <- "prl"
est_f2_inv_alg <- "rank_reg"
sample_sizes <- c(1000)

all_accs <- c()

num_datasets <- 100

numCores <- detectCores() - 1
print(paste("num coress  --- ", numCores))

for(n in sample_sizes) {
  dummy_fun <- function(i) {
    return(run_bivariate_for_sim_data(beta_alg=est_beta_alg, 
                                      f2_inv_alg=est_f2_inv_alg, sample_size=n))
  }

  system.time(
    results <- mclapply(1:num_datasets, dummy_fun, mc.cores = numCores)
  )

  results <- unlist(results)
  acc <- sum(results)/length(results)
  print(paste("finished, accuracy is ", acc))
  all_accs <- c(all_accs, acc)
}

res <- list("all_accs"=all_accs, "est_beta_alg"=est_beta_alg, 
            "est_f2_inv_alg"=est_f2_inv_alg)

res

save_file <- "../res/bivariate_pnl/results_idf"
save_file <- paste(save_file, est_beta_alg, sep="_")
save_file <- paste(save_file, est_f2_inv_alg, sep="_")
save_file <- paste(save_file, num_datasets, sep="_")
save_file <- paste(save_file, sample_sizes[1], sep="_")
save_file <- paste(save_file, ".json", sep="")

json_data <- toJSON(res)
write(json_data, save_file)


# ress <- fromJSON(file = "../res/bivariate_pnl/results.json")
# ress

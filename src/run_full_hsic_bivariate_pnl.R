source("multivariate_pnl/full_hsic_pnl.R")
source("utils.R")
library(rjson)
library(parallel)

run_bivariate_hsic <- function(sample_size) {
  data <- simulate_bivariate_pnl_idf(n=sample_size)
  
  df <- data.frame(data$data)
  gt_beta <- data$beta
  gt_beta
  
  head(df)
  
  res <- order_recovery_by_last_node(df)
  print(res)
  return(res[1] == "X1")
}

# run_bivariate_hsic(100)

alg <- "full_hsic" 
sample_sizes <- c(500, 1000, 3000)

all_accs <- c()



num_datasets <- 100

numCores <- detectCores() - 1
print(paste("num coress  --- ", numCores))

for(n in sample_sizes) {
  dummy_fun <- function(i) {
    return(run_bivariate_hsic(sample_size=n))
  }
  
  system.time(
    results <- mclapply(1:num_datasets, dummy_fun, mc.cores = numCores)
  )
  
  results <- unlist(results)
  acc <- sum(results)/length(results)
  print(paste("for alg ", alg))
  print(paste("and sample size ", n))
  print(paste("finished, accuracy is ", acc))
  all_accs <- c(all_accs, acc)
}

res <- list("sample_sizes"=sample_sizes, "all_accs"=all_accs, "alg"=alg)

res

save_file <- "../res/bivariate_pnl/results_idf"
save_file <- paste(save_file, alg, sep="_")
save_file <- paste(save_file, num_datasets, sep="_")
save_file <- paste(save_file, ".json", sep="")

json_data <- toJSON(res)
write(json_data, save_file)


# ress <- fromJSON(file = "../res/bivariate_pnl/results.json")
# ress


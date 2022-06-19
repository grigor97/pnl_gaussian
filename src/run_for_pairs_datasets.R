source("bivariate_pnl/ltm_biv_pnl.R")

library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores-1)
registerDoParallel(cl) 

beta_alg = "prl"
f2_inv_alg = "rank_reg"

load_pairs_run_bivariate <- function(i, files) {
  f_name <- files[i]
  print(paste("alg starts for ", f_name))
  full_path <- paste("../data/pairs/", f_name, sep = '')
  
  
  dat <- scan(full_path)
  dat <- matrix(dat, ncol=2, byrow = T)
  
  print("number of rows in a dataset is ")
  print(nrow(dat))
  
  if(nrow(dat) > 1000) {
    random_inds <- sample(1:nrow(dat), 1000)
    dat <- dat[random_inds, ]
  }
  
  res <- find_bivariate_direction(dat, beta_alg, f2_inv_alg)
  
  return(c(f_name, res$est_direction))
}

# run for pairs
files <- list.files(path="../data/pairs/", pattern = "*\\d.txt")
files

# print(files[-c(52, 53, 54, 55, 71, 81, 82, 83, 105)])
files <- files[-c(52, 53, 54, 55, 71, 81, 82, 83, 105)]
files
# res <- load_pairs_run_bivariate(1, files)

res <- foreach(i=1:length(files), .combine="rbind", .packages = c("EnvStats", "dHSIC")) %dopar% {
  .GlobalEnv$cdf_z <- cdf_z
  load_pairs_run_bivariate(i, files)
}

res

stopCluster(cl)




# gt_dir <- readLines("../data/pairs/answers.txt")
# gt_dir <- gt_dir[1:52]
# 
# sum(gt_dir == df$direction)
# 
# write.csv(df, file="../data/predictions_on_pairs_exp_rank_1_52.csv", row.names = F)
# 
# d <- read.csv("../data/predictions_on_pairs_exp_rank_1_52.csv", header = T)



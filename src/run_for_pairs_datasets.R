source("rank_reg_bivariate_pnl.R")

# library(foreach)
# library(doParallel)

# myCluster <- makeCluster(6, type = "PSOCK")
# registerDoParallel(myCluster)

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
  
  res <- find_bivariate_direction(dat, rank_alg="expected_rank", f_name)
  
  print(res$est_direction)
  print(res$res1$dhsic_val)
  print(res$res2$dhsic_val)
  
  return(res)
}

# run for pairs
files <- list.files(path="../data/pairs/", pattern = "*\\d.txt")
files

# load_pairs_run_bivariate(1, files)
# 
# foreach_res <- foreach(i=1:length(files)) %dopar% 
# {
#   load_pairs_run_bivariate(i, files)
# }
# foreach_res
# stopCluster(myCluster)

print(files[-c(81:84)])
files <- files[-c(81:84)]
files[1]
df <- data.frame(f_name=character(), direction=character())
for(i in 1:length(files)) {
  res <- load_pairs_run_bivariate(i, files)
  df[nrow(df)+1, ] <- c(res$f_name, res$est_direction)
}

df

gt_dir <- readLines("../data/pairs/answers.txt")
gt_dir <- gt_dir[1:52]

sum(gt_dir == df$direction)

# df <- data.frame(f_name=character(), direction=character())
# df[nrow(df)+1, ] <- c(f_name, res$est_direction)
# print(df1)

# df[df$direction == "1 -> 2", ]
# df[df$direction == "2 -> 1", ]

write.csv(df, file="../data/predictions_on_pairs_exp_rank_1_52.csv", row.names = F)

d <- read.csv("../data/predictions_on_pairs_exp_rank_1_52.csv", header = T)

sum(d$direction == d$gt)

library(doParallel)
library(rjson)
source("utils.R")
source("linear_transofrmation_models/noise_estimation_ltm.R")

no_cores <- detectCores()
cl <- makeCluster(no_cores-1)
registerDoParallel(cl) 

n <- 1000
m <- 1
num_datasets <- 100
save_file_name <- "../res/rank_regression/"
alg_name <- "prl.rankreg"
alg <- function(Y, X) {
  res = noise.est(Y, X, "prl", "rank_reg")
  return(res)
}

run.noise.est.alg <- function(n, m) {
  data <- simulate_rank_regression_data(n, m, beta_max = 100)
  res <- alg(data$Y, data$X)
  
  return(matrix(c(data$Y, res$est_f2_inv), nrow=length(data$Z)))
}

res <- foreach(i=1:num_datasets, .combine = 'rbind') %dopar% {
  .GlobalEnv$cdf_z <- cdf_z
  run.noise.est.alg(n, m)
}

stopCluster(cl)

res

library(MASS)
write.csv(data.frame(res), file="../res/h_esitmation_rank_reg_1000_100.csv", row.names = F)
head(data.frame(res))
df <- read.csv("../res/h_esitmation_rank_reg_1000_100.csv")
head(df)
dim(df)

library(ggplot2)

f <- function(x) (x - 4.7)^3
xs <- seq(min(res[, 1])-0.5, max(res[, 1])+0.5, 0.1)
ggplot() + geom_point(aes(x = res[, 1], y=res[, 2], colour="estimated h")) + 
  geom_line(aes(x = xs, y=sapply(xs, f), colour="true h")) +
  scale_color_manual(name='', breaks = c("estimated h", "true h"), values = c("black", "red")) +
  labs(x='y', y='h(y)') +
  theme(legend.position = 'top')





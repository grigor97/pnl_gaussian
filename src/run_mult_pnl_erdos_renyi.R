source("utils.R")
source("multivariate_pnl/ltm_mult_pnl.R")

library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores-1)
registerDoParallel(cl) 

n = 100
d = 3
num_datasets <- 4

beta_alg <- "prl"
f2_inv_alg <- "rank_reg"

res <- foreach(i=1:num_datasets, .combine="rbind", .packages = c("dHSIC")) %dopar% {
  .GlobalEnv$cdf_z <- cdf_z
  data <- simulate.mult.pnl.erdos.renyi(n, d)
  A <- data$A
  df <- data.frame(data$X)
  ord <- order_recovery_ltm(df, beta_alg, f2_inv_alg)
  for(i in 1:(d-1)) {
    for(j in (i+1):d) {
      v1 <- strtoi(substr(ord[i], 2, 100))
      v2 <- strtoi(substr(ord[j], 2, 100))
      if(A[v2, v1]) {
        return(F)
      }
    }
  }
  return(T)
}

res
print("true orders are")
print(sum(res))

stopCluster(cl)

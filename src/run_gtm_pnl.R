source("utils.R")
source("multivariate_pnl/gtm_mult_pnl.R")

library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores-1)
registerDoParallel(cl) 

n = 30
num_datasets <- 4

res <- foreach(i=1:num_datasets, .combine="rbind", .packages = c("dHSIC")) %dopar% {
  X <- simulate.bv.pnl.gtm(n)
  df <- data.frame(X)
  ord <- order_recovery_gtm(df)
  if(strtoi(substr(ord[1], 2, 100))==1) {
    return(T)
  } 
  return(F)
}

res
print("true orders are")
print(sum(res))

stopCluster(cl)

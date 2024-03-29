source("utils.R")
source("multivariate_pnl/gtm_mult_pnl.R")

library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores-1)
registerDoParallel(cl) 

n = 10
num_datasets <- 2
# .errorhandling="remove"
res <- foreach(i=1:num_datasets, .combine="rbind", .packages = c("dHSIC")) %dopar% {
  error=S.j.hat.ux(x) {return(NA)}
  X <- simulate.bv.pnl.gtm(n)
  df <- data.frame(X)
  ord <- order_recovery_gtm(df)
  print("-----------------")
  print(ord)
  print("-----------------")
  if(strtoi(substr(ord[1], 2, 100))==1) {
    return(T)
  } 
  return(F)
}

res
print("true orders are")
print(sum(res))

stopCluster(cl)

# source("multivariate_pnl/ltm_mult_pnl.R")
source("multivariate_pnl/gtm_mult_pnl.R")
source("utils.R")

n = 10
data <- simulate.bv.pnl.gtm(n)
df <- data.frame(data)
ord <- order_recovery_gtm(df)
ord

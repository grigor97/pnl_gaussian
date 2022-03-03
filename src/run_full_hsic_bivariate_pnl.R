source("full_hsic_pnl.R")
source("rank_reg_bivariate_pnl.R")

correct_num <- 0
for(i in 1:100) {
  data <- simulate_bivariate_pnl(n=100)
  
  df <- data.frame(data$data)
  gt_beta <- data$beta
  gt_beta
  
  head(df)
  
  res <- order_recovery_by_last_node(df)
  print(res)
  if(res[1] == "X1") {
    correct_num <- correct_num + 1
  }
}

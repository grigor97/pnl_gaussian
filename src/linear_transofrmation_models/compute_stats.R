library(rjson)
source("utils.R")

res <- fromJSON(file = "../res/rank_regression/1000/all_betas_mc_0_1000_1_100_13.json")

gt_betas <- res$betas
est_betas <- est_betas <- matrix(res$est_betas, res$num_datasets, res$num_betas)

gt_betas
dim(est_betas)

comute.RB.Var.MSE(est_betas, gt_betas)


mean_ests <- matrix(colMeans(est_betas), nrow = nrow(est_betas), ncol = ncol(est_betas), byrow = T)
colMeans((est_betas - mean_ests)^2)
colMeans(est_betas)


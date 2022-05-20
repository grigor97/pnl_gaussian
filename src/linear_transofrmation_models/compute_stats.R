library(rjson)
source("utils.R")

res <- fromJSON(file = "../res/rank_regression/all_betas_prl_1000_1_100_13.json")
res$est_betas
gt_betas <- res$betas
est_betas <- matrix(res$est_betas, res$num_datasets, res$num_betas)

gt_betas
dim(est_betas)

comute.RB.Var.MSE(est_betas, gt_betas)



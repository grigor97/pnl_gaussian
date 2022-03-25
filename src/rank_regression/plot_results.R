# for local run
library(ggplot2)
library(rjson)

save_plots <- function(est_betas1, est_betas2, name1, name2, gt_betas, alg_name, 
                       n, m, lamb, from=1, to=100, file_n='../plots/rank_regression/1000/'){
  if(to > length(gt_betas)) {
    to <- length(gt_betas)
  }
  estimated_betas_cut1 <- est_betas1[, from:to]
  estimated_betas_cut2 <- est_betas2[, from:to]
  gt_betas <- gt_betas[from:to]
  number_of_datasets <- nrow(estimated_betas_cut1)
  stacked_vals1 <- stack(as.data.frame(estimated_betas_cut1))
  stacked_vals2 <- stack(as.data.frame(estimated_betas_cut2))
  
  # title_name <- paste(alg_name, 'rank regression for')
  # title_name <- paste(title_name, number_of_datasets)
  # title_name <- paste(title_name, "datasets")
  pl <- ggplot() + 
    geom_boxplot(aes(x=stacked_vals1$ind, y=stacked_vals1$values, colour=name1)) +
    geom_boxplot(aes(x=stacked_vals2$ind, y=stacked_vals2$values, colour=name2), alpha=0.3) +
    geom_point(aes(x=unique(stacked_vals1$ind), y=gt_betas, colour='ground truth')) +
    labs(x="",y="betas") + # labs(title=title_name, x="",y="betas") +
    scale_color_manual(name='',
                       breaks=c(name1, name2, 'ground truth'),
                       values=c('blue4', 'black', 'red')) +
    guides(colour = guide_legend(override.aes = list(
      linetype = c("solid", "solid", "blank"),
      color = c("blue4", 'black', "red")
    ))) +
    # theme(legend.position=c(0.15,0.91), plot.title = element_text(hjust = 0.5))
    theme(legend.position='top', plot.title = element_text(hjust = 0.5))
  
  file_name <- paste(file_n, alg_name, sep='')
  file_name <- paste(file_name, "lamb_", sep='')
  file_name <- paste(file_name, lamb, sep='')
  file_name <- paste(file_name, "_", sep='')
  file_name <- paste(file_name, n, sep='')
  file_name <- paste(file_name, "_", sep='')
  file_name <- paste(file_name, m, sep='')
  file_name <- paste(file_name, "_", sep='')
  file_name <- paste(file_name, number_of_datasets, sep='')
  file_name <- paste(file_name, "_", sep='')
  file_name <- paste(file_name, from, sep='')
  file_name <- paste(file_name, "_", sep='')
  file_name <- paste(file_name, to, sep='')
  
  file_name <- paste(file_name, '.png', sep='')
  ggsave(filename = file_name, plot = pl)
  
  return(pl)
}

res10 <- fromJSON(file = "../res/rank_regression/1000/all_betas_1000_1_no_regularization")
res10

exp_l2_betas <- matrix(res10$exp_betas, res10$num_datasets, res10$num_betas)
fixed_betas <- matrix(res10$fixed_betas, res10$num_datasets, res10$num_betas)
# exp_l1_betas <- matrix(res10$exp_l1_betas, res10$num_datasets, res10$num_betas)
betas <- res10$betas

save_plots(exp_l2_betas, fixed_betas, "exp no-reg", "fixed_point", betas, alg_name='', 
           res10$lamb, from=1, to = 100, n=res10$n, m=res10$m)


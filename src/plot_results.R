# for local run
library(ggplot2)
library(rjson)

ress <- fromJSON(file = "/Users/grigorkeropyan/pnl_gaussian/res/all_betas_l2_lamb1_1000_1_100_13")
ress

exp_betas <- matrix(ress$exp_betas, ress$num_datasets, ress$num_betas)
fixed_betas <- matrix(ress$fixed_betas, ress$num_datasets, ress$num_betas)
betas <- ress$betas

save_plots <- function(estimated_betas, gt_betas, alg_name, n, m, from=1, to=100,
                       file_n='/Users/grigorkeropyan/pnl_gaussian/plots/'){
  if(to > length(gt_betas)) {
    to <- length(gt_betas)
  }
  estimated_betas_cut <- estimated_betas[, from:to]
  gt_betas <- gt_betas[from:to]
  number_of_datasets <- nrow(estimated_betas_cut)
  stacked_vals <- stack(as.data.frame(estimated_betas_cut))
  
  title_name <- paste(alg_name, 'rank regression for')
  title_name <- paste(title_name, number_of_datasets)
  title_name <- paste(title_name, "datasets")
  pl <- ggplot() + geom_boxplot(aes(x=stacked_vals$ind, y=stacked_vals$values, colour='estimated betas')) +
    geom_point(aes(x=unique(stacked_vals$ind), y=gt_betas, colour='ground truth betas')) +
    labs(title=title_name, x="",y="betas") +
    scale_color_manual(name='',
                       breaks=c('estimated betas', 'ground truth betas'),
                       values=c('black', 'red')) +
    guides(colour = guide_legend(override.aes = list(
      linetype = c("solid", "blank"),
      color = c("black","red")
    ))) +
    # theme(legend.position=c(0.15,0.91), plot.title = element_text(hjust = 0.5))
    theme(legend.position='top', plot.title = element_text(hjust = 0.5))
  file_name <- paste(file_n, alg_name, sep='')
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

save_plots(exp_betas, betas, alg_name='exp_rank_lamb_1888_', from=1, to = 7, 
           n=ress$n, m=ress$m)
save_plots(fixed_betas, betas, alg_name='fix_point_lamb_1888_', from=1, to = 7,
           n=ress$n, m=ress$m)

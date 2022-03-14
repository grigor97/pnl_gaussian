library(rjson)
library(ggplot2)
library(latex2exp)

res <- fromJSON(file = "../res/bivariate_pnl/results_idf_expected_l2_rank_100_1000.json")

lambdas <- res$lambdas
accuracy <- res$all_accs
accuracy
pl <- ggplot() + geom_point(aes(x=lambdas, y = accuracy), colour="red") +
  labs(title = TeX(r'(expected rank ell2 algorithm, model $Y = (beta X^2 + N)^{1/3}$)')) +
  theme(plot.title = element_text(hjust = 0.5))

pl

ggsave(filename = "../plots/bivariate_pnl/results_exp_l2_idf.png", plot = pl)

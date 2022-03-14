library(rjson)
library(ggplot2)
library(latex2exp)

res_l2 <- fromJSON(file = "../res/bivariate_pnl/results_exp_l2.json")
res_l1 <- fromJSON(file = "../res/bivariate_pnl/results_exp_l1.json") 

res_fixed <- fromJSON(file = "../res/bivariate_pnl/results_fixed_point.json")
res_fixed

lambdas_l2 <- res_l2$lambdas
accuracy_l2 <- res_l2$all_accs
accuracy_l2

lambdas_l1 <- res_l1$lambdas
accuracy_l1 <- res_l1$all_accs
accuracy_l1

pl <- ggplot() + geom_line(aes(x=lambdas_l2, y = accuracy_l2, colour="ell2")) +
  geom_line(aes(x=lambdas_l1, y = accuracy_l1, colour="ell1")) +
  geom_line(aes(x=lambdas_l1, y = rep(res_fixed$all_accs[2], length(lambdas_l1)), colour="fixed point")) +
  labs(title = TeX(r'(model $X_2 = (beta X_1 + N)^{1/3}$)'), x=TeX(r'(\lambda)'), y="accuracy") +
  scale_color_manual(name='',
                     breaks=c('ell2', 'ell1', "fixed point"),
                     values=c('blue', 'red', "black")) +
  theme(legend.position='top', plot.title = element_text(hjust = 0.5))

pl

ggsave(filename = "../plots/bivariate_pnl/results_exp_non_idf.png", plot = pl)

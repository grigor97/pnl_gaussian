source('rank_regression/rank_regression.R')
source("utils.R")
library(ggplot2)

library(stochQN)


simu_rr <- function(n, m) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  
  noise <- rnorm(n)
  X <- matrix(rnorm(n*m), n, m)
  
  beta <- runif(n=m, min=-10000, max=10000)
  # beta <- runif(n=m, min=-1, max=1)
  # beta[1] <- 0
  # beta[3] <- 0
  # beta[4] <- 0
  
  Y <- X %*% beta + noise
  Y <- exponent(Y, 1/3) + 4.7
  
  res <- list("X"=X, "Y"=Y, "beta"=beta)
  return(res)
}
data <- simu_rr(500, 1)
gt_beta <- data$beta
X <- data$X
Y <- data$Y

prl_res <- rank.reg.prl.gaussian(Y, X)

prl_res


g <- function(x) {
  res <- (-x*dnorm(x)*pnorm(x) - dnorm(x)^2)/(2*pnorm(x)^2)
  return(res)
}
g(-100:100)
n <- 30
ggplot() + geom_line(aes(x=0:n, y=g(0:n)))

# Monte Carlo methods
# data <- simulate_rank_regression_data(1000, 1)
# gt_beta <- data$beta
# X <- data$X
# Y <- data$Y
# 
# gt_beta
# 
# est_beta <- find_f1_coefs_cond_monte_carlo_algorithm(Y, X, M=1000, max_iter = 100)
# est_beta
# gt_beta

G_j_beta <- function(beta, j, X) {
  val <- 0
  for(i in 1:nrow(X)) {
    if(i == j) {
      next
    }
    val <- val + pnorm(sum((X[j, ] - X[i, ])*beta)/2**0.5)
  }
  return (1 + val)
}

S_beta <- function(beta, X, ranks_Y, lamb=0, penalty='e') {
  val <- 0
  for(j in 1:nrow(X)) {
    val <- val + (ranks_Y[j] - G_j_beta(beta, j, X))**2
  }
  
  if (penalty=="ell1") {
    val <- val + lamb*sum(abs(beta))
  } else if(penalty=="ell2") {
    val <- val + lamb*sum((beta)**2)
  } 
  
  return (val)
}

data <- simulate_rank_regression_data(1000, 1)
gt_beta <- data$beta
X <- data$X
Y <- data$Y

est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X)

est_beta_0 <- find_f1_coefs_expected_rank_algorithm(Y, X, lamb=0, penalty = 'e')

S_beta(est_beta, X, rank(Y), 0, 'e')
S_beta(est_beta_0, X, rank(Y), 0, 'e')
S_beta(gt_beta, X, rank(Y), 0, 'e')
est_beta
est_beta_0
gt_beta

S_beta(est_beta, X, rank(Y), 10, 'ell2')
S_beta(est_beta_0, X, rank(Y), 10, 'ell2')
S_beta(gt_beta, X, rank(Y), 10, 'ell2')

R <- rank(Y)

S_beta(gt_beta-4, X, R)
S_beta(gt_beta, X, R)

dumm <- function(beta) {
  return(S_beta(beta, X, R))
}

gt_beta
val_betas <- seq(9990, 10000, 1)
S_vals <- sapply(val_betas, dumm)
S_vals
df <- data.frame(val_betas=val_betas, S_vals=S_vals)

ggplot(df, aes(x= val_betas, y=S_vals)) + geom_line(color='red')

ggplot(df, aes(x= val_betas, y=S_vals)) + 
  geom_line(color='red') + xlim(c(-2, 3)) #+ 
  #geom_point(color='blue', size=0.5, alpha=0.1) #+
  #scale_x_continuous(breaks = val_betas)

ggplot(df, aes(x= val_betas, y=S_vals)) + geom_line(color='red') + xlim(c(-1000, 0)) + ylim(c(50, 125))

ggplot(df, aes(x= val_betas, y=S_vals)) + geom_line(color='red') + xlim(c(0, 50)) + ylim(c(80000, 333000))


gt_beta
min_beta <- val_betas[which.min(S_vals)]
min_beta
min(S_vals)


num <- 1000
smp <- rnorm(num)
st <- sort(smp)
min(st[2:num] - st[1:(num-1)])

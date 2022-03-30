source('rank_regression/rank_regression.R')
source("utils.R")
library(ggplot2)


# Monte Carlo methods
data <- simulate_rank_regression_data(1000, 1)
gt_beta <- data$beta
X <- data$X
Y <- data$Y

gt_beta

est_beta <- find_f1_coefs_cond_monte_carlo_algorithm(Y, X, M=1000, max_iter = 100)
est_beta
gt_beta

expected_rank_algorithm <- function(Y, X, lamb=10, penalty="ell2") {
  print(paste("starting expected rank ", penalty))
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
  
  grad_G_j_beta <- function(beta, j, X) {
    grad <- 0
    for(i in 1:nrow(X)) {
      if(i == j) {
        next
      }
      grad <- grad + dnorm(sum((X[j, ] - X[i, ])*beta)/2**0.5)*(X[j, ] - X[i, ])/2**0.5
    }
    return (grad)
  }
  
  S_beta <- function(beta, X, ranks_Y, lamb, penalty) {
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
  
  grad_S_beta <- function(beta, X, ranks_Y, lamb, penalty) {
    grad <- 0
    for(j in 1:nrow(X)) {
      grad <- grad - 2*(ranks_Y[j] - G_j_beta(beta, j, X))*grad_G_j_beta(beta, j, X)
    }
    
    if (penalty=="ell1") {
      grad <- grad + lamb*sign(beta)
    } else if(penalty=="ell2") {
      grad <- grad + 2*lamb*beta
    } 
    
    return (grad)
  }
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  
  m <- ncol(X)
  n <- nrow(X)
  coefs <- matrix(runif(m, min=-10, max=10), m, 1)
  ranks_Y <- rank(Y)
  
  est_beta <- optim(par=coefs, fn=S_beta, gr=grad_S_beta, method = "BFGS", 
                    control = list(trace=T, maxit=4, REPORT=1),
                    X=X, ranks_Y=ranks_Y, lamb=lamb, penalty=penalty)
  
  return(est_beta)
}

data <- simulate_rank_regression_data(1000, 1)
gt_beta <- data$beta
X <- data$X
Y <- data$Y

gt_beta

system.time(
  res_noreg <- expected_rank_algorithm(Y, X, lamb = 10, penalty = "dd")
)

res_noreg

system.time(
  res_reg <- expected_rank_algorithm(Y, X, lamb = 10, penalty = "ell2")
)

res_reg


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

S_beta <- function(beta, X, ranks_Y, lamb=0, penalty="eeeee") {
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

data <- simulate_rank_regression_data(100, 1)
gt_beta <- data$beta
X <- data$X
Y <- data$Y

gt_beta
R <- rank(data$Y)

S_beta(gt_beta+6.5, X, R)
S_beta(gt_beta, X, R)

dumm <- function(beta) {
  return(S_beta(beta, X, R))
}

gt_beta
val_betas <- seq(-1000, 1000, 1)
S_vals <- sapply(val_betas, dumm)
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




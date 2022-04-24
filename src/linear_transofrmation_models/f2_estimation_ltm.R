source("linear_transofrmation_models/beta_estimation_ltm.R")

# Works for beta around 10 since otherwise log \Phi becomes NaN
# name is rl
f2.inv.est.rl <- function(Y, X, est_beta) {
  m <- ncol(X)
  n <- nrow(X)
  X <- X - matrix(rep(colMeans(X), n), n, m, byrow = T)

  rl_y <- function(t, y, Y, X, est_beta) {
    m_X <- X %*% est_beta
    lik <- 0
    for(j in 1:n) {
      lik <- lik + 1*(Y[j] <= y)*log(pnorm(t - m_X[j])) + 
        1*(Y[j] > y)*log(pnorm(-t + m_X[j])) 
    }
    return(-lik)
  }
  
  grad_rl_y <- function(t, y, Y, X, est_beta) {
    m_X <- X %*% est_beta
    grad <- 0
    for(j in 1:n) {
      grad <- grad + 1*(Y[j] <= y)*dnorm(t - m_X[j])/pnorm(t - m_X[j]) + 
        1*(Y[j] > y)*dnorm(-t + m_X[j])/pnorm(-t + m_X[j])
    }
    return(-grad)
  }
  
  f2_inv_est <- c()
  
  for(y in Y) {
    res <- optim(0, fn=rl_y, gr=grad_rl_y,  method = "BFGS", 
                 control = list(trace=T, REPORT=1),
                 y=y, Y=Y, X=X, est_beta=est_beta)
    f2_inv_est <- c(f2_inv_est, res$par)
  }
  
  
  return(f2_inv_est)
}

# works the best: name is rank_reg
f2.inv.est.rank.reg <- function(Y, X, est_beta) {
  ranks_Y <- rank(Y)
  empirical_cdf_Y <- ranks_Y/(length(Y)+1)
  m_X <- X %*% est_beta
  
  f2_inv_y <- as.vector(sapply(empirical_cdf_Y, inverse_cdf_z, arg=m_X))
  return(f2_inv_y)
}



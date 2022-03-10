set.seed(12)

simulate_bivariate_pnl <- function(n) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  X <- matrix(rnorm(n), n, 1)
  beta <- runif(1, -100, 100)
  noise <- rnorm(n)
  Y <- X %*% beta + noise
  Y <- exponent(Y, 1/3) + 4.7
  data <- cbind(X, Y)
  res <- list("data"=data, "beta"=beta)
  
  return(res)
}

simulate_rank_regression_data <- function(n, m) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  
  noise <- rnorm(n)
  X <- matrix(rnorm(n*m), n, m)
  
  beta <- runif(n=m, min=-100, max=100)
  # beta[1] <- 0
  # beta[3] <- 0
  # beta[4] <- 0
  
  Y <- X %*% beta + noise
  Y <- exponent(Y, 1/3) + 4.7
  
  res <- list("X"=X, "Y"=Y, "beta"=beta)
  return(res)
}
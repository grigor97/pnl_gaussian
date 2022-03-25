source("utils.R")

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
                    control = list(trace=T),
                    X=X, ranks_Y=ranks_Y, lamb=lamb, penalty=penalty)
  
  return(est_beta$par)
}

data <- simulate_rank_regression_data(200, 3)
gt_beta <- data$beta
X <- data$X
Y <- data$Y

gt_beta

system.time(
  res <- expected_rank_algorithm(Y, X, lamb = 10, penalty = "dd")
)

res


cdf_z <- function(y, subtract_values) {
  res = 0
  for(i in subtract_values) {
    res = res + pnorm(y - i)
  }
  res = res / length(subtract_values)
  return(res)
}

inverse <- function(f, lower, upper){
  function(y, arg){
    # if(f(lower, arg) - y > 0) {
    #   return(lower)
    # }
    # if (f(upper, arg) - y < 0) {
    #   return(upper)
    # }
    uniroot(function(x){f(x, arg) - y}, lower = lower, upper = upper, 
            extendInt="upX", tol=1e-3,)[1]$root
  }
}

inverse_cdf_z <- inverse(cdf_z, -100, 100)

find_root <- function(val, sub_vals) {
  dummy <- function(x, val, sub_vals) {
    return((cdf_z(x, sub_vals) - val)**2)
  }
  
  grad_cdf_z <- function(y, sub_vals) {
    grad = 0
    for(i in sub_vals) {
      grad = grad + dnorm(y - i)
    }
    grad = grad / length(sub_vals)
    return(grad)
  }
  
  grad_dummy <- function(x, val, sub_vals) {
    grad <- 2*(cdf_z(x, sub_vals) - val)*grad_cdf_z(x, sub_vals)
    return (grad)
  }
  
  res <- optim(par=runif(1, -10, 10), fn=dummy, gr=grad_dummy, method = "BFGS", 
               # control = list(trace=T),
               val=val, sub_vals=sub_vals)
  return(res$par)
}


n <- 2000
vals <- runif(n, -10000, 10000) 
vals <- rank(vals)/(length(vals) + 1)
sub_vals <- runif(n, -10000, 10000) 

system.time(
  res_uniroot <- sapply(vals, inverse_cdf_z, arg=sub_vals)
)

system.time(
  res_optim <- sapply(vals, find_root, sub_vals=sub_vals)
)



zeros_optim <- sapply(res_optim, cdf_z, subtract_values=sub_vals) 
zeros_uniroot <- sapply(res_uniroot, cdf_z, subtract_values=sub_vals)
sum(zeros_optim < 0.5)
sum(zeros_optim < zeros_uniroot)

sum((zeros_optim - zeros_uniroot)**2)**0.5

zeros_optim
zeros_uniroot



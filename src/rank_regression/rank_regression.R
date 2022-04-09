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

find_f1_coefs_fixed_point_stochastic <- function(Y, X, batch_size=64, 
                                                 max_iter=100, tol=1e-9) {
  print("starting fixed point algorithm")
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  
  m <- ncol(X)
  n <- nrow(X)
  coefs <- matrix(runif(m, min=-10, max=10), m, 1)
  # coefs <- beta
  
  err <- 1
  for(iter in 1:max_iter) {
    random_inds <- sample(1:n, batch_size)
    batch_X <- X[random_inds, ]
    batch_Y <- Y[random_inds, ]
    
    A <- t(batch_X) %*% batch_X
    if(nrow(A) > 1) {
      P <- solve(A) %*% t(batch_X)
    } else {
      P <- 1/A[1, 1] * t(batch_X)
    }
    ranks_Y <- rank(batch_Y)
    empirical_cdf_Y <- ranks_Y/(batch_size+1)
    
    # print(paste("iter for fixed point alg  ", iter))
    sub_vals <- batch_X %*% coefs
    yhat <- as.vector(sapply(empirical_cdf_Y, inverse_cdf_z, arg=sub_vals))
    
    a <- P %*% yhat
    err <- sum((a - coefs)**2)**0.5
    
    coefs <- a
    
    # print(paste("error is ", err))
    # print("coefs are ")
    # print(coefs)
    if(err < tol) {
      print("converged")
      break
    }
  }
  
  return(coefs)
}

# possible penalties are ell1 and ell2
find_f1_coefs_expected_rank_algorithm <- function(Y, X, lamb=10, penalty="ell2") {
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
                    control = list(trace=T, REPORT=1),
                    X=X, ranks_Y=ranks_Y, lamb=lamb, penalty=penalty)
  
  return(est_beta$par)
}

softmax_stable <- function(input_x) {
  # softmax(x) = softmax(x+c)
  x_shifted <- input_x - max(input_x)
  numerator <- exp(x_shifted)
  denominator <- sum(numerator)
  
  softmax <- numerator/denominator
  return(softmax)
}

find_f1_coefs_monte_carlo_algorithm <- function(Y, X, M=1000, max_iter=100, tol=1e-9) {
  ranks_Y <- rank(Y)
  n <- nrow(X)
  m <- ncol(X)
  # centering a matrix column wise
  X <- X - matrix(rep(colMeans(X), n), n, m, byrow = T)
  
  Z <- matrix(data=rnorm(M*n), nrow = n, ncol = M)
  Z <- apply(Z, 2, sort)
  Z <- Z[ranks_Y, ]
  
  inv_XTX <- solve(t(X) %*% X)
  B <-  inv_XTX %*% t(X) %*% Z
  Z_hat <- X %*% B
  ssf <- colSums(Z_hat**2)
  # v <- exp(ssf/2)
  
  coefs <- runif(m, min=-10, max=10)
  for(iter in 1:max_iter) {
    u_beta <- X %*% (B - coefs)
    log_u_beta <- -colSums(u_beta**2)/2
    log_w <- ssf/2 + log_u_beta
    re_weights <- softmax_stable(log_w)
    
    # print('re_weigths')
    # print(sum(re_weights == 1))
    
    next_coefs <- c(B %*% re_weights)
    if(sum((next_coefs - coefs)**2)**0.5 < tol) {
      coefs <- next_coefs
      break
    }
    coefs <- next_coefs
    # print("coefs")
    # print(coefs)
  }
  
  return (coefs)
}

find_f1_coefs_cond_monte_carlo_algorithm <- function(Y, X, M=1000, max_iter=100, tol=1e-9) {
  n <- nrow(X)
  m <- ncol(X)
  ranks_Y <- rank(Y)
  # centering a matrix column wise
  X <- X - matrix(rep(colMeans(X), n), n, m, byrow = T)
  
  inv_XTX <- solve(t(X) %*% X)
  M_mult <- inv_XTX %*% t(X)
  
  B_hat <- matrix(0, nrow=m, ncol=M)
  rss_hats <- c()
  
  for(i in 1:M) {
    # Generate standard normal sample, standardize by sample mean and standard
    # deviation, order and put into the same rank order as the original sample
    max_val <- 1000000
    while(T) {
      z <- rnorm(n)
      z_hat <- (z - mean(z))/sd(z)
      z <- sort(z_hat)
      z_hat <- z_hat[ranks_Y]
      
      # finding b hat and ssf hat
      b_hat <- M_mult %*% z_hat
      ssf_hat <- sum((X %*% b_hat)**2)
      rss_hat <- n-1-ssf_hat
      
      if(rss_hat > max_val) {
        next
      }
      
      B_hat[, i] <- b_hat
      rss_hats <- c(rss_hats, rss_hat)
      break
    }
  }
  
  t <- ((n-3)/rss_hats)**0.5
  TT <- matrix(t, m, M, byrow = T)
  log_v <- -(n-1)*log(rss_hats)/2
  
  coefs <- runif(m, min=-10, max=10)
  # print("intial coefs")
  # print(coefs)
  for(i in 1:max_iter) {
    U_beta <- X %*% (TT*B_hat - coefs)
    u_beta <- colSums(U_beta**2)
    log_u_beta <- -u_beta/2
    
    log_w <- log_v + log_u_beta
    re_weights <- softmax_stable(log_w)
  
    next_coefs <- c(B_hat %*% re_weights)
    
    if(sum((next_coefs - coefs)**2)**0.5 < tol) {
      coefs <- next_coefs
      print(paste("converged in iteration ", i))
      break
    }
    
    coefs <- next_coefs
    # print("coefs")
    # print(coefs)
  }
  
  return (coefs)
}

# https://arxiv.org/pdf/2103.13435.pdf adopted version
rank.reg.prl.gaussian <- function(Y, X, gt_beta=NA) {
  m <- ncol(X)
  n <- nrow(X)
  ord <- order(Y)
  Y <- Y[ord]
  # FIXME centering a matrix column wise 
  X <- X - matrix(rep(colMeans(X), n), n, m, byrow = T)
  
  if(m > 1) {
    X <- X[ord, ]
  } else {
    X <- matrix(X[ord], n, 1)
  }
  
  prl <- function(beta, Y, X) {
    n <- nrow(X)
    lik <- 0
    m_X <- X %*% beta
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        lik <- lik + log(pnorm((m_X[j] - m_X[i])/sqrt(2)))
      }
    }
    
    return(-lik)
  }
  
  grad_prl <- function(beta, Y, X) {
    n <- nrow(X)
    grad <- 0
    m_X <- X %*% beta
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        grad <- grad + dnorm((m_X[j] - m_X[i])/sqrt(2))*(X[j, ]- X[i, ])/(pnorm((m_X[j] - m_X[i])/sqrt(2))*sqrt(2))
      }
    }
    
    return(-grad)
  }
  
  #coefs <- runif(m, -10, 10)
  coefs <- rep(0, m)
  res <- optim(coefs, fn=prl, gr=grad_prl, method = "BFGS", 
               control = list(trace=T, REPORT=1),
               Y=Y, X=X)
  
  prl_lik <- res$value
  gt_lik <- NA
  if(!is.na(gt_beta)) {
    gt_lik <- prl(gt_beta, Y, X)
  }
  
  return(list("est_beta"=res$par, "est_obj"=prl_lik, 
              "gt_beta"=gt_beta, "gt_obj"=gt_lik))
}

rank.reg.oprl.gaussian <- function(Y, X, gt_beta=NA) {
  m <- ncol(X)
  n <- nrow(X)
  ord <- order(Y)
  Y <- Y[ord]
  # FIXME centering a matrix column wise 
  X <- X - matrix(rep(colMeans(X), n), n, m, byrow = T)
  
  if(m > 1) {
    X <- X[ord, ]
  } else {
    X <- matrix(X[ord], n, 1)
  }
  
  oprl <- function(beta, Y, X) {
    n <- nrow(X)
    lik <- 0
    m_X <- X %*% beta
    for(i in 1:(n-1)) {
      lik <- lik + log(pnorm((m_X[i+1] - m_X[i])/sqrt(2)))
    }
    
    return(-lik)
  }
  
  grad_oprl <- function(beta, Y, X) {
    n <- nrow(X)
    grad <- 0
    m_X <- X %*% beta
    for(i in 1:(n-1)) {
      grad <- grad + dnorm((m_X[i+1] - m_X[i])/sqrt(2))*(X[i+1, ]- X[i, ])/(pnorm((m_X[i+1] - m_X[i])/sqrt(2))*sqrt(2))
    }
    
    return(-grad)
  }
  
  #coefs <- runif(m, -10, 10)
  coefs <- rep(0, m)
  res <- optim(coefs, fn=oprl, gr=grad_oprl, method = "BFGS", 
               control = list(trace=T, REPORT=1),
               Y=Y, X=X)
  
  oprl_lik <- res$value
  gt_lik <- NA
  if(!is.na(gt_beta)) {
    gt_lik <- oprl(gt_beta, Y, X)
  }
  
  return(list("est_beta"=res$par, "est_obj"=oprl_lik, 
              "gt_beta"=gt_beta, "gt_obj"=gt_lik))
}

# No scale of beta, no Gaussian noise specification! Kendall rank correlation!
# Maybe change BFGS
# FIXME check only unit circle for beta: Consistency result under some assumptions
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.454.6938&rep=rep1&type=pdf
lin.tr.models.mrc <- function(Y, X, gt_beta=NA) {
  m <- ncol(X)
  n <- nrow(X)
  ord <- order(Y)
  Y <- Y[ord]
  # centering a matrix column wise 
  X <- X - matrix(rep(colMeans(X), n), n, m, byrow = T)

  if(m > 1) {
    X <- X[ord, ]
  } else {
    X <- matrix(X[ord], n, 1)
  }
  
  Sn_beta <- function(beta, Y, X) {
    n <- nrow(X)
    m_X <- X %*% beta
    Sn <- 0
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        Sn <- Sn + 1*(Y[i] < Y[j])*(m_X[i] < m_X[j])
      }
    }
    Sn <- 2*Sn/(n*(n-1))
    
    return(-Sn)
  }
  
  coefs <- rep(0, m)
  res <- optim(coefs, fn=Sn_beta, method = "BFGS", 
               control = list(trace=T, REPORT=1),
               Y=Y, X=X)
  
  mrc_Sn <- res$value
  gt_Sn <- NA
  if(!is.na(gt_beta)) {
    gt_Sn <- Sn_beta(gt_beta, Y, X)
  }
  
  return(list("est_beta"=res$par, "est_obj"=mrc_Sn, 
              "gt_beta"=gt_beta, "gt_obj"=gt_Sn))
}

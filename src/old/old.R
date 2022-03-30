print("here are the codes that do not run currently")

Gbeta <- function(beta, ecdf_Y, X) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  sub_vals <- X %*% beta
  zhat <- as.vector(sapply(ecdf_Y, inverse_cdf_z, arg=sub_vals))
  
  res <- t(X) %*% (zhat - X %*% beta)
  return (res)
}
# FIXME values of matrix in Z are too large
Gbeta_grad <-function(beta, ecdf_Y, X, pert=1e-4) {
  n <- nrow(X)
  m <- ncol(X)

  X_bar <- X^(-1)

  # print(X_bar)

  Z_bar <- matrix(0, nrow=n, ncol=n)
  sub_vals <- X %*% beta
  zhat <- as.vector(sapply(ecdf_Y, inverse_cdf_z, arg=sub_vals))

  for(i in 1:n) {
    for(j in 1:n) {
      Z_bar[i, j] <- 1/(dnorm(zhat[j] - t(X[i, ]) %*% beta) + pert)
    }
  }
  print(Z_bar)
  grad_G_beta <- -n*t(X) %*% (Z_bar %*% X_bar - X)

  return (grad_G_beta)
}

compute_G_beta <- function(n) {
  data <- simulate_rank_regression_data(n, 6)
  
  X <- data$X
  Y <- data$Y
  
  ranks_Y <- rank(Y)
  ecdf_Y <- ranks_Y/(n+1)
  ecdf_Y
  
  beta <- data$beta
  g_beta <- Gbeta(ecdf_Y, X, beta)
  
  res <- list("gt_beta"=beta, "G_beta"=g_beta)
  return (res)
}

newton_root_finding <- function(f, start_point, grad_f, ecdf_Y, X, tol=1e-9, max_iter=1000) {
  cur_point <- start_point
  if(sum(f(cur_point, ecdf_Y, X) != rep(0, length(start_point))) == 0) {
    return (list("root_approx"=cur_point, "iterations"=0))
  }

  for (i in 1:max_iter) {
    cur_grad <- grad_f(cur_point, ecdf_Y, X)
    print(cur_grad)
    next_point <- cur_point - inv(cur_grad, ecdf_Y, X) %*% f(cur_point, ecdf_Y, X)

    print(next_point)
    if (sum((next_point - cur_point)**2)**0.5 < tol) {
      return (list("root_approx"=c(next_point), "iterations"=i))
    }

    print(c(next_point))
    cur_point <- c(next_point)
  }

  print("max iterations exceeded ...")
  return (list("root_approx"=cur_point, "iterations"=max_iter))
}

find_f1_coefs_fixed_point <- function(Y, X, tol=1e-9, max_iter=300) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  A <- t(X) %*% X
  P <- inv(A) %*% t(X)

  m <- ncol(X)
  n <- nrow(X)
  coefs <- runif(m, min=-10, max=10)

  ranks_Y <- rank(Y)
  empirical_cdf_Y <- ranks_Y/(n+1)

  err <- 1

  for(iter in 1:max_iter) {
    print(paste("iter  ", iter))

    sub_vals <- X %*% coefs
    yhat <- as.vector(sapply(empirical_cdf_Y, inverse_cdf_z, arg=sub_vals))

    a <- P %*% yhat
    err <- sum((a - coefs)**2)**0.5

    coefs <- a

    print(paste("error is ", err))
    print("coefs are ")
    print(coefs)
    if(err < tol) {
      print("converged")
      break
    }
  }

  return(coefs)
}
f <- function(x) pbeta(x, shape1=2, shape2=3)
f.inv <- inverse(f,lower=0,upper=1)
f.inv(.2)

# FIXME works only for small n(<40), otherwise u_beta vanishes
find_f1_coefs_conditional_monte_carlo_algorithm <- function(Y, X, M=10, max_iter=100,
                                                            tol=1e-9) {
  n <- nrow(X)
  m <- ncol(X)
  ranks_Y <- rank(Y)
  # centering a matrix column wise
  X <- X - matrix(rep(colMeans(X), n), n, m, byrow = T)

  inv_XTX <- NaN
  if (m == 1) {
    inv_XTX <- 1/(t(X) %*% X)
  } else {
    inv_XTX <- inv(t(X) %*% X)
  }
  M_mult <- inv_XTX %*% t(X)

  B_hat <- matrix(0, nrow=m, ncol=M)
  rss_hats <- c()

  for(i in 1:M) {
    # Generate standard normal sample, standardize by sample mean and standard
    # deviation, order and put into the same rank order as the original sample
    max_val <- 1000
    while(T) {
      z <- rnorm(n)
      z_hat <- (z - mean(z))/sd(z)
      z <- sort(z_hat)
      z_hat <- z_hat[ranks_Y]

      # finding b hat and ssf hat
      b_hat <- M_mult %*% z_hat
      ssf_hat <- sum((X %*% b_hat)**2)
      rss_hat <- n-1-ssf_hat
      print(rss_hat)
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
  v <- rss_hats**(-(n-1)/2)

  coefs <- runif(m, min=-10, max=10)
  print("rss_hats")
  print(rss_hats)
  print("B_hat")
  print(B_hat)
  print("TT")
  print(TT)
  print("v ")
  print(v)
  for(i in 1:max_iter) {
    U_beta <- X %*% (TT*B_hat - coefs)
    u_beta <- colSums(U_beta**2)

    u_beta <- exp(-u_beta/2)
    print("u_beta")
    print(u_beta)
    w <- v * u_beta
    coefs_next <- 1/sum(w) * (matrix(w, m, M, byrow = T)*B_hat)
    coefs_next <- rowSums(coefs_next)

    print("coefs next")
    print(coefs_next)
    err <- (sum(coefs_next - coefs)**2)**0.5
    coefs <- coefs_next

    print(paste("err ", err))
    if(err < tol) {
      print("converged ...")
      break
    }
  }

  return (coefs)
}
# 
# data <- simulate_rank_regression_data(30, 2)
# coefs <- find_f1_coefs_conditional_monte_carlo_algorithm(data$Y, data$X, M = 100, max_iter = 10)
# 
# coefs
# data$beta

# sub_vals <- c(2, 5, 6)
# val <- cdf_z(10.95, sub_vals)
# val
# sub_vals
# inverse_cdf_z(val, sub_vals)

# data <- simulate_rank_regression_data(1000, 6)
# X <- data$X
# head(X)
# ss <- X %*% data$beta
# hist(X %*% data$beta)
# ss <- sort(ss)
# lags <- ss[-1] - ss[1: (length(ss)-1)] 
# summary(lags)
# 
# res <- run_fixed_point(1000)
# res
# system.time(
#   res_betas <- mclapply(c(100, 300), run_fixed_point, mc.cores = 2)
# )

# numCores <- detectCores() - 1
# numCores
# system.time(
#   results <- mclapply(c(1000, 3000, 5000, 10000, 15000), compute_G_beta, mc.cores = 5)
# )
# results

# data <- simulate_rank_regression_data(10, 6)
# 
# X <- data$X
# 
# Y <- data$Y
# m <- ncol(X)
# n <- nrow(X)
# ranks_Y <- rank(Y)
# ecdf_Y <- ranks_Y/(n+1)
# beta <- data$beta
# 
# Gbeta(beta, ecdf_Y, X)
# Gbeta_grad(beta, ecdf_Y, X)
# 
# root_G_beta <- newton_root_finding(Gbeta, rep(1, m), Gbeta_grad, ecdf_Y, X)

# Gbeta(ecdf_Y, X, beta)
# zzz <- solve(t(X), t(X) %*% X %*% beta)
# t(X) %*% zzz - t(X) %*% X %*% beta


# Full hsic parts
kernel_matrix_func <- function(x, sigma_sq=1) {
  num_row <- nrow(x)
  K <- matrix(nrow = num_row, ncol = num_row)
  for(i in 1:num_row) {
    for(j in 1:i) {
      dist_sq <- sum((x[i, ]-x[j, ])^2)
      K[i, j] <- exp(-dist_sq/(2*sigma_sq))
      K[j, i] <- K[i, j]
    }
  }
  
  return(K)
}

dhsic_func <- function(x, y) {
  if(!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if(!is.matrix(y)) {
    y <- as.matrix(y)
  }
  Kx <- kernel_matrix_func(x)
  Ky <- kernel_matrix_func(y)
  num_row <- nrow(x)
  
  dist_hsic <- 1/num_row^2*sum(Kx * Ky) + 1/num_row^4*sum(Kx)*sum(Ky) - 
    2/num_row^3*sum(colSums(Kx)*colSums(Ky))
  
  return(dist_hsic)
}

learn_funcs_grad <- function(params, xj, x_rem) {
  nj <- get_noise(params, xj, x_rem)
  
  # to be moved out to the parent function
  if(!is.matrix(nj)) {
    nj <- as.matrix(nj)
  }
  if(!is.matrix(x_rem)) {
    x_rem <- as.matrix(x_rem)
  }
  Knj <- kernel_matrix_func(nj)
  Kx_rem <- kernel_matrix_func(x_rem)
  # end moved out part
  
  diff_kernel <- function(u, v) {
    return(u-v)
  }
  
  kernel_matrix_withfunc_func <- function(x, func){
    if(!is.matrix(x)) {
      x <- as.matrix(x)
    }
    num_row <- nrow(x)
    K <- matrix(nrow=num_row, ncol=num_row)
    for (i in 1:num_row) {
      for (j in 1:num_row) {
        K[i,j] <- match.fun(func)(x[i,], x[j,])
      }
    }
    return(K)
  }
  
  Nj_diff_kerenl <- kernel_matrix_withfunc_func(nj, diff_kernel)
  
  Xj_diff_kernel_1 <- kernel_matrix_withfunc_func(xj, diff_kernel)
  Xj_diff_kernel_2 <- kernel_matrix_withfunc_func(xj^2, diff_kernel)
  Xj_diff_kernel_3 <- kernel_matrix_withfunc_func(xj^3, diff_kernel)
  
  num_row <- nrow(xj)
  
  mat_for_grad_1 <- -2*Knj*Xj_diff_kernel_1*Nj_diff_kerenl
  param1_grad <- 1/num_row^2*sum(Kx_rem * mat_for_grad_1) + 
    1/num_row^4*sum(Kx_rem)*sum(mat_for_grad_1) - 
    2/num_row^3*sum(colSums(Kx_rem)*colSums(mat_for_grad_1))
  
  mat_for_grad_2 <- -2*Knj*Xj_diff_kernel_2*Nj_diff_kerenl
  param2_grad <- 1/num_row^2*sum(Kx_rem * mat_for_grad_2) + 
    1/num_row^4*sum(Kx_rem)*sum(mat_for_grad_2) - 
    2/num_row^3*sum(colSums(Kx_rem)*colSums(mat_for_grad_2))
  
  mat_for_grad_3 <- -2*Knj*Xj_diff_kernel_3*Nj_diff_kerenl
  param3_grad <- 1/num_row^2*sum(Kx_rem * mat_for_grad_3) + 
    1/num_row^4*sum(Kx_rem)*sum(mat_for_grad_3) - 
    2/num_row^3*sum(colSums(Kx_rem)*colSums(mat_for_grad_3))
  
  params_grad <- c(param1_grad, param2_grad, param3_grad)
  num_rem_nodes <- ncol(x_rem)
  for(k in 1:num_rem_nodes) {
    Xk_diff_kernel <- kernel_matrix_withfunc_func(x_rem[, k], diff_kernel)
    mat_for_grad_3plusk <- 2*Knj*Xk_diff_kernel*Nj_diff_kerenl
    param3plusk_grad <- 1/num_row^2*sum(Kx_rem * mat_for_grad_3plusk) + 
      1/num_row^4*sum(Kx_rem)*sum(mat_for_grad_3plusk) - 
      2/num_row^3*sum(colSums(Kx_rem)*colSums(mat_for_grad_3plusk))
    
    params_grad <- c(params_grad, param3plusk_grad)
  }
  
  for(k1 in 1:num_rem_nodes) {
    for(k2 in 1:num_rem_nodes) {
      Xk1k2_diff_kernel <- kernel_matrix_withfunc_func(x_rem[, k1]*x_rem[, k2], diff_kernel)
      mat_for_grad_3pluskk1k2 <- 2*Knj*Xk1k2_diff_kernel*Nj_diff_kerenl
      param3pluskk1k2_grad <- 1/num_row^2*sum(Kx_rem * mat_for_grad_3pluskk1k2) + 
        1/num_row^4*sum(Kx_rem)*sum(mat_for_grad_3pluskk1k2) - 
        2/num_row^3*sum(colSums(Kx_rem)*colSums(mat_for_grad_3pluskk1k2))
      
      params_grad <- c(params_grad, param3pluskk1k2_grad)
    }
  }
  
  return(params_grad)
}

# library(MonteCarlo)
#n_grid <- c(1000, 5000, 10000)
#fj2_func_grid <- c("square", "other")
#param_list=list("n"=n_grid, "fj2_func"=fj2_func_grid)
#mc_result <- MonteCarlo(func = simulate_recover, nrep = 100,  ncpus = 6,
#                        param_list = param_list)
#mc_result
# MakeTable(output = mc_result, rows = "n", cols = "fj2_func")

find_f1_coefs_expected_rank_algorithm <- function(Y, X, lamb=10, penalty="ell2") {
  print(paste("starting expected rank ", penalty))
  G_j_beta <- function(j, beta, X) {
    val <- 0
    for(i in 1:nrow(X)) {
      val <- val + pnorm(sum((X[j, ] - X[i, ])*beta)/2**0.5)
    }
    return (1/2 + val)
  }
  
  S_beta <- function(beta, X, ranks_Y) {
    val <- 0
    for(j in 1:nrow(X)) {
      val <- val + (ranks_Y[j] - G_j_beta(j, beta, X))**2
    }
    
    if (penalty=="ell1") {
      val <- val + lamb*sum(abs(beta))
    } else if(penalty=="ell2") {
      val <- val + lamb*sum((beta)**2)
    } else {
      print("no penalty")
    }
    
    return (val)
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
  
  est_beta <- optim(par=coefs, fn=S_beta, method = "BFGS", X=X, ranks_Y=ranks_Y)
  
  return(est_beta$par)
}

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


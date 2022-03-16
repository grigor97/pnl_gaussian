print("here are the codes that do not run currently")

# TODO maybe standardize the data
run_fixed_point <- function(n, batch_size=64, max_iter=100, tol = 1e-5) {
  data <- simulate_rank_regression_data(n, 6)

  X <- data$X
  Y <- data$Y
  beta <- data$beta

  ranks_Y <- rank(Y)
  ecdf_Y <- ranks_Y/(n+1)
  G_beta <- Gbeta(beta, ecdf_Y, X)

  pred_beta <- find_f1_coefs_fixed_point_stochastic(Y, X, batch_size = batch_size,
                                                    max_iter = max_iter, tol=tol)

  G_pred_beta <- Gbeta(pred_beta, ecdf_Y, X)

  res <- list("gt_beta"=beta, "pred_beta"=c(pred_beta), "G_beta"=c(G_beta),
              "G_pred_beta"=c(G_pred_beta),
              "n"=nrow(X), "m"=ncol(X), "max_iter"=max_iter,
              "batch_size"=batch_size, "tol"=tol,
              "l2_dist_betas"=sum((pred_beta - beta)**2)**0.5)
  return (res)
}

# run_fixed_point(100)

save_plots <- function(res_betas, file_n='/Users/grigorkeropyan/Desktop/Stat/codes/plots/batch_64/') {
  lgth <- length(res_betas)
  for(i in 1:lgth) {
    pred <- c(res_betas[[i]]$pred_beta)
    actual <- res_betas[[i]]$gt_beta
    title_name <- paste('sample size(n) ', res_betas[[i]]$n)
    title_name <- paste(title_name, ", batch size ", sep='')
    title_name <- paste(title_name, res_betas[[i]]$batch_size)
    title_name <- paste(title_name, ", max iter ", sep='')
    title_name <- paste(title_name, res_betas[[i]]$max_iter)
    pl <- ggplot() +
      geom_point(aes(x=1:6, y=actual, colour="ground truth"), size=1.5) +
      geom_point(aes(x=1:6, y=pred, colour="estimated"), size=1) +
      labs(title=title_name,
           x="beta1:6", y="ground truth and estimated beta's" ) +
      scale_color_manual(name="legend", breaks=c("ground truth", "estimated"),
                         values=c("ground truth"="black", "estimated"="red")) +
      theme(legend.position = c(0.1,0.91), plot.title = element_text(hjust = 0.5))

    file_name <- paste(file_n, 'betas_')
    file_name <- paste(file_name, res_betas[[i]]$n, sep='')
    file_name <- paste(file_name, '_')
    file_name <- paste(file_name, res_betas[[i]]$batch_size, sep='')
    file_name <- paste(file_name, '.png', sep='')
    ggsave(filename = file_name, plot = pl)

    pl_G <- ggplot() +
      geom_point(aes(x=1:6, y=res_betas[[i]]$G_beta, colour="G_beta"), size=1.5) +
      geom_point(aes(x=1:6, y=res_betas[[i]]$G_pred_beta, colour="G_pred_beta"), size=1) +
      labs(title=title_name,
           x="beta1:6", y="ground truth and estimated G_beta's" ) +
      scale_color_manual(name="legend", breaks=c("G_beta", "G_pred_beta"),
                         values=c("G_beta"="black", "G_pred_beta"="red")) +
      theme(legend.position = c(0.1,0.91), plot.title = element_text(hjust = 0.5))

    file_name_G <- paste(file_n, 'G_betas_')
    file_name_G <- paste(file_name_G, res_betas[[i]]$n, sep='')
    file_name_G <- paste(file_name_G, '_')
    file_name_G <- paste(file_name_G, res_betas[[i]]$batch_size, sep='')
    file_name_G <- paste(file_name_G, '.png', sep='')
    ggsave(filename = file_name_G, plot = pl_G)
  }
}

# save_plots(res_betas)

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

get_ranks <- function(df) {
  ranks <- data.frame(matrix(ncol=0, nrow=nrow(df)))
  for(col_name in colnames(df)) {
    ranks[col_name] <- rank(df[[col_name]])
  }
  return(ranks)
}

#' Title computes the value of G(beta)
#'
#' @param beta coefficients
#' @param ecdf_Y adjusted empirical CDF of Y
#' @param X matrix X with rows X_i corresponding to response variable Y_i
#'
#' @return value of function G in beta
#' @export
#'
#' @examples
Gbeta <- function(beta, ecdf_Y, X) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  sub_vals <- X %*% beta
  zhat <- as.vector(sapply(ecdf_Y, inverse_cdf_z, arg=sub_vals))

  res <- t(X) %*% (zhat - X %*% beta)
  return (res)
}

#FIXME works only for small n and M, otherwise the u_beta become zero (vanishing)
find_f1_coefs_monte_carlo_algorithm <- function(Y, X, M=10, max_iter=100) {
  ranks_Y <- rank(Y)
  n <- nrow(X)
  m <- ncol(X)
  coefs <- runif(m, min=-10, max=10)

  Z <- matrix(data=rnorm(M*n), nrow = n, ncol = M)
  Z <- apply(Z, 2, sort)
  Z <- Z[ranks_Y, ]

  inv_XTX <- NaN
  if (m == 1) {
    inv_XTX <- 1/(t(X) %*% X)
  } else {
    inv_XTX <- inv(t(X) %*% X)
  }

  B <-  inv_XTX %*% t(X) %*% Z
  Z_hat <- X %*% B

  ssf <- colSums(Z_hat**2)
  v <- exp(ssf/2)

  for(iter in 1:max_iter) {
    u_beta <- X %*% (B - coefs)
    u_beta <- colSums(u_beta**2)
    u_beta <- exp(-u_beta/2)
    w <- v * u_beta
    coefs <- 1/sum(w) * (matrix(w, m, M, byrow = T)*B)
    coefs <- rowSums(coefs)
  }

  return (coefs)
}

# data <- simulate_rank_regression_data(40, 2)
# coefs <- find_f1_coefs_monte_carlo_algorithm(data$Y, data$X, M = 20, max_iter = 3)
# 
# coefs
# data$beta

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

func2 <- function(x) {
  return(x^3 - 2^x - 5)
}

grad_func2 <- function(x) {
  diag <- 3*x^2 - 2
  return(diag(diag))
}

grad_func2(c(9, 7))
func2(c(9, 7))

uniroot(func2, c(2,3))

uniroot(func2, lower = 2, upper = 3)
newton_root_finding(func2, c(0, 0), grad_func2)

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


simulate_recover <- function(n, fj2_func) {
  print("start of data simulation ...")
  x1 <- rexp(n)
  noise1 = runif(n)
  noise2 = runif(n)
  noise3 = runif(n)
  x2 <- 200 + 100*sin(x1) + noise1
  x2 <- x2^(1/3) + 1
  x3 <- 200 + 100*cos(x1) + noise2
  x3 <- x3^(1/5) + 1
  if(fj2_func == "square") {
    x4 <- (x2 + x3 + noise3)^2
  } else {
    x4 <- x2 + x3 + noise3
  }
  
  data <- cbind(x1, x2, x3, x4)
  data <- data.frame(data)
  print("end of data simulation")
  
  order <- order_recovery_by_last_node(data)
  res <- list()
  i <- 0
  for(ord in order) {
    i = i+1
    node <- substring(ord, 2, 10)
    res[paste("node", toString(i), sep="")] = strtoi(node)
    
  }
  print("result of the simulation is")
  print(res)
  return(res)
}


ord <- simulate_recover(1000, "other")
ord

cn <- c(2000, 4000, 6000)
ress <- c()
for(n in cn) {
  ord <- simulate_recover(n, "other")
  print(ord)
  ress <- c(ress, ord)
}

ress 

# library(MonteCarlo)
#n_grid <- c(1000, 5000, 10000)
#fj2_func_grid <- c("square", "other")
#param_list=list("n"=n_grid, "fj2_func"=fj2_func_grid)
#mc_result <- MonteCarlo(func = simulate_recover, nrep = 100,  ncpus = 6,
#                        param_list = param_list)
#mc_result
# MakeTable(output = mc_result, rows = "n", cols = "fj2_func")

find_f1_coefs_expected_rank_l1_algorithm <- function(Y, X, lamb=10) {
  print("starting expected rank l1 algorithm")
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
    
    val <- val + lamb*sum(abs(beta))
    
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

source("rank_regression.R")
n = 500
for(beta in c(0.1, 1, 10, 100)) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  X <- matrix(rnorm(n), n, 1)
  noise <- rnorm(n)
  Y <- X %*% beta + noise
  Y <- exponent(Y, 1/3) + 4.7
  Y <- as.matrix(Y)

  acc <- sum(c(rank(X %*% beta)) == rank(Y))/length(Y)
  print("beta")
  print(beta)
  print("expected")
  print(acc)

  est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, lamb=0, penalty = "ell2")
  acc <- sum(c(rank(X %*% beta)) == rank(Y))/length(Y)
  print("expected rank l2 lambda 0")
  print(acc)

  est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, lamb=10, penalty = "ell2")
  acc <- sum(c(rank(X %*% beta)) == rank(Y))/length(Y)
  print("expected rank l2 lambda 10")
  print(acc)
}

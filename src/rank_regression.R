library(matlib)
library(parallel)

library(ggplot2)

set.seed(12)


#' Title CDF  function for variable z 
#'
#' @param y value of y
#' @param subtract_values values that are subtracted from y 
#'
#' @return value of the CDF function
#' @export
#'
#' @examples
cdf_z <- function(y, subtract_values) {
  res = 0
  for(i in subtract_values) {
    res = res + pnorm(y - i)
  }
  res = res / length(subtract_values)
  return(res)
}

#' Title
#'
#' @param f function for root finding
#' @param lower lower bound of the root
#' @param upper upper bound of the root
#'
#' @return root corresponding to f(x, arg) - y = 0
#' @export
#'
#' @examples 
inverse <- function(f, lower, upper){
  function(y, arg){
    uniroot(function(x){f(x, arg) - y}, lower = lower, upper = upper, tol=1e-3,)[1]$root
  }
}


inverse_cdf_z<- inverse(cdf_z, -1000, 1000)

#' Title finds beta coefficients using a fixed point iteration for a 
#' whole data. Should be called only for small sizes of n(sample size) < 1000.
#'
#' @param Y n i.i.d. response variables Y_i = f_2^{-1}(beta^T x_i + error_i)
#' @param X matrix X with rows X_i corresponding to Y_i
#' @param tol tolerance, stop the algorithm when th error is smaller than @tol
#' @param max_iter maximum number of iteration for the algorithm
#'
#' @return returns the values of coefficients
#' @export 
#'
#' @examples 
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

#' Title dinds beta coefficients using a stochastic fixed point iteration. 
#' Can be called for large sizes of n(sample size).
#'
#' @param Y n i.i.d. response variables Y_i = f_2^{-1}(beta^T x_i + error_i)
#' @param X matrix X with rows X_i corresponding to Y_i
#' @param batch_size batch size for every iteration. Algorithm chooses randomly 
#' number of batch size elements and do one step fixed point iteration
#' @param tol tolerance, stop the algorithm when th error is smaller than @tol
#' @param max_iter maximum number of iteration for the algorithm
#'
#' @return returns the obtained values of coefficients
#' @export
#'
#' @examples
find_f1_coefs_fixed_point_stochastic <- function(Y, X, batch_size=64, tol=1e-9, 
                                                 max_iter=100) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
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
      P <- inv(A) %*% t(batch_X)
    } else {
      P <- 1/A[1, 1] * t(batch_X)
    }
    ranks_Y <- rank(batch_Y)
    empirical_cdf_Y <- ranks_Y/(batch_size+1)
    
    # print(paste("iter  ", iter))
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

#' Title simulate rank regression data by random coefficients
#'
#' @param n number of samples
#' @param m dimension of the coefficients
#'
#' @return returns a dataset 
#' @export
#'
#' @examples
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

#FIXME check the algorithm
find_f1_coefs_expected_rank_algorithm <- function(Y, X, max_iter=100) {
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
    
    return (val)
  }
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  m <- ncol(X)
  n <- nrow(X)
  coefs <- matrix(runif(m, min=-10, max=10), m, 1)
  ranks_Y <- rank(Y)
  
  est_beta <- optim(par=coefs, fn=S_beta, method = "BFGS", 
                     X=X, ranks_Y=ranks_Y)
  
  return(est_beta$par)
}

data <- simulate_rank_regression_data(1000, 1)
system.time(
  coefs <- find_f1_coefs_expected_rank_algorithm(data$Y, data$X)
)
coefs
data$beta


find_f1_coefs_monte_carlo_algorithm <- function(Y, X, M=10, max_iter=100) {
  ranks_Y <- rank(Y)
  n <- nrow(X)
  m <- ncol(X)
  coefs <- matrix(runif(m, min=-10, max=10), m, 1)
  for(iter in 1:max_iter) {
    
  }
}


#' Title run rank regression algorithms 
#' for specific arguments and number of datasets
#'
#' @param n number of samples in each dataset
#' @param m number of coefficients
#' @param number_of_datasets number of datasets
#'
#' @return ground truth betas and estimated betas
#' @export
#'
#' @examples
run_rank_regression_algorithms <- function(n=1000, m=1, number_of_datasets=100){
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  betas <- c(0.01, 0.1, 1, 10, 30, 100, 300, 1000)
  
  rank_reg_est_betas <- matrix(0, nrow=number_of_datasets, ncol=length(betas))
  for(i in 1:number_of_datasets) {
    print(paste("dataset ", i))
    X <- matrix(rnorm(n*m), n, m)
    for (j in 1:length(betas)) {
      print(paste("beta ", j))
      beta <- betas[j]
      noise <- rnorm(n)
      Y <- X %*% beta + noise
      Y <- exponent(Y, 1/3) + 4.7
      
      ranks_Y <- rank(Y)
      ecdf_Y <- ranks_Y/(n+1)
      
      pred_beta <- find_f1_coefs_fixed_point_stochastic(Y, X, batch_size = 64, 
                                                        max_iter = 100, tol=1e-9)
      
      rank_reg_est_betas[i, j] <- pred_beta
    }
    EST_betas <- rank_reg_est_betas
  }
  return (list("rank_reg_est_betas"=rank_reg_est_betas, "betas"=betas))
}

EST_betas <- 0

#' Title run rank regression algorithms and 
#' save a boxplots of the estimated betas
#'
#' @param file_n file name where the plots should be saved
#'
#' @return ground truth betas and estimated betas
#' @export
#'
#' @examples
run_algorithms_and_save_plots <- function(file_n='/Users/grigorkeropyan/pnl_gaussian/plots/'){
  res <- run_rank_regression_algorithms()
  rank_reg_est_betas <- res$rank_reg_est_betas
  betas <- res$betas
  
  dat <- stack(as.data.frame(rank_reg_est_betas))
  
  pl <- ggplot() + geom_boxplot(aes(x=dat$ind, y=dat$values, colour='estimated betas')) + 
    geom_point(aes(x=unique(dat$ind), y=betas, colour='ground truth betas')) + 
    labs(title="rank regression", x="",y="betas") +
    scale_color_manual(name='legend', 
                       breaks=c('estimated betas', 'ground truth betas'),
                       values=c('black', 'red')) +
    guides(colour = guide_legend(override.aes = list(
      linetype = c("solid", "blank"),
      color = c("black","red")
    ))) +
    # theme(legend.position=c(0.15,0.91), plot.title = element_text(hjust = 0.5))
    theme(legend.position='top', plot.title = element_text(hjust = 0.5))
  
  # file_name <- paste(file_n, 'rank_reg_alg/', sep='')
  file_name <- paste(file_n, 'rank_reg_box_plots.png', sep='')
  ggsave(filename = file_name, plot = pl)
  
  plot(pl)
  
  return (res)
}

res <- run_algorithms_and_save_plots()
res

# data <- simulate_rank_regression_data(10000, 6)
# saveRDS(data, '/Users/grigorkeropyan/pnl_gaussian/datasets/data1')
# data <- readRDS('/Users/grigorkeropyan/pnl_gaussian/datasets/data1')
# data

# TODO maybe standardize the data
# run_fixed_point <- function(n, batch_size=64, max_iter=100, tol = 1e-5) {
#   data <- simulate_rank_regression_data(n, 6)
#   
#   X <- data$X
#   Y <- data$Y
#   beta <- data$beta
#   
#   ranks_Y <- rank(Y)
#   ecdf_Y <- ranks_Y/(n+1)
#   G_beta <- Gbeta(beta, ecdf_Y, X)
#   
#   pred_beta <- find_f1_coefs_fixed_point_stochastic(Y, X, batch_size = batch_size, 
#                                                     max_iter = max_iter, tol=tol)
#   
#   G_pred_beta <- Gbeta(pred_beta, ecdf_Y, X)
#   
#   res <- list("gt_beta"=beta, "pred_beta"=c(pred_beta), "G_beta"=c(G_beta), 
#               "G_pred_beta"=c(G_pred_beta),
#               "n"=nrow(X), "m"=ncol(X), "max_iter"=max_iter, 
#               "batch_size"=batch_size, "tol"=tol,
#               "l2_dist_betas"=sum((pred_beta - beta)**2)**0.5)
#   return (res)
# }

# run_fixed_point(100)
# 
# system.time(
#   res_betas <- mclapply(c(1000, 3000, 5000, 10000, 15000), run_fixed_point, mc.cores = 5)
# )
# 
# res_betas

# save_plots <- function(res_betas, file_n='/Users/grigorkeropyan/Desktop/Stat/codes/plots/batch_64/') {
#   lgth <- length(res_betas)
#   for(i in 1:lgth) {
#     pred <- c(res_betas[[i]]$pred_beta)
#     actual <- res_betas[[i]]$gt_beta
#     title_name <- paste('sample size(n) ', res_betas[[i]]$n)
#     title_name <- paste(title_name, ", batch size ", sep='')
#     title_name <- paste(title_name, res_betas[[i]]$batch_size)
#     title_name <- paste(title_name, ", max iter ", sep='')
#     title_name <- paste(title_name, res_betas[[i]]$max_iter)
#     pl <- ggplot() + 
#       geom_point(aes(x=1:6, y=actual, colour="ground truth"), size=1.5) +
#       geom_point(aes(x=1:6, y=pred, colour="estimated"), size=1) + 
#       labs(title=title_name, 
#            x="beta1:6", y="ground truth and estimated beta's" ) +
#       scale_color_manual(name="legend", breaks=c("ground truth", "estimated"),
#                          values=c("ground truth"="black", "estimated"="red")) +
#       theme(legend.position = c(0.1,0.91), plot.title = element_text(hjust = 0.5))
#     
#     file_name <- paste(file_n, 'betas_')
#     file_name <- paste(file_name, res_betas[[i]]$n, sep='')
#     file_name <- paste(file_name, '_')
#     file_name <- paste(file_name, res_betas[[i]]$batch_size, sep='')
#     file_name <- paste(file_name, '.png', sep='')
#     ggsave(filename = file_name, plot = pl)
#     
#     pl_G <- ggplot() + 
#       geom_point(aes(x=1:6, y=res_betas[[i]]$G_beta, colour="G_beta"), size=1.5) +
#       geom_point(aes(x=1:6, y=res_betas[[i]]$G_pred_beta, colour="G_pred_beta"), size=1) + 
#       labs(title=title_name, 
#            x="beta1:6", y="ground truth and estimated G_beta's" ) +
#       scale_color_manual(name="legend", breaks=c("G_beta", "G_pred_beta"),
#                          values=c("G_beta"="black", "G_pred_beta"="red")) +
#       theme(legend.position = c(0.1,0.91), plot.title = element_text(hjust = 0.5))
#     
#     file_name_G <- paste(file_n, 'G_betas_')
#     file_name_G <- paste(file_name_G, res_betas[[i]]$n, sep='')
#     file_name_G <- paste(file_name_G, '_')
#     file_name_G <- paste(file_name_G, res_betas[[i]]$batch_size, sep='')
#     file_name_G <- paste(file_name_G, '.png', sep='')
#     ggsave(filename = file_name_G, plot = pl_G)
#   }
# }
# 
# save_plots(res_betas)

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

# FIXME values of matrix in Z are too large
# Gbeta_grad <-function(beta, ecdf_Y, X, pert=1e-4) {
#   n <- nrow(X)
#   m <- ncol(X)
#   
#   X_bar <- X^(-1)
#   
#   # print(X_bar)
#   
#   Z_bar <- matrix(0, nrow=n, ncol=n)
#   sub_vals <- X %*% beta
#   zhat <- as.vector(sapply(ecdf_Y, inverse_cdf_z, arg=sub_vals))
#   
#   for(i in 1:n) {
#     for(j in 1:n) {
#       Z_bar[i, j] <- 1/(dnorm(zhat[j] - t(X[i, ]) %*% beta) + pert)
#     }
#   }
#   print(Z_bar)
#   grad_G_beta <- -n*t(X) %*% (Z_bar %*% X_bar - X)
#   
#   return (grad_G_beta)
# }
# 
# newton_root_finding <- function(f, start_point, grad_f, ecdf_Y, X, tol=1e-9, max_iter=1000) {
#   cur_point <- start_point
#   if(sum(f(cur_point, ecdf_Y, X) != rep(0, length(start_point))) == 0) {
#     return (list("root_approx"=cur_point, "iterations"=0))
#   }
#   
#   for (i in 1:max_iter) {
#     cur_grad <- grad_f(cur_point, ecdf_Y, X)
#     print(cur_grad)
#     next_point <- cur_point - inv(cur_grad, ecdf_Y, X) %*% f(cur_point, ecdf_Y, X)
#     
#     print(next_point)
#     if (sum((next_point - cur_point)**2)**0.5 < tol) {
#       return (list("root_approx"=c(next_point), "iterations"=i))
#     }
#     
#     print(c(next_point))
#     cur_point <- c(next_point)
#   }
#   
#   print("max iterations exceeded ...")
#   return (list("root_approx"=cur_point, "iterations"=max_iter))
# }

# compute_G_beta <- function(n) {
#   data <- simulate_rank_regression_data(n, 6)
#   
#   X <- data$X
#   Y <- data$Y
#   
#   ranks_Y <- rank(Y)
#   ecdf_Y <- ranks_Y/(n+1)
#   ecdf_Y
#   
#   beta <- data$beta
#   g_beta <- Gbeta(ecdf_Y, X, beta)
#   
#   res <- list("gt_beta"=beta, "G_beta"=g_beta)
#   return (res)
# }

# func2 <- function(x) {
#   return(x^3 - 2^x - 5)
# }
# 
# grad_func2 <- function(x) {
#   diag <- 3*x^2 - 2
#   return(diag(diag))
# }
# 
# grad_func2(c(9, 7))
# func2(c(9, 7))
# 
# uniroot(func2, c(2,3))
# 
# uniroot(func2, lower = 2, upper = 3)
# newton_root_finding(func2, c(0, 0), grad_func2)

# sub_vals <- c(2, 5, 6)
# val <- cdf_z(10.95, sub_vals)
# val
# sub_vals
# inverse_cdf_z(val, sub_vals)

# get_ranks <- function(df) {
#   ranks <- data.frame(matrix(ncol=0, nrow=nrow(df)))
#   for(col_name in colnames(df)) {
#     ranks[col_name] <- rank(df[[col_name]])
#   }
#   return(ranks)
# }

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
# Gbeta <- function(beta, ecdf_Y, X) {
#   if (!is.matrix(X)) {
#     X <- as.matrix(X)
#   }
#   
#   sub_vals <- X %*% beta
#   zhat <- as.vector(sapply(ecdf_Y, inverse_cdf_z, arg=sub_vals))
#   
#   res <- t(X) %*% (zhat - X %*% beta)
#   return (res)
# }

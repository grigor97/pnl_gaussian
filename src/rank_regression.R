library(parallel)
library(rjson)
# library(ggplot2)

set.seed(12)

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
find_f1_coefs_fixed_point_stochastic <- function(Y, X, batch_size=64, 
                                                 max_iter=100, tol=1e-9) {
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
  
  est_beta <- optim(par=coefs, fn=S_beta, method = "BFGS", X=X, ranks_Y=ranks_Y)
  
  return(est_beta$par)
}



run_rank_regression_algorithms <- function(n=1000, m=1){
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  betas <- c(0.01, 0.1, 1)#, 0.3, 0.5) #, 0.7, 0.9, 1, 10, 30, 100)
  
  fixed_est_betas <- matrix(0, nrow=m, ncol=length(betas))
  expected_est_betas <- matrix(0, nrow=m, ncol=length(betas))
  
  # generate data
  X <- matrix(rnorm(n*m), n, m)
  
  for (j in 1:length(betas)) {
    print(paste("beta ", j))
    beta <- betas[j]
    
    # generate noise and corresponding y
    noise <- rnorm(n)
    Y <- X %*% beta + noise
    Y <- exponent(Y, 1/3) + 4.7
    
    fixed_pred_beta <- find_f1_coefs_fixed_point_stochastic(Y, X, batch_size=4,
                                                      max_iter=10, tol=1e-9)
    
    expected_pred_beta <- find_f1_coefs_expected_rank_algorithm(Y, X)
    
    fixed_est_betas[ , j] <- fixed_pred_beta
    expected_est_betas[ , j] <- expected_pred_beta
  }
  return (list("fixed_est_betas"=fixed_est_betas, 
               "expected_est_betas"=expected_est_betas,
               "betas"=betas))
}

dummy_fun <- function(i) {
  res <- run_rank_regression_algorithms(n=100, m=1)
  return(res)
}

# res <- run_rank_regression_algorithms(n=100, m=1)
# res

numCores <- detectCores() - 1
system.time(
  results <- mclapply(c(seq(1, 7)), dummy_fun, mc.cores = numCores)
)
print("finished algorithms")

parse_result_and_save <- function(results) {
  betas <- results[[1]]$betas

  fixed_betas <- matrix(0, length(results), length(betas))
  exp_betas <- matrix(0, length(results), length(betas))

  for(i in 1:length(results)) {
    fixed_betas[i, ] <- results[[i]]$fixed_est_betas
    exp_betas[i, ] <- results[[i]]$expected_est_betas
  }
  
  res <- list("betas"=betas, "fixed_betas"=fixed_betas, "exp_betas"=exp_betas,
              "num_datasets"=nrow(fixed_betas), "num_betas"=length(betas))
  json_data <- toJSON(res)
  write(json_data, "../res/all_betas")
  
  return(res)
}

results

parse_result_and_save(results)




# for local run
# ress <- fromJSON(file = "/Users/grigorkeropyan/pnl_gaussian/res/all_betas")
# ress
# 
# exp_betas <- matrix(ress$exp_betas, ress$num_datasets, ress$num_betas)
# fixed_betas <- matrix(ress$fixed_betas, ress$num_datasets, ress$num_betas)
# betas <- ress$betas
# # 
# # all.equal(exp_betas, res$exp_betas)
# # all.equal(fixed_betas, res$fixed_betas)
# 
# save_plots <- function(estimated_betas, gt_betas, alg_name,
#                        file_n='/Users/grigorkeropyan/pnl_gaussian/res/'){
#   number_of_datasets <- nrow(estimated_betas)
#   stacked_vals <- stack(as.data.frame(estimated_betas))
# 
#   title_name <- paste(alg_name, 'rank regression for')
#   title_name <- paste(title_name, number_of_datasets)
#   title_name <- paste(title_name, "datasets")
#   pl <- ggplot() + geom_boxplot(aes(x=stacked_vals$ind, y=stacked_vals$values, colour='estimated betas')) +
#     geom_point(aes(x=unique(stacked_vals$ind), y=gt_betas, colour='ground truth betas')) +
#     labs(title=title_name, x="",y="betas") +
#     scale_color_manual(name='',
#                        breaks=c('estimated betas', 'ground truth betas'),
#                        values=c('black', 'red')) +
#     guides(colour = guide_legend(override.aes = list(
#       linetype = c("solid", "blank"),
#       color = c("black","red")
#     ))) +
#     # theme(legend.position=c(0.15,0.91), plot.title = element_text(hjust = 0.5))
#     theme(legend.position='top', plot.title = element_text(hjust = 0.5))
#   file_name <- paste(file_n, alg_name, sep='')
#   file_name <- paste(file_name, 'reg_box_plots.png', sep='')
#   ggsave(filename = file_name, plot = pl)
# 
#   return(pl)
# }
# 
# library(ggplot2)
# save_plots(fixed_betas, betas, 'fixed_point')
# save_plots(exp_betas, betas, 'expected_rank')





# res$betas
# means_est_betas <- colMeans(res$rank_reg_est_betas)
# pl1 <- ggplot() + geom_line(aes(x=res$betas[1:10], y=res$betas[1:10] - means_est_betas[1:10]))+
#   labs(title = "rank regression for 100 datasets", 
#        x="ground truth betas", y="l_2 dist of betas") +
#   theme(plot.title = element_text(hjust = 0.5))
# pl1
# 
# ggsave(filename = '/Users/grigorkeropyan/pnl_gaussian/plots/rank_reg_diff_plot_100.png', plot = pl1)

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
# 
# #FIXME works only for small n and M, otherwise the u_beta become zero (vanishing)
# find_f1_coefs_monte_carlo_algorithm <- function(Y, X, M=10, max_iter=100) {
#   ranks_Y <- rank(Y)
#   n <- nrow(X)
#   m <- ncol(X)
#   coefs <- runif(m, min=-10, max=10)
#   
#   Z <- matrix(data=rnorm(M*n), nrow = n, ncol = M)
#   Z <- apply(Z, 2, sort)
#   Z <- Z[ranks_Y, ]
#   
#   inv_XTX <- NaN
#   if (m == 1) {
#     inv_XTX <- 1/(t(X) %*% X)
#   } else {
#     inv_XTX <- inv(t(X) %*% X)
#   }
#   
#   B <-  inv_XTX %*% t(X) %*% Z
#   Z_hat <- X %*% B
#   
#   ssf <- colSums(Z_hat**2)
#   v <- exp(ssf/2)
#   
#   for(iter in 1:max_iter) {
#     u_beta <- X %*% (B - coefs)
#     u_beta <- colSums(u_beta**2)
#     u_beta <- exp(-u_beta/2)
#     w <- v * u_beta
#     coefs <- 1/sum(w) * (matrix(w, m, M, byrow = T)*B)
#     coefs <- rowSums(coefs)
#   }
#   
#   return (coefs)
# }
# 
# data <- simulate_rank_regression_data(40, 2)
# coefs <- find_f1_coefs_monte_carlo_algorithm(data$Y, data$X, M = 20, max_iter = 3)
# 
# coefs
# data$beta


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
# find_f1_coefs_fixed_point <- function(Y, X, tol=1e-9, max_iter=300) {
#   if (!is.matrix(X)) {
#     X <- as.matrix(X)
#   }
#   A <- t(X) %*% X
#   P <- inv(A) %*% t(X)
# 
#   m <- ncol(X)
#   n <- nrow(X)
#   coefs <- runif(m, min=-10, max=10)
# 
#   ranks_Y <- rank(Y)
#   empirical_cdf_Y <- ranks_Y/(n+1)
# 
#   err <- 1
# 
#   for(iter in 1:max_iter) {
#     print(paste("iter  ", iter))
# 
#     sub_vals <- X %*% coefs
#     yhat <- as.vector(sapply(empirical_cdf_Y, inverse_cdf_z, arg=sub_vals))
# 
#     a <- P %*% yhat
#     err <- sum((a - coefs)**2)**0.5
# 
#     coefs <- a
# 
#     print(paste("error is ", err))
#     print("coefs are ")
#     print(coefs)
#     if(err < tol) {
#       print("converged")
#       break
#     }
#   }
# 
#   return(coefs)
# }
# f <- function(x) pbeta(x, shape1=2, shape2=3)
# f.inv <- inverse(f,lower=0,upper=1)
# f.inv(.2)

# FIXME works only for small n(<40), otherwise u_beta vanishes
# find_f1_coefs_conditional_monte_carlo_algorithm <- function(Y, X, M=10, max_iter=100, 
#                                                             tol=1e-9) {
#   n <- nrow(X)
#   m <- ncol(X)
#   ranks_Y <- rank(Y)
#   # centering a matrix column wise
#   X <- X - matrix(rep(colMeans(X), n), n, m, byrow = T)
#   
#   inv_XTX <- NaN
#   if (m == 1) {
#     inv_XTX <- 1/(t(X) %*% X)
#   } else {
#     inv_XTX <- inv(t(X) %*% X)
#   }
#   M_mult <- inv_XTX %*% t(X)
#   
#   B_hat <- matrix(0, nrow=m, ncol=M)
#   rss_hats <- c()
#   
#   for(i in 1:M) {
#     # Generate standard normal sample, standardize by sample mean and standard 
#     # deviation, order and put into the same rank order as the original sample
#     max_val <- 1000
#     while(T) {
#       z <- rnorm(n)
#       z_hat <- (z - mean(z))/sd(z)
#       z <- sort(z_hat)
#       z_hat <- z_hat[ranks_Y]
#       
#       # finding b hat and ssf hat
#       b_hat <- M_mult %*% z_hat
#       ssf_hat <- sum((X %*% b_hat)**2)
#       rss_hat <- n-1-ssf_hat
#       print(rss_hat)
#       if(rss_hat > max_val) {
#         next
#       }
#       
#       B_hat[, i] <- b_hat
#       rss_hats <- c(rss_hats, rss_hat)
#       break
#     }
#   }
#   
#   t <- ((n-3)/rss_hats)**0.5
#   TT <- matrix(t, m, M, byrow = T)
#   v <- rss_hats**(-(n-1)/2)
#   
#   coefs <- runif(m, min=-10, max=10)
#   print("rss_hats")
#   print(rss_hats)
#   print("B_hat")
#   print(B_hat)
#   print("TT")
#   print(TT)
#   print("v ")
#   print(v)
#   for(i in 1:max_iter) {
#     U_beta <- X %*% (TT*B_hat - coefs)
#     u_beta <- colSums(U_beta**2)
#     
#     u_beta <- exp(-u_beta/2)
#     print("u_beta")
#     print(u_beta)
#     w <- v * u_beta
#     coefs_next <- 1/sum(w) * (matrix(w, m, M, byrow = T)*B_hat)
#     coefs_next <- rowSums(coefs_next)
#     
#     print("coefs next")
#     print(coefs_next)
#     err <- (sum(coefs_next - coefs)**2)**0.5
#     coefs <- coefs_next
#     
#     print(paste("err ", err))
#     if(err < tol) {
#       print("converged ...")
#       break
#     }
#   }
#   
#   return (coefs)
# }
# 
# data <- simulate_rank_regression_data(30, 2)
# coefs <- find_f1_coefs_conditional_monte_carlo_algorithm(data$Y, data$X, M = 100, max_iter = 10)
# 
# coefs
# data$beta
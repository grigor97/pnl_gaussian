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
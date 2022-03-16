source("rank_regression/rank_regression.R")
library(dHSIC)

find_estimates <- function(Y, X, rank_alg, lamb) {
  est_beta <- NaN
  if(rank_alg=="expected_l2_rank") {
    est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, 
                                                      lamb = lamb, penalty = "ell2")
  } else if(rank_alg=="expected_l1_rank") {
    est_beta <- find_f1_coefs_expected_rank_algorithm(Y, X, 
                                                      lamb = lamb, penalty = "ell1")
  } else if(rank_alg=="fixed_point") {
    est_beta <- find_f1_coefs_fixed_point_stochastic(Y, X)
  } else {
    print("no such rank regression algorithm. Options are expected_l2_rank, 
          expected_l1_rank and fixed_point")
    return()
  }
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  ranks_Y <- rank(Y)
  empirical_cdf_Y <- ranks_Y/(length(Y)+1)
  sub_vals <- X %*% est_beta
  # print("estimate f_2_inv_y")
  
  f_2_inv_y <- as.vector(sapply(empirical_cdf_Y, inverse_cdf_z, arg=sub_vals))
  est_errors <- f_2_inv_y - sub_vals
  
  # print("compute dhsic value")
  dhsic_val <- dhsic(X, est_errors, kernel = "gaussian.fixed")$dHSIC
  
  res <- list("est_beta"=est_beta, "est_errors"=est_errors, "dhsic_val"=dhsic_val)
  
  return(res)
}

find_last_node_rank_reg <- function(data, nodes, rank_alg, lamb) 
{
  min_dhsic <- 1000
  last_node <- 0
  for(node in nodes) {
    rem_nodes <- nodes[!nodes %in% node]
    
    print("----------------")
    print(node)
    print("versus")
    print(rem_nodes)
    
    Y <- data[node]
    X <- data[rem_nodes]
    
    res <- find_estimates(Y, X, rank_alg, lamb)
    
    dhsic_val <- res$dhsic_val
    
    if(dhsic_val < min_dhsic) {
      min_dhsic <- dhsic_val
      print("min dhsic is updated")
      last_node <- node
    }
    print(dhsic_val)
  }
  
  return(last_node)
}

order_recovery_rank_reg <- function(data, rank_alg, lamb) 
{
  print("start of recovering the order from data using rank regression")
  order_of_nodes <- c()
  cur_nodes <- names(data)
  
  while(length(cur_nodes) > 1) {
    last_node <- find_last_node_rank_reg(data, cur_nodes, rank_alg, lamb)
    order_of_nodes <- c(last_node, order_of_nodes)
    cur_nodes <- cur_nodes[!cur_nodes %in% last_node]
  }
  
  order_of_nodes <- c(cur_nodes[1], order_of_nodes)
  print("end of recovering the order from data")
  return(order_of_nodes)
}
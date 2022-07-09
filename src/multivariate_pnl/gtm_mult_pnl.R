source("general_transformation_models/fully_non_param_estimation_gtm.R")
library(dHSIC)

min.dhsic.val <- function(Y, X) {
  res_noise_est <- noise.est.gtm(Y, X)
  est_noise <- res_noise_est$est_noise
  est_f2_inv <- res_noise_est$est_f2_inv
  est_g <- res_noise_est$est_g
  
  dhsic_val <- dhsic(X, est_noise, kernel = "gaussian.fixed")$dHSIC
  
  res <- list("est_g"=est_g, "est_noise"=est_noise, 
              "est_f2_inv"=est_f2_inv, "dhsic_val"=dhsic_val)
  
  return(res)
}

find_last_node_gtm <- function(data, nodes) 
{
  min_dhsic <- 1000
  last_node <- 0
  for(node in nodes) {
    rem_nodes <- nodes[!nodes %in% node]
    
    print("----------------")
    print(node)
    print("versus")
    print(rem_nodes)
    
    Y <- as.matrix(data[node])
    X <- as.matrix(data[rem_nodes])
    
    res <- min.dhsic.val(Y, X)
    
    dhsic_val <- res$dhsic_val
    
    if(dhsic_val < min_dhsic) {
      min_dhsic <- dhsic_val
      print("min dhsic is updated")
      last_node <- node
    }
    print("dhsic_val")
    print(dhsic_val)
  }
  
  return(last_node)
}

order_recovery_gtm <- function(data) 
{
  print("start of recovering the order from data using gtm")
  order_of_nodes <- c()
  cur_nodes <- names(data)
  
  while(length(cur_nodes) > 1) {
    last_node <- find_last_node_gtm(data, cur_nodes)
    order_of_nodes <- c(last_node, order_of_nodes)
    cur_nodes <- cur_nodes[!cur_nodes %in% last_node]
  }
  
  order_of_nodes <- c(cur_nodes[1], order_of_nodes)
  print("end of recovering the order from data")
  return(order_of_nodes)
}

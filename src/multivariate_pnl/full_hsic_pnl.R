library(dHSIC)

get_noise <- function(params, xj, x_rem) {
  num_cols = ncol(x_rem)
  fj2_inverse_xj = params[1]*xj + params[2]*xj^2 + params[3]*xj^3
  fj1_others_1 = 0
  fj1_others_2 = 0
  for(i in 1:num_cols) {
    fj1_others_1 = fj1_others_1 + params[3+i]*x_rem[i]
  }
  for(i in 1:num_cols) {
    for (j in 1:num_cols) {
      fj1_others_2 = fj1_others_2 + params[3+num_cols+ (i-1)*num_cols + j]*x_rem[i]*x_rem[j]
    }
  }
  nj = fj2_inverse_xj - fj1_others_1 - fj1_others_2
  
  return(nj)
}

learn_funcs <- function(params, xj, x_rem) 
{
  nj <- get_noise(params, xj, x_rem)
  res = dhsic(x_rem, nj, kernel = "gaussian.fixed")
  if(res$dHSIC < 0) {
    print("Erorrrrr ------ dHSIC is negativeeeeeee")
  }
  return(res$dHSIC)
}

find_last_node_naive <- function(data, nodes) 
{
  min_dhsic <- 1000
  last_node <- 0
  for(node in nodes) {
    rem_nodes <- nodes[!nodes %in% node]
    
    print("----------------")
    print(node)
    print("versus")
    print(rem_nodes)
    
    # FIXME check the parameters
    num_cols <- length(rem_nodes)
    num_params <- 3 + num_cols*(num_cols + 1)
    
    params <- runif(num_params)
    dhsic_val <- optim(par=params, fn=learn_funcs, method = "BFGS", 
                       xj=data[node], x_rem=data[rem_nodes])
    
    if(dhsic_val$value < min_dhsic) {
      min_dhsic <- dhsic_val$value
      print("min dhsic is updated")
      last_node <- node
    }
    print(dhsic_val)
  }
  return(last_node)
}

order_recovery_by_last_node <- function(data) 
{
  print("start of recovering the order from data")
  order_of_nodes <- c()
  cur_nodes <- names(data)
  
  while(length(cur_nodes) > 1) {
    last_node <- find_last_node_naive(data, cur_nodes)
    order_of_nodes <- c(last_node, order_of_nodes)
    cur_nodes <- cur_nodes[!cur_nodes %in% last_node]
  }
  
  order_of_nodes <- c(cur_nodes[1], order_of_nodes)
  print("end of recovering the order from data")
  return(order_of_nodes)
}

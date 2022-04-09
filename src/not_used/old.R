print("here are the codes that do not run currently")

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

# library(MonteCarlo)
#n_grid <- c(1000, 5000, 10000)
#fj2_func_grid <- c("square", "other")
#param_list=list("n"=n_grid, "fj2_func"=fj2_func_grid)
#mc_result <- MonteCarlo(func = simulate_recover, nrep = 100,  ncpus = 6,
#                        param_list = param_list)
#mc_result
# MakeTable(output = mc_result, rows = "n", cols = "fj2_func")

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


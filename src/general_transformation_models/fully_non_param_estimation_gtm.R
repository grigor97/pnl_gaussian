# implementation of https://projecteuclid.org/journals/bernoulli/volume-25/issue-4B/Estimation-of-fully-nonparametric-transformation-models/10.3150/19-BEJ1110.short
# here assuming the model f_2^{-1}(Y) = f_1(X) + epsilon, 
# where f_2, f_1 and epsilon are completly unknown
simulate.gtm.data <- function(n) {
  # as stated in the paper
  # FIXME add other transformation functions
  f2_inv <- function(y, i=1) {
    if(i == 1) {
      if(y >= 0) {
        return(log(y + 1)/log(2))
      } else {
        return((1 - (1 - y)^2)/(2*log(2)))
      }
    } 
  }
  
  f2 <- function(z, i=1) {
    if(i == 1) {
      if(z >= 0) {
        return(exp(log(2)*z) - 1)
      } else {
        return(1 - sqrt(1 - 2*log(2)*z))
      }
    }
  }
  
  X <- runif(n)
  noise <- rnorm(n)
  noise[noise < -3] <- -3
  noise[noise > 3] <- 3
  
  Z <- 6*X - 3 + noise
  Y <- sapply(Z, f2)
  
  res <- list("X"=as.matrix(X), "Y"=Y, "noise"=noise, "Z"=Z)
  return(res)
}

ep.kernel <- function(u) {
  3*(abs(u) <= 1)*(1 - u^2)/4
}

int.ep.kernel <- function(u) {
  # if(u <= -1) {
  #   return(0)
  # } else if(u >= 1) {
  #   return(1)
  # } else {
  #   return((3*u - u^3 + 2)/4)
  # }
  return(ifelse(u >= 1, 1, 1*(u > -1)*(3*u - u^3 + 2)/4))
}

grad.ep.kernel <- function(u) {
  -3*(abs(u) <= 1)*u/2
}

# plot(seq(-3, 3, 0.001), sapply(seq(-3, 3, 0.001), grad.ep.kernel))

classical.bandwidth <- function(X) {
  ((40*sqrt(pi))^(1/5))*sd(X)/(length(X)^(1/5))
}

get.U.hat <- function(Y) {
  n <- length(Y)
  F_hat_Y <- rank(Y)/n
  U_hat <- (F_hat_Y - sum(Y <= 0)/n)/(sum(Y <= 1)/n - sum(Y <= 0)/n)
  return(U_hat)
}

T.hat.y <- function(y, Y) {
  ui <- (sum(Y <= y) - sum(Y <= 0))/(sum(Y <= 1) - sum(Y <= 0))
  return(ui)
}

# approx. 
p.hat.ux <- function(u, x, U_hat, X, K, int.K) {
  hx <- classical.bandwidth(X)
  hu <- classical.bandwidth(U_hat)
  
  n <- nrow(X)
  s <- 0
  for(i in 1:n) {
    s <- s + int.K((u - U_hat[i])/hu)*prod(sapply((X[i]-x)/hx, K)/hx)
  }
  return(s/n)
}

# approx. density of X
fX.hat.x <- function(x, X, K) {
  hx <- classical.bandwidth(c(X))
  n <- nrow(X)
  s <- 0
  for(i in 1:n) {
    s <- s + prod(sapply((X[i]-x)/hx, K)/hx)
  }
  return(s/n)
}

# approx. conditional distribution function of U given X
phi.hat.ux <- function(u, x, U_hat, X, K, int.K) {
  p_hat_ux <- p.hat.ux(u, x, U_hat, X, K, int.K)
  fX_hat_x <- fX.hat.x(x, X, K)
  return(p_hat_ux/fX_hat_x)
}

grad.u.p.hat.ux <- function(u, x, U_hat, X, K) {
  hx <- classical.bandwidth(X)
  hu <- classical.bandwidth(U_hat)
  
  n <- nrow(X)
  grad <- 0
  for(i in 1:n) {
    grad <- grad + K((u - U_hat[i])/hu)*prod(sapply((X[i]-x)/hx, K)/hx)/hu
  }
  return(grad/n)
} 

grad.j.p.hat.ux <- function(u, x, U_hat, X, K, int.K, grad.K, j) {
  hx <- classical.bandwidth(X)
  hu <- classical.bandwidth(U_hat)
  
  n <- nrow(X)
  grad <- 0
  for(i in 1:n) {
    A <- int.K((u - U_hat[i])/hu)
    B <- 1
    for(l in 1:ncol(X)) {
      if(l == j) {
        B <- -B*grad.K((X[i][j]-x[j])/hx)/(hx^2)
      } else {
        B <- B*K((X[i][l]-x[l])/hx)/hx
      }
    }
    grad <- grad + A*B
  }
  return(grad/n)
}

grad.j.fX.hat.x <- function(x, X, K, grad.K, j) {
  hx <- classical.bandwidth(c(X))
  n <- nrow(X)
  grad <- 0
  for(i in 1:n) {
    B <- 1
    for(l in 1:ncol(X)) {
      if(l == j) {
        B <- -B*grad.K((X[i][j]-x[j])/hx)/(hx^2)
      } else {
        B <- B*K((X[i][l]-x[l])/hx)/hx
      }
    }
    grad <- grad + B
  }
  return(grad/n)
}

grad.u.phi.hat.ux <- function(u, x, U_hat, X, K) {
  grad_u_p_hat_ux <- grad.u.p.hat.ux(u, x, U_hat, X, K)
  fX_hat_x <- fX.hat.x(x, X, K)
  return(grad_u_p_hat_ux/fX_hat_x)
}

grad.j.phi.hat.ux <- function(u, x, U_hat, X, K, int.K, grad.K, j) {
  fX_hat_x <- fX.hat.x(x, X, K)
  grad_j_fX_hat_x <- grad.j.fX.hat.x(x, X, K, grad.K, j)
  p_hat_ux <- p.hat.ux(u, x, U_hat, X, K, int.K)
  grad_j_p_hat_ux <- grad.j.p.hat.ux(u, x, U_hat, X, K, int.K, grad.K, j)
  
  numerator <- grad_j_p_hat_ux*fX_hat_x - p_hat_ux*grad_j_fX_hat_x
  denom <- fX_hat_x^2
  
  return(numerator/denom)
}

s.j.hat.ux <- function(u, x, U_hat, X, K, int.K, grad.K, j) {
  grad_u_phi_hat_ux <- grad.u.phi.hat.ux(u, x, U_hat, X, K)
  grad_j_phi_hat_ux <- grad.j.phi.hat.ux(u, x, U_hat, X, K, int.K, grad.K, j)
  s <- grad_u_phi_hat_ux/grad_j_phi_hat_ux
  return(s)
}

S.j.hat.ux <- function(u, x, U_hat, X, K, int.K, grad.K, j) {
  res <- tryCatch(
    {
      integrate(s.j.hat.ux, lower=0, upper=u, x=x, U_hat=U_hat, X=X, K=ep.kernel, int.K=int.ep.kernel, grad.K=grad.ep.kernel, j)$value
    },
    error=function(cond) {
      message(cond)
      print("error")
      return(NA)
    }
  )
  return(res)
}

lambda.j.hat.ux <- function(u, x, U_hat, X, K, int.K, grad.K, j) {
  numerator <- S.j.hat.ux(u, x, U_hat, X, K, int.K, grad.K, j)
  denom <- S.j.hat.ux(1, x, U_hat, X, K, int.K, grad.K, j)
  if(is.na(numerator) | is.na(denom) | denom == 0) {
    return(NA)
  }
  return(numerator/denom)
}

# FIXME choose proper arguments for Nx and kernels
lambda.vals <- function(u, U_hat, X, K=ep.kernel, int.Kint.ep.kernel, 
                        grad.K=grad.ep.kernel, j=1, Nx=100) {
  points <- seq(min(X), max(X), (max(X) - min(X))/Nx)
  
  lambda_vals <- c()
  point_vals <- c()
  for(point in points) {
    val <- lambda.j.hat.ux(u, point, U_hat, X, K, int.K, grad.K, j)
    if(is.na(val)) {
      next
    }
    lambda_vals <- c(lambda_vals, val)
    point_vals <- c(point_vals, point)
  }
  
  return(lambda_vals)
}

f2.inv.est.LS <- function(y, Y, X) {
  U_hat <- get.U.hat(Y)
  res <- lambda.vals(T.hat.y(y, Y), U_hat, X)
  return(mean(res))
}

f2.inv.est.LAD <- function(y, Y, X) {
  U_hat <- get.U.hat(Y)
  res <- lambda.vals(T.hat.y(y, Y), U_hat, X)
  return(median(res))
}


# Nadaraya–Watson estimator for f_1 after estimating f_2
f1.est.NW <- function(x, f2_inv_Y, X, K) {
  hx <- classical.bandwidth(c(X))
  
  W <- sapply((x - X)/hx, K)/hx
  W <- W/sum(W)
  f1_x <- sum(W*Y)
  return(f1_x)
}

noise.est.gtm <- function(Y, X) {
  dumy.f1 <- function(x) {
    f1.est.NW(x, est_Z, X, ep.kernel)
  }
  
  dummy.f2.inv <- function(y) {
    est_Z <- f2.inv.est.LS(y, Y, X)
  }
  est_Z <- sapply(Y, dummy.f2.inv)
  est_noise <- est_Z - sapply(X ,dumy_f1)
  
  return(est_noise)
}

# data <- simulate.gtm.data(1000)
# X <- data$X
# Y <- data$Y
# 
# X <- as.matrix(X)
# U_hat <- get.U.hat(Y)
# 
# U_hat
# 
# library(ggplot2)

# p.hat.ux(0.001, min(X), U_hat, X, ep.kernel, int.ep.kernel)
# dummy <- function(l, f=max) {
#   p.hat.ux(l, f(X), U_hat, X, ep.kernel, int.ep.kernel)
# }
# inds <- seq(0, 1, 0.01)
# vals <- sapply(inds, dummy)
# ggplot() + geom_line(aes(x=inds, y=vals), colour="red")


# fX.hat.x(median(X), X, ep.kernel)
# dummy <- function(l) {
#   fX.hat.x(l, X, ep.kernel)
# }
# inds <- seq(min(X), max(X), 0.01)
# vals <- sapply(inds, dummy)
# ggplot() + geom_line(aes(x=inds, y=vals), colour="red")


# phi.hat.ux(0.001, max(X), U_hat, X, ep.kernel, int.ep.kernel)
# dummy <- function(l, f=min) {
#   phi.hat.ux(l, f(X), U_hat, X, ep.kernel, int.ep.kernel)
# }
# inds <- seq(0, 1, 0.01)
# vals <- sapply(inds, dummy)
# ggplot() + geom_line(aes(x=inds, y=vals), colour="red")

# grad.u.p.hat.ux(0.01, min(X), U_hat, X, ep.kernel)
# dummy <- function(l, f=min) {
#   grad.u.p.hat.ux(l, f(X), U_hat, X, ep.kernel)
# }
# inds <- seq(0, 1, 0.01)
# vals <- sapply(inds, dummy)
# ggplot() + geom_line(aes(x=inds, y=vals), colour="red")

# grad.j.p.hat.ux(0.01, max(X), U_hat, X, ep.kernel, int.ep.kernel, grad.ep.kernel, 1)
# dummy <- function(l, f=mean) {
#   grad.j.p.hat.ux(l, f(X), U_hat, X, ep.kernel, int.ep.kernel, grad.ep.kernel, 1)
# }
# inds <- seq(0, 1, 0.01)
# vals <- sapply(inds, dummy)
# ggplot() + geom_line(aes(x=inds, y=vals), colour="red")

# grad.i.fX.hat.x(max(X), X, ep.kernel, grad.ep.kernel, 1)
# dummy <- function(l) {
#   grad.i.fX.hat.x(l, X, ep.kernel, grad.ep.kernel, 1)
# }
# inds <- seq(min(X), max(X), 0.01)
# vals <- sapply(inds, dummy)
# ggplot() + geom_line(aes(x=inds, y=vals), colour="red")

# grad.u.phi.hat.ux(0.8, min(X), U_hat, X, ep.kernel)
# dummy <- function(l, f=min) {
#   grad.u.phi.hat.ux(l, f(X), U_hat, X, ep.kernel)
# }
# inds <- seq(0, 1, 0.01)
# vals <- sapply(inds, dummy)
# ggplot() + geom_line(aes(x=inds, y=vals), colour="red")

# grad.j.phi.hat.ux(0.01, min(X), U_hat, X, ep.kernel, int.ep.kernel, grad.ep.kernel, 1)
# dummy <- function(l, f=min) {
#   grad.j.phi.hat.ux(l, f(X), U_hat, X, ep.kernel, int.ep.kernel, grad.ep.kernel, 1)
# }
# inds <- seq(0, 1, 0.01)
# vals <- sapply(inds, dummy)
# ggplot() + geom_line(aes(x=inds, y=vals), colour="red")


# s.j.hat.ux(0.8, max(X), U_hat, X, ep.kernel, int.ep.kernel, grad.ep.kernel, 1)
# dummy <- function(l, f=max) {
#   s.j.hat.ux(l, f(X), U_hat, X, ep.kernel, int.ep.kernel, grad.ep.kernel, 1)
# }
# inds <- seq(0, 1, 0.01)
# vals <- sapply(inds, dummy)
# ggplot() + geom_line(aes(x=inds, y=vals), colour="red")

# f2_inv <- function(y, i=1) {
#   if(i == 1) {
#     if(y >= 0) {
#       return(log(y + 1)/log(2))
#     } else {
#       return((1 - (1 - y)^2)/(2*log(2)))
#     }
#   }
# }

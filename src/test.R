# library(dHSIC)
# 
# # non identifiability
# sample_size <- 10000
# beta <- runif(1, min = -100, max = 100)
# x <- rnorm(sample_size)
# noise_2 <- rnorm(sample_size)
# y = sign(beta*x + noise_2)*abs((beta*x + noise_2))^(1/3)
# 
# dhsic(x, noise_2, kernel = "gaussian.fixed")$dHSIC
# cor(x, noise_2)
# 
# noise_1 <- (x - beta*noise_2)/(beta^2+1)**0.5
# dhsic(y^3, noise_1, kernel = "gaussian.fixed")$dHSIC
# cor(y, noise_1)
# 
# plot(x, y)
# plot(y, x)

args <- commandArgs(trailingOnly = T)
print(args[1])
print(args[2])



?mgcv::b.spline

source("utils.R")
library(splines2)

data <- simulate_rank_regression_data(5, 2)
X <- data$X
Y <- data$Y

data$beta
dim(X)
g_x1 <- bSpline(X[, 1], intercept=T, degree=10)
dg_x1 <- deriv(g_x1)
g_x2 <- bSpline(X[, 2], intercept=T)
dg_x2 <- deriv(g_x2)

predict(g_x1, X[, 1][1])
g_x1

predict(g_x2, X[, 2][1])
g_x2

?bSpline
dg_x1
X[, 1]

ml_optim <- function(YY, XX, dXX) {
  dim_x <- dim(XX)[2]
  dim_y <- dim(YY)[2]
  
  neg_lik <- function(beta, YY, XX, dXX) {
    dim_y <- dim(YY)[2]
    
    beta_y <- beta[1:dim_y]
    beta_x <- beta[(dim_y+1):(length(beta))]
    res <- sum((YY %*% beta_y - XX %*% beta_x)^2)/2 - sum(log(abs(dXX %*% beta_x)))
    return(res)
  }
  
  # first part for YY second part for XX
  coefs <- runif(dim_y+dim_x, -10, 10)
  print(coefs)
  
  res <- optim(coefs, fn=neg_lik, method = "BFGS", 
               control = list(trace=T, REPORT=1),
               YY=YY, XX=XX, dXX=dXX)
  
  est <- res$par
  est_y <- est[1:dim_y]
  est_x <- est[(dim_y+1):(length(est))]
  h_y <- YY %*% est_y
  g_x <- XX %*% est_x
  noise <- h_y - g_x
  return(list("est_beta"=est, "h_y"=h_y, "g_x"=g_x, "noise"=noise))
}

data <- simulate_rank_regression_data(1000, 1)
X <- data$X
Y <- data$Y

XX <- bSpline(X, df=6, intercept=F)
dXX <- deriv(XX)
YY <- bSpline(Y, df=7, intercept=T)

res <- ml_optim(YY, XX, dXX)
res$est_beta
plot(data$noise, res$noise)

plot(data$noise, -res$noise)

library(ggplot2)
f <- function(x) (x - 4.7)^3
xs <- seq(min(data$Y)-0.5, max(data$Y)+0.5, 0.1)
ggplot() + geom_point(aes(x = data$Y, y=res$h_y, colour="estimated h")) + 
  geom_line(aes(x = xs, y=sapply(xs, f), colour="true h")) +
  scale_color_manual(name='', breaks = c("estimated h", "true h"), values = c("black", "red")) +
  labs(x='y', y='h(y)') +
  theme(legend.position = 'top')

# XX <- bSpline(X, intercept=F)
# dXX <- deriv(XX)
# YY <- bSpline(Y, intercept=T)
# 
# res <- ml_optim(cbind(Y, YY), cbind(X, XX), cbind(1, dXX))
# 
# plot(data$noise, res$noise)
# 
# plot(data$noise, -res$noise)
# 
# mean((data$noise - res$noise)^2)
# mean((data$noise + res$noise)^2)
# 
# 
# 
# XX <- naturalSpline(X, intercept=F)
# dXX <- deriv(XX)
# YY <- naturalSpline(Y, intercept=T)
# 
# res <- ml_optim(cbind(Y, YY), cbind(X, XX), cbind(1, dXX))
# 
# plot(data$noise, res$noise)
# 
# plot(data$noise, -res$noise)
# 
# mean((data$noise - res$noise)^2)
# mean((data$noise + res$noise)^2)

# simple polynomials, maybe ortogonalization will improve
data <- simulate_rank_regression_data(1000, 1)
X <- data$X
Y <- data$Y

XX <- cbind(X, X^2, X^3)
dXX <- cbind(1, 2*X, 3*X^2)
YY <- cbind(1, Y, Y^2, Y^3, Y^4)
head(XX)
head(dXX)
res <- ml_optim(YY, XX, dXX)
res$est_beta
data$beta
plot(data$noise, res$noise)

plot(data$noise, -res$noise)

f <- function(x) (x - 4.7)^3
xs <- seq(min(data$Y)-0.5, max(data$Y)+0.5, 0.1)
ggplot() + geom_point(aes(x = data$Y, y=res$h_y, colour="estimated h")) + 
  geom_line(aes(x = xs, y=sapply(xs, f), colour="true h")) +
  scale_color_manual(name='', breaks = c("estimated h", "true h"), values = c("black", "red")) +
  labs(x='y', y='h(y)') +
  theme(legend.position = 'top')

ggplot() + geom_point(aes(x = data$Y, y=-res$h_y, colour="estimated h")) + 
  geom_line(aes(x = xs, y=sapply(xs, f), colour="true h")) +
  scale_color_manual(name='', breaks = c("estimated h", "true h"), values = c("black", "red")) +
  labs(x='y', y='h(y)') +
  theme(legend.position = 'top')

res$est_beta

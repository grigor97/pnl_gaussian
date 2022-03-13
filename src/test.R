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


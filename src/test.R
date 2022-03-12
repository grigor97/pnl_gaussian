sample_size <- 1000
beta <- 4.3 
x <- rnorm(sample_size)
noise <- rnorm(sample_size)
y = sign(beta*x + noise)*abs((beta*x + noise))^(1/3)

plot(x, y)
plot(y, x)

# shifted lognormal analysis

n <- 100000000
alpha <- 0
beta <- 1
RF <- 0
MKT <- 0
min.mkt <- 0 # investigate why I did this!
stdv <- 0.3

# simple linear
normal <- alpha + RF + beta * MKT + rnorm(n, 0, stdv)

# shifted lognormal (goal: same mean and standard deviation for normal & shifted.lognormal)
lower.bound <- -1 # set to a total return of -100%
sigma <- sqrt(log(1 + stdv^2 / lower.bound^2))
mu <- log(-lower.bound) - sigma^2 / 2
shifted.lognormal <- alpha + RF + beta * MKT + exp(rnorm(n, mu, sigma)) + lower.bound

mean(normal) - mean(shifted.lognormal) # = 0
sd(normal) - sd(shifted.lognormal) # = 0

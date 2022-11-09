### Simple Averaged Pricing Error Example
# Prologue ----

if(sys.nframe() == 0L) rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Model -----

sim.fun <- function(input.data=NA) {
  
  beta <- 1 # paramter we want to estimate
  n <- 5000 # number of cross-sectional units
  periods <- 1 # number of time-periods for discounting
  
  sigma.sdf <- 0.1
  sigma.err <- 0.3
  
  # Fill df
  df <- data.frame(sdf = rep(1, n))
  for (i in 1:periods) {
    # public market factor
    # x <- rnorm(n, 0.08, sigma.sdf)
    x <- rlnorm(n, 0.08-sigma.sdf^2/2, sigma.sdf) - 1
    
    df[, paste0("x", i)] <- x
    # error term
    error <- "multiplicative"
    if (error == "additive") {
      e <- rnorm(n, 0, sigma.err)
      y <- 1 / (1 + beta * x + e)
    } else {
      e <- rlnorm(n, 0-sigma.err^2/2, sigma.err)
      
      # error in denominator
      # y <- 1 / ((1 + beta * x)  * e)
      # error in numerator (yields worse estimation results)
      y <- e / (1 + beta * x)
    }

    df$sdf <- df$sdf * y
  }
  
  cf <- 1 / df$sdf
  
  if(periods > 1) {
    # Define loss functions
    loss.npv <- function(b) {
      df.x <- 1 + df[, -1] * b
      sdf <- 1 / apply(df.x, 1, prod)
      npv <- cf * sdf - 1
      error <- sum((npv)^2)
      return(error)
    }
    
    loss.fv <- function(b) {
      df.x <- 1 + df[, -1] * b
      sdf <- 1 / apply(df.x, 1, prod)
      fv <- cf - 1 / sdf
      error <- sum((fv)^2)
      return(error)
    }
    
    loss.averaged <- function(b) {
      df.x <- 1 + df[, -1] * b
      sdf <- 1 / apply(df.x, 1, prod)
      npv <- cf * sdf - 1
      fv <- cf- 1 / sdf
      error <- sum((npv)^2) + sum((fv)^2)
      return(error)
    }
  } else {
    # Define loss functions
    loss.npv <- function(b) {
      df.x <- 1 + df[, -1] * b
      sdf <- 1 / df.x
      
      # E[M/Me] (prices the cash flow by the sdf)
      #npv <-  1 - cf * sdf
      npv <- cf*sdf - 1
      
      # E[Me/M]
      # npv <- 1 / (sdf * cf) - 1
      
      # both
      #npv <- 1 / (sdf * cf) + cf*sdf - 2
      
      error <- sum((npv)^2)
      return(error)
    }
    
    loss.fv <- function(b) {
      df.x <- 1 + df[, -1] * b
      sdf <- 1 / df.x
      
      # good (replicates the data generating process)
      #fv <- cf - 1 / sdf
      fv <- 1 / sdf - cf
      
      # bad
      # fv <- 1 / cf - sdf
      
      error <- sum((fv)^2)
      return(error)
    }
    
    loss.averaged <- function(b) {
      df.x <- 1 + df[, -1] * b
      sdf <- 1 / df.x
      npv <- 1 - cf * sdf
      fv <- cf - 1 / sdf
      error <- sum((npv)^2) + sum((fv)^2)
      return(error)
    }
  }
  

  
  # Optimize
  res <- optimise(loss.npv, 0:10)
  df.res <- data.frame(NPV = res$minimum)
  res <- optimise(loss.fv, 0:10)
  df.res$FV <- res$minimum
  res <- optimise(loss.averaged, 0:10)
  df.res$AVG <- res$minimum
  
  return(df.res)
}
sim.fun()

set.seed(55)
system.time(
  df.out <- lapply(1:10000, sim.fun)
)
df.out <- data.frame(do.call(rbind, t(df.out)))
colMeans(df.out)
apply(df.out, 2, sd)
summary(df.out)

hist(df.out$NPV)
hist(df.out$FV)

'
> colMeans(df.out)
     NPV       FV      AVG 
2.927163 2.701741 2.757635 
> apply(df.out, 2, sd)
NPV        FV       AVG 
0.5686698 0.8367402 0.5904145 
'
# Multiplicative example ----
n <- 1000 * 1000 * 10
s1 <- 1
sdf <- rlnorm(n, -0.1 - s1^2/2, s1)
mean(sdf)
mean(1/sdf)
cor(sdf, 1/sdf)

s2 <- 1
e <- rlnorm(n, 0 - s2^2/2, s2)
mean(e)
cor(sdf*e, 1/(sdf*e))

summary(sdf/(sdf*e))
summary((sdf*e)/sdf)
summary(1/sdf - 1/(sdf*e))

x <- (sdf * e) / sdf
mean(x)
sd(x)
cor(sdf, x)
hist(x)
hist(1/sdf - 1/(sdf*e))
# s1=1, s2=1: -0.191
# s1=1, s2=2: -0.01226564
# s1=2, s2=1: -0.04248953
# s1=0.1, s2=2:-0.01795819
# s1=0.1, s2=0.1: -0.6983986
# s1=0.001, s2=0.001: -0.7070273
# s1=0.0001, s2=3: -0.0003168856
# s1=1, s2=0.01: -0.00765646


# Another multiplicative example -----

n <- 1000000

m1 <- rlnorm(n, -0.1,0.2)
m2 <- rlnorm(n, -0.1,0.2)
mean(m1)

s1 <- 0.1
e1 <- rlnorm(n, 0 - s1^2/2, s1)
e2 <- rlnorm(n, 0 - s1^2/2, s1)
e3 <- rlnorm(n, 0 - s1^2/2, s1)
mean(c(e1, e2, e3))
mean((e1)^2)

cf1 <- 1 / m1 * e1
cf2 <- 1 / (m1 * m2) * e2 * e3

cov(m1, cf1 + m2 * cf2)
cov(m1 * m2, cf1 /m2 + cf2)

r <- rlnorm(n, 0.1, 0.2)
summary(r)
mean(r) * mean(1/r) + cor(1/r, r) * sd(1/r) * sd(r)

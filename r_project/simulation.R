# simulation study ----
if(sys.nframe() == 0L) rm(list = ls())
library(data.table)
library(parallel)
library(doParallel)
library(foreach)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(100)

df.q5 <- read.csv("data_prepared/q_factors.csv")
df.q5$Date <- as.Date(df.q5$Date)

sim.sd <- function(stdv = 0.312) {
  start <- as.Date("1985-12-31")
  end <- as.Date("2005-12-31")
  mkt <- df.q5$MKT[(df.q5$Date >= start) & (df.q5$Date <= end)]
  rf <- df.q5$RF[(df.q5$Date >= start) & (df.q5$Date <= end)]
  print(paste("min MKT:", min(mkt), "max RF", max(rf)))
  y <- replicate(1000, mkt + rnorm(length(mkt), 0, stdv))
  print(paste(mean(y) * 12, sd(y) * sqrt(12)))
}
sim.sd()


# sim data -----
make.fund <- function(no=0, vintage=1990, stdv = 0.01, exp.aff = FALSE, df=df.q5) {
  # parameter
  N <- 15
  T <- 10 # 5
  months <- 12
  alpha <- -0.0025 # per month
  beta <- 2.5
  
  # market
  df <- df[as.integer(format(df$Date, "%Y")) >= vintage, ]
  df$Vintage <- vintage
  df$CF <- 0
  sdf <- ifelse(exp.aff, "ea", "sl")
  df$type <- paste0("sdf:", sdf, " stdv:", stdv, " #d:", N)
  df$Fund.ID <- paste(vintage, no, sep = "_")
  
  # shifted lognormal (mu, sigma)
  min.mkt <- 0.25
  foo <- function(sigma) {
    mu <- log((1-min.mkt) / exp(0.5*sigma^2))
    y <- exp(mu+0.5*sigma^2) * sqrt(exp(sigma^2)-1)
    return((stdv - y)^2)
  }
  sigma <- optimize(foo, c(0,10))$minimum
  mu <- log((1-min.mkt) / exp(0.5*sigma^2))
  
  # simulate N deals
  for (i in 1:N) {
    start <- sample(1:T/2 * months, 1)
    end <- start + sample(1:T * months, 1)
    end <- min(end, nrow(df))
    if(exp.aff) {
      df$deal <- alpha + log(1 + df$RF) + beta * log(1 + df$MKT) + rnorm(nrow(df), 0, stdv) - stdv^2/2
      m <- exp(sum(df$deal[start:end])) # exp affine
    } else {
      # simple linear
      df$deal <- alpha + df$RF + beta * df$MKT + rnorm(nrow(df), 0, stdv)
      #df$deal <- alpha + df$RF + beta * df$MKT + exp(rnorm(nrow(df), 0, stdv) - stdv^2/2) - 1
      
      # shifted lognormal
      #df$deal <- alpha + df$RF + beta * df$MKT + exp(rnorm(nrow(df), mu, sigma)) - (1-min.mkt)
      
      path <- cumprod(1 + df$deal[start:end])
      if(any(path < 0)) { # default
        m <- 0
      } else { # no default
        m <- tail(path, 1)
      }
    }
    df$CF[start] <- df$CF[start] - 1
    df$CF[end] <- df$CF[end] + m
  }
  
  df <- df[df$Date <= max(df$Date[df$CF != 0]), ]
  df <- df[df$CF != 0, ]
  cols <- c("Date", "CF", "type", "Vintage", "Fund.ID")
  df <- df[, cols]
  return(df)
}
# View(make.fund())

cl <- parallel::makeForkCluster(2)
doParallel::registerDoParallel(cl)
output <- foreach (x = 1:1000, .verbose = FALSE) %do% {
  print(x)
  set.seed(x)
  l <- list()
  no.funds <- 20 # 20
  min.vin <- 1986  # 1986 # 1967 # 1996
  max.vin <- 2005 # 2005 # 1995 
  stdvs <- 0.2 # 0.2 # 0.312
  for(i in 1:no.funds) {
    for(v in min.vin:max.vin) {
      for(s in stdvs) {
        l[[paste(i, v, s)]] <- make.fund(i, v, s)
      }
    }
  }
  df.out <- as.data.frame(data.table::rbindlist(l))
  df.out$type <- paste0(df.out$type, " #f:", no.funds, " #y:", max.vin - min.vin, " .", x)
  df.out$Fund.ID <- paste0(df.out$type, "+", df.out$Fund.ID)
  #table(df.out$Vintage[!(duplicated(df.out$Fund.ID))])
  return(df.out)
}
parallel::stopCluster(cl)
system.time(df.out <- data.table::rbindlist(output))
rm(output)

length(levels(as.factor(df.out$Fund.ID)))
length(1986:2005) * 20 * length(1:100)

system.time(df.vyp <- df.out[, .(CF=sum(CF)), by = .(Date, type, Vintage)])
df.vyp$Fund.ID <- paste0(df.vyp$type, " v:", df.vyp$Vintage)
head(df.vyp)

system.time(write.csv(df.vyp, "data_prepared/simulated_cashflows_EW_VYP.csv", row.names=FALSE))
#system.time(write.csv(df.out, "data_prepared/simulated_cashflows_EW.csv", row.names=FALSE))

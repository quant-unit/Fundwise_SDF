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


# sim data function -----

if (FALSE) {
  # parameter
  no.deals <- 15
  investment.period <- 5 # (in years)
  max.holding.period <- 10 # 5 # (in years)
  alpha <- -0.0025 # per month
  beta <- 2.5
  
  # parameter 2
  no.samples <- 1000
  no.funds <- 20 # 20
  min.vin <- 1986  # 1986 # 1967 # 1996
  max.vin <- 2005 # 2005 # 1995 
  stdvs <- 0.2 # 0.2 # 0.312
  exp.aff.sdf <- FALSE
}

create.simulation <- function(
    no.deals,
    investment.period,
    max.holding.period,
    alpha,
    beta,
    no.samples,
    no.funds,
    min.vin,
    max.vin,
    stdvs,
    exp.aff.sdf
) {
  
  months <- 12 # constant (no parameter to vary)
  
  make.fund <- function(no=0, vintage=1990, stdv = 0.01, exp.aff = exp.aff.sdf, df=df.q5) {
    # market
    df <- df[as.integer(format(df$Date, "%Y")) >= vintage, ]
    df$Vintage <- vintage
    df$CF <- 0
    sdf <- ifelse(exp.aff, "ea", "sl")
    df$type <- paste0("sdf:", sdf, " stdv:", stdv, " #d:", no.deals)
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
    for (i in 1:no.deals) {
      if (investment.period == 0) {
        start <- 1
      } else {
        start <- sample(1:(investment.period * months), 1)
      }
      end <- start + sample(1:(max.holding.period * months), 1)
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

  multi <- FALSE
  
  if (multi) {
    
    cl <- parallel::makeForkCluster(12)
    doParallel::registerDoParallel(cl)
    output <- foreach (x = 1:no.samples, .verbose = FALSE) %do% {
      print(x)
      set.seed(x)
      l <- list()
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
    
  } else {

    output <- list()
    for (x in 1:no.samples) {
      #print(x)
      set.seed(x)
      l <- list()
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
      output[[paste0("i", x)]] <- df.out
    }
    
  }

  system.time(df.out <- data.table::rbindlist(output))
  rm(output)
  
  length(levels(as.factor(df.out$Fund.ID)))
  length(min.vin:max.vin) * 20 * length(1:100)
  
  system.time(df.vyp <- df.out[, .(CF=sum(CF)), by = .(Date, type, Vintage)])
  df.vyp$Fund.ID <- paste0(df.vyp$type, " v:", df.vyp$Vintage)
  head(df.vyp)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename.vyp <- paste0("data_prepared_sim/", timestamp,"_simulated_cashflows_EW_VYP.csv")
  filename.fund <- paste0("data_prepared_sim/", timestamp,"_simulated_cashflows_EW.csv")
  filename.meta <- paste0("data_prepared_sim/", timestamp,"_simulated_cashflows_EW_meta.csv")
  
  df.meta <- data.frame(
    investment.period=investment.period,
    max.holding.period=max.holding.period,
    no.deals=no.deals,
    alpha=alpha,
    beta=beta,
    no.samples=no.samples,
    no.funds=no.funds,
    min.vin=min.vin,
    max.vin=max.vin,
    stdvs=stdvs,
    exp.aff.sdf=exp.aff.sdf
  )
  df.meta
  
  system.time(write.csv(df.vyp, filename.vyp, row.names=FALSE))
  system.time(write.csv(df.out, filename.fund, row.names=FALSE))
  system.time(write.csv(df.meta, filename.meta, row.names=FALSE))
  
}

# sim runs ----

## base case: cross-sectional unit
system.time(
  create.simulation(
  no.deals=15,
  investment.period=5,
  max.holding.period=10,
  alpha=0,
  beta=1,
  no.samples=1000,
  no.funds=20,
  min.vin=1986,
  max.vin=2005,
  stdvs=0.2,
  exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=5 + max.holding=5)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=5,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## High Beta, Negative Alpha
system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=10,
    alpha=-0.0025, # per month
    beta=2.5,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)


## Exp aff models

### exp.aff, beta=1
system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=10,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=TRUE
  )
)


### exp.aff, negative alpha, high beta
system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=10,
    alpha=-0.0025, # per month
    beta=2.5,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=TRUE
  )
)


## Double half models

### Big n/V: 40 Funds, vintages 1986-2005
system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=10,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=40,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)


### Big V: 10 Funds, vintages 1967-2005
system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=10,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=10,
    min.vin=1967,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)


### Big V: 20 Funds, vintages 1967-2005


system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=10,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1967,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)


### Small V: vintages 1986-1995
system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=10,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=1995,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)


### Small V: vintages 1996-2005
system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=10,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1996,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=2 + max.holding=2)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=2,
    max.holding.period=2,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=4 + max.holding=4)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=4,
    max.holding.period=4,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=2 + max.holding=4)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=2,
    max.holding.period=4,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=4 + max.holding=2)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=4,
    max.holding.period=2,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=1 + max.holding=5)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=1,
    max.holding.period=5,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=5 + max.holding=1)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=5,
    max.holding.period=1,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=3 + max.holding=3)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=3,
    max.holding.period=3,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=1 + max.holding=1)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=1,
    max.holding.period=1,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=6 + max.holding=6)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=6,
    max.holding.period=6,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=7 + max.holding=7)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=7,
    max.holding.period=7,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=8 + max.holding=8)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=8,
    max.holding.period=8,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=9 + max.holding=9)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=9,
    max.holding.period=9,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=0 + max.holding=1)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=0,
    max.holding.period=1,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=0 + max.holding=3)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=0,
    max.holding.period=3,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=0 + max.holding=5)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=0,
    max.holding.period=5,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=0 + max.holding=7)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=0,
    max.holding.period=7,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)

## Shorter fund lifetime max(investing=0 + max.holding=10)
system.time(
  create.simulation(
    no.deals=15,
    investment.period=0,
    max.holding.period=10,
    alpha=0,
    beta=1,
    no.samples=1000,
    no.funds=20,
    min.vin=1986,
    max.vin=2005,
    stdvs=0.2,
    exp.aff.sdf=FALSE
  )
)
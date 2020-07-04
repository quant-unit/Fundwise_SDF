#### estimated model
# 0) Prologue -----
if(FALSE) {
  files <- c("prep_preqin.R", "prep_public.R")
  for(file in files) {
    source(file)
    if(sys.nframe() == 0L) rm(list = ls())
  }
}

if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 1.1) PARAMETERS ----
weighting <- "EW"
use.vintage.year.pfs <- TRUE
use.simulation <- FALSE
do.cross.validation <- TRUE
do.cache <- TRUE
do.parallel <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")

public.filename <- "public_returns"
public.filename <- "msci_market_factors"
public.filename <- "q_factors"

error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

max.months <- 180 # c(12.5, 15, 17.5) * 12 # c(1/12, 5 , 10, 15, 20, 25, 30) * 12
lambdas <- 0
kernel.bandwidth <- 12

# 1.2) load data -----
df.public <- read.csv(paste0("data_prepared/", public.filename, ".csv"))
df.public$Date <- as.Date(df.public$Date)

if(!use.simulation) {
  df.preqin0 <- read.csv(paste0("data_prepared/preqin_cashflows_", weighting, ".csv"))
} else {
  df.preqin0 <- read.csv(paste0("data_prepared/simulated_cashflows_", weighting, ".csv"))
}
df.preqin0$Date <- as.Date(df.preqin0$Date)

to.monthly <- function(df.ss) {
  # fill zero cash flows (necessary for estimation)
  max.date <- max(df.ss$Date) + 1
  #max.date <- as.Date("2020-01-01")
  df.m <- data.frame(Date = seq(min(df.ss$Date) + 1, max.date, by = "month") - 1)
  df.ss <- merge(df.ss, df.m, by = "Date", all = TRUE)
  df.ss$CF[is.na(df.ss$CF)] <- 0
  for(col in c("type", "Vintage", "Fund.ID")) df.ss[, col] <- df.ss[1, col]
  return(df.ss)
}

df.preqin <- as.data.frame(data.table::rbindlist(lapply(split(df.preqin0, df.preqin0$Fund.ID), to.monthly)))
rm(df.preqin0)
df.preqin$type <- as.factor(as.character(df.preqin$type))
df.preqin$Fund.ID <- as.factor(as.character(df.preqin$Fund.ID))
length(levels(df.preqin$type))

df.ss <- df.preqin[df.preqin$Fund.ID == levels(df.preqin$Fund.ID)[5], ]
#View(df.ss)
rm(df.ss)

df0 <- merge(df.preqin, df.public, by="Date")
df0$Fund.ID <- as.factor(df0$Fund.ID)
df0$Alpha <- 1
rm(df.preqin)

#df0 <- df0[df0$Vintage <= 2011, ]
aggregate(Vintage ~ type , df0, min)
no.of.vin <- function(x) length(table(x))
aggregate(Vintage ~ type , df0, no.of.vin)

# 1.3) rbind.all ----
rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}

# 2.1) getNPVs function ----
library(Rcpp)

Rcpp::cppFunction("
    double getNPVs(NumericVector cf, NumericVector ret, int max_month){
    double DCF = sum(cf / ret);

    int n = cf.length();
    int i_max = std::min(max_month, n);
    NumericVector VectorOut(i_max, 0.0);
    
    for(int i = 0; i < i_max; i++) {
      VectorOut[i] = DCF * ret[i];
    }

    double y = mean(VectorOut);
    return y;
    }
                      ")
n <- 13
getNPVs(rep(1,n), seq(1,n,1), n)

# 2.2) err.sqr.calc function ----

if(sdf.model == "linear") {
  f1 <- function(df.ss, max.month, par0) {
    return(getNPVs(df.ss$CF, 
                   exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0)))), 
                   max.month))
  }
}

if(sdf.model == "exp.aff") {
  f1 <- function(df.ss, max.month, par0) {
    df.ss$Alpha <- exp(1) - 1
    return(getNPVs(df.ss$CF, 
                   exp(cumsum(as.matrix(log(1 + df.ss[, names(par0)])) %*% par0)), 
                   max.month))
  }
}

if(error.function == "L2_Lasso") {
  err.sqr.calc <- function(par, max.month, lambda, df) {
    dfx <- split(df, df$Fund.ID)
    npvs <- sapply(dfx, f1, max.month=max.month, par0=c("RF" = 1, par))
    p <- ifelse(length(par) == 1, par, par[-1])
    return(
      sqrt(sum(npvs^2)) / length(npvs) + lambda / length(p) * sum(abs(p))
    )
  }
}

if(error.function == "L1_Ridge") {
  err.sqr.calc <- function(par, max.month, lambda, df) {
    dfx <- split(df, df$Fund.ID)
    npvs <- sapply(dfx, f1, max.month=max.month, par0=c("RF" = 1, par))
    p <- ifelse(length(par) == 1, par, par[-1])
    return(
      sum(abs(npvs)) / length(npvs) + lambda / length(p) * sum(p^2)
    )
  }
}

max.month <- 180
lambda <- 0
par <- c(MKT = 1, Alpha = 0)
type <- levels(df0$type)[1]
df.in <- df0[df0$type == type, ]

df.in$Fund.ID <- (as.character(df.in$Fund.ID))
system.time(
  print(err.sqr.calc(par, df=df.in, lambda=lambda, max.month=max.month))
)

# 2.3 Asymptotic gradient hessian -----
if(TRUE) {
  res <- optimx::optimx(par, err.sqr.calc, 
                        lambda = lambda,
                        max.month = max.month,
                        df = df.in,
                        #method = "nlminb"
                        method = "Nelder-Mead"
                        #method = c("Nelder-Mead", "L-BFGS-B", "nlminb", "nlm", "ucminf")
                        #control = list(all.methods=TRUE)
  )
  res # 2.156218 2.826169 135.062
  names.par <- names(par)
  par <- as.numeric(res[1, names.par])
  names(par) <- names.par
}

# Test Data
dfx <- split(df.in, df.in$Fund.ID)
df.ss <- dfx[[1]]
delta <- 0.0001


# GRADIENT
get.grad.per.fund <- function(df.ss, par) {
  # by central difference
  f2 <- function(x) {
    p.plus <- p.minus <- par
    p.plus <- par[x] + delta
    p.minus <- par[x] - delta
    return(
      (err.sqr.calc(p.plus, max.month, lambda, df.ss) - err.sqr.calc(p.minus, max.month, lambda, df.ss)) / (2 * delta)
    )
  }
  return(sapply(names(par), f2))
}
get.grad.per.fund(df.ss, par)

Avg.Grad <- colMeans(
  Reduce(rbind, lapply(dfx, get.grad.per.fund, par=par))
  )
Avg.Grad

weighter <- function(x, y, bw = kernel.bandwidth) {
  abs.dis <- abs( (x - y) / bw )
  ifelse(abs.dis <= 1, 1 - abs.dis, 0) # Bartlett kernel
  }
weighter(2010, 2015)

get.psi.matrix.dep <- function(dfx, par) {
  n <- length(dfx)
  A <- as.matrix(data.frame(do.call(rbind, lapply(dfx, get.grad.per.fund, par=par))))
  dim.GG <- length(A[1, ])
  M <- matrix(0, dim.GG, dim.GG)
  
  for (i in 1:n) {
    for (j in 1:n) {
      w <- weighter(dfx[[i]]$Vintage[1], dfx[[j]]$Vintage[1])
      if (w == 0) next
      GG <- A[i, ] %*% t(A[j, ])
      M <- M + w * GG
    }
  }
  return(M / n)
}

get.psi.matrix.indep <- function(dfx, par) {
  n <- length(dfx)
  A <- as.matrix(data.frame(do.call(rbind, lapply(dfx, get.grad.per.fund, par=par))))
  dim.GG <- length(A[1, ])
  M <- matrix(0, dim.GG, dim.GG)
  
  for (i in 1:n) {
    GG <- A[i, ] %*% t(A[i, ])
    M <- M + GG
  }
  return(M / n)
}

#system.time(psi <- get.psi.matrix2(dfx, par))

# HESSIAN
get.hess.per.fund <- function(df.ss, par) {
  n <- length(par)
  H <- matrix(0, n, n, dimnames = list(names(par), names(par)))
  for (i in names(par)) {
    for (j in names(par)) {
      
      if(i == j) {
        # https://en.wikipedia.org/wiki/Finite_difference
        p.p <- p.m <- par
        p.p[i] <- p.p[i] + delta
        p.m[j] <- p.m[j] - delta
        
        #out <- f1(df.ss, max.month, p.p) + f1(df.ss, max.month, p.m) - 2 * f1(df.ss, max.month, par)
        out <- err.sqr.calc(p.p, max.month, lambda, df.ss) + err.sqr.calc(p.m, max.month, lambda, df.ss) - 2 * err.sqr.calc(par, max.month, lambda, df.ss)
        out <- out / (delta * delta)
      }
      
      if(FALSE) {
        # https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
        p.pp <- p.mm <- p.p <- p.m <- par
        p.pp[i] <- p.pp[i] + 2 * delta
        p.mm[i] <- p.mm[i] - 2 * delta
        p.p[i] <- p.p[i] + delta
        p.m[j] <- p.m[j] - delta
        
        out <- - f1(df.ss, max.month, p.pp)
        out <- out + 16 * err.sqr.calc(p.p, max.month, lambda, df.ss)
        out <- out - 30 * err.sqr.calc(par, max.month, lambda, df.ss)
        out <- out + 16 * err.sqr.calc(p.m, max.month, lambda, df.ss)
        out <- out - err.sqr.calc(p.mm, max.month, lambda, df.ss)
        
        out <- out / (12 * delta * delta)
      }
      
      if(i != j) {
        # https://en.wikipedia.org/wiki/Finite_difference
        p.pp <- par + delta
        p.mm <- par - delta
        p.pm <- p.mp <- par
        p.pm[i] <- p.pm[i] + delta
        p.pm[j] <- p.pm[j] - delta
        p.mp[i] <- p.mp[i] - delta
        p.mp[j] <- p.mp[j] + delta
        
        out <- err.sqr.calc(p.pp, max.month, lambda, df.ss) - err.sqr.calc(p.pm, max.month, lambda, df.ss) - err.sqr.calc(p.mp, max.month, lambda, df.ss) + err.sqr.calc(p.mm, max.month, lambda, df.ss)
        out <- out  / (4 * delta * delta)
      }
      
      H[i, j] <-  out 
      
    }
  }
  return(H)
}
get.hess.per.fund(df.ss, par)


estim.CoVM <- function(dfx, par, indep=FALSE) {
  if(indep) {
    psi <- get.psi.matrix.indep(dfx, par)
  } else {
    psi <- get.psi.matrix.dep(dfx, par)
    
  }
  hessian <- Reduce('+', lapply(dfx, get.hess.per.fund, par=par)) / length(dfx)
  
  CoVM <- solve(hessian) %*% psi %*% t(solve(hessian))  / length(dfx)
  std.err <- sqrt(diag(CoVM))
  names(std.err) <- paste0("SE.", names(std.err))
  
  par.test <- par - c(1, rep(0, length(par)-1))
  wald <- 1 - pchisq(t(par.test) %*% CoVM %*% par.test, df = length(par.test))
  wald <- as.numeric(wald)
  
  out <- list(
    Psi = psi,
    Hessian = hessian,
    CoVM = CoVM,
    par = par,
    SE = std.err,
    Wald.p.value.MKT_1 = wald
  )
  return(out)
}
system.time(
  out <- estim.CoVM(dfx, par, indep=FALSE)
)
system.time(
  out <- estim.CoVM(dfx, par, indep=TRUE)
)
out



# 3.1) cross validation -------

vintage.blocks <- list()

if (do.cross.validation) {
  for(year in seq(1988, 2015, by=3)) {
    validate <- year:(year + 2)
    block <- (year - 3):(year - 1)
    block <- c(block, (year + 3):(year + 5))
    
    estimate <- 1980:2020
    estimate <- estimate[!(estimate %in% validate)]
    estimate <- estimate[!(estimate %in% block)]
    
    vintage.blocks[[paste(year)]] <- list(validate=validate, block=block, estimate=estimate)
  }
  
  for(year in names(vintage.blocks)) {
    for(key in c("estimate", "validate")) {
      vintage.blocks[[year]][[paste0("df.", key)]] <- df0[df0$Vintage %in% vintage.blocks[[year]][[key]], ]
      vintage.blocks[[year]][["cv.year"]] <- year
    }
  }
  
  viblo.list <- list()
  for(block in vintage.blocks) {
    viblo.list[[paste(block$validate[1])]] <- data.frame(
      training.before = paste0("start-", block$block[1]-1),
      h.block.before = paste(block$block[1:3], collapse = ","), 
      validation = paste(block$validate, collapse = ","), 
      h.block.after =paste(block$block[4:6], collapse = ","),
      training.after = paste0(block$block[6]+1, "-end")
      
    )
  }
  print(xtable::xtable(data.frame(do.call(rbind, viblo.list))), include.rownames = FALSE)
  
}

vintage.blocks[["ALL"]] <- list(df.estimate = df0, df.validate = df0, cv.year = "ALL")

names(vintage.blocks)

# 3.2) iterative run ----
library(optimx)
library(parallel)
library(doParallel)
library(foreach)

iter.run <- function(input.list) {
  
  # 1) INITIALIZE
  df.estimate <- input.list$df.estimate
  df.validate <- input.list$df.validate
  
  if (public.filename == "msci_market_factors") {
    factors <- c("HML", "SMB", "HDY", "QLT", "MOM", "ESG", "LOV")
    types <- levels(df0$type)
  }
  if (public.filename == "q_factors") {
    factors <- colnames(df.public)[2:6] # q_factors
    #factors <- "MKT"
    #factors <- "Alpha"
    factors <- c(factors, "Alpha")
    #factors <- c("ALL", factors)
    
    types <- c("PE", "VC", "PD", "RE", "NATRES", "INF") # asset classes
    #types <- levels(df0$type)
  }
  
  # 2) RUNNER
  l <- list()
  for (lambda in lambdas) {
    for(max.month in max.months) {
      print(paste("Max.M", max.month, "Lambda", lambda))
      
      if(do.parallel) {
        cl <- parallel::makeForkCluster(3)
        doParallel::registerDoParallel(cl)    
      }

      output <- foreach (type = types, .combine = "rbind", .verbose = FALSE) %dopar% {
        df.optim.in <- df.estimate[df.estimate$type == type, ]
        df.optim.in$Fund.ID <- (as.character(df.optim.in$Fund.ID))
        df.val <- df.validate[df.validate$type == type, ]
        df.val$Fund.ID <- (as.character(df.val$Fund.ID))

        factor.list <- list()
        for (factor in factors) {
          if(nrow(df.optim.in) == 0 | nrow(df.val) == 0) next
          
          par <- c("MKT" = 1)
          
          if(factor == "ALL") {
            for(fac in setdiff(factors, "ALL")) {
              par[fac] <- 0
            }
          } else if (factor == "Alpha") {
            par <- c("MKT" = 1, "Alpha"=0)
          } else {
            if ("ALL" %in% factors) break
            if(!(factor %in% names(par))) par[factor] <- 0
          }
          print(type)
          print(par)

          if (TRUE) {
            res <- optimx::optimx(par, err.sqr.calc, 
                                  lambda = lambda,
                                  max.month = max.month,
                                  df = df.optim.in,
                                  method = "Nelder-Mead" # "nlminb"
            )
            
            # if optimx does not terminate
            if(sum(res[1, names(par)] - par) < 0.001) {
              res <- optimx::optimx(par, err.sqr.calc, 
                                    lambda = lambda,
                                    max.month = max.month,
                                    df = df.optim.in,
                                    itnmax = 1000,
                                    method = "Nelder-Mead"
              )
            }
          } else {
            res <- data.frame(MKT=1, HML=1, SMB=1, HDY=1, QLT=1, MOM=1, ESG=1, LOV=1)
          }

          # Standard Error calc
          for(v in names(par)) par[v] <- res[1, v]
          dfx <- split(df.optim.in, df.optim.in$Fund.ID)
          
          # dependent - kernel bandwidht
          cov <- estim.CoVM(dfx, par=par, indep = FALSE)
          print(cov)
          for(v in names(cov$SE)) res[, v] <- cov$SE[v]
          res$Wald.p.value.MKT_1 <- cov$Wald.p.value.MKT_1
          
          # independet with other funds
          cov.indep <- estim.CoVM(dfx, par=par, indep = TRUE)
          for(v in names(cov.indep$SE)) res[, paste0(v, ".indep")] <- cov.indep$SE[v]
          
          print(res)
          res$Factor <- factor
          res$Type <- type
          res$max.month <- max.month
          res$lambda <- lambda
          res$kernel.bandwidth <- kernel.bandwidth
          res$CV.key <- input.list$cv.year
          res$weighting <- weighting
          res$error.fun <- error.function
          res$sdf.model <- sdf.model
          res$datetime <- Sys.time() + 2 * 60 * 60
          
          # validation error
          res$validation.error <- err.sqr.calc(par, 
                                               max.month = max.month, 
                                               lambda = 0, # set to zero
                                               df = df.val)
          
          factor.list[[factor]] <- res
        }
        return(data.frame(Reduce(rbind.all.columns, factor.list)))
      }
      
      if(do.cache & nrow(output) > 0) {
        cache.path <- paste0("data_out/cache_", public.filename, "/")
        if(!dir.exists(cache.path)) dir.create(cache.path)
        path.cache <- paste0(cache.path, 
                             format(Sys.time() + 2 * 60 * 60, "%Y-%m-%d_%H%M%S"), 
                             "_cached_res.csv")
        write.csv(output, path.cache)
      }
      
      print(output)
      l[[paste0(max.month, lambda)]] <- output
      
      if(do.parallel) parallel::stopCluster(cl)
    }
  }
  df.res <- data.frame(Reduce(rbind.all.columns, l))
  return(df.res)
}

if(do.cross.validation) {
  for(v in names(vintage.blocks)) df.res <- iter.run(vintage.blocks[[v]])
} else {
  system.time(df.res <- iter.run(vintage.blocks$ALL))
}


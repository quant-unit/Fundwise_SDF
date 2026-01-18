#### estimate model
# 0) Prologue -----

if(!exists("source.internally", envir = .GlobalEnv)) {
  source.internally <- TRUE
}

if (source.internally) {
  
  if(FALSE) {
    files <- c("prep_preqin.R", "prep_public.R")
    for(file in files) {
      source(file)
      if(sys.nframe() == 0L) rm(list = ls())
    }
  }
  
  if(sys.nframe() == 0L) rm(list = ls())
  source.internally <- TRUE
  
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  getwd()
  
}

# 1.1) PARAMETERS ----
if (source.internally) {
  
  export.data <- FALSE
  
  use.vintage.year.pfs <- TRUE
  use.simulation <- FALSE
  do.cross.validation <- FALSE
  do.cache <- TRUE
  do.parallel <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)
  
  private.source <- "pitchbook"
  private.source <- "preqin"
  cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
  cutoff <- "_cutoff_2021"
  
  #public.filename <- "public_returns" # outdated
  #public.filename <- "msci_market_factors"
  public.filename <- "q_factors"
  #public.filename <- "DebtFactorsEURUSD" # outdated
  #public.filename <- "iBoxxFactorsMIX"
  #public.filename <- "iBoxxFactorsUSD" 
  #public.filename <- "iBoxxFactorsEUR"
  
  
  # CHOICES
  weighting <- "EW"
  weighting <- "FW"
  #error.function <- "L1_Ridge"
  error.function <- "L2_Lasso"
  
  sdf.model <- "linear"
  #sdf.model <- "exp.aff"
  
  max.months <- c(1, 60, 120, 150, 180, 210, 240, 300, 360) # c(10, 20) * 12 # c(12.5, 15, 17.5) * 12
  max.months <- c(1, 60, 120, 150, 180)
  
  include.alpha.term <- FALSE
  lambdas <- 0
  kernel.bandwidth <- 12
  if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
  cache.folder.tag <- paste0(private.source, "_")
  # cache.folder.tag <- paste0(private.source, ifelse(do.cross.validation, "_cv_", "_"))
  cache.folder.tag <- paste0(cache.folder.tag, ifelse(include.alpha.term, "alpha_", ""))
  cache.folder.tag <- paste0(cache.folder.tag, weighting)
  cache.folder.tag
  
  # cache.folder.tag <- "20250808_222540_simulated_cashflows_EW"
  simulation.filename <- paste0("20250808_222540/", cache.folder.tag, ".csv")
  
  part.to.keep <- 1
  no.partitions <- 1 # 10
  
  data.out.folder <- "results/data_out_2026-emp-B"
  if(!dir.exists(data.out.folder)) dir.create(data.out.folder)
  
  factors.to.use <- ""
}

# 1.2) load data -----

data.prepared.folder <- "data_prepared"

# load public data
if (public.filename == "q_factors") {
  df.public <- read.csv(paste0("empirical/", data.prepared.folder, "/", public.filename, ".csv"))
  
} else {
  df.public <- read.csv2(paste0("empirical/", data.prepared.folder, "/", public.filename, ".csv"))
}
colnames(df.public) <- gsub("_World", "", colnames(df.public))
df.public$Date <- as.Date(df.public$Date)
bond.files <- c(
  "DebtFactorsUSD", "DebtFactorsEUR", "DebtFactorsEURUSD",
  "iBoxxFactorsEUR", "iBoxxFactorsUSD", "iBoxxFactorsMIX"
  )
if (public.filename %in% bond.files ) {
  df.public <- df.public[complete.cases(df.public[, 2:5]), ]
  df.public[is.na(df.public)] <- 0
  public.equity <- "msci_market_factors"
  df.pubeq <- read.csv2(paste0(data.prepared.folder, "/", public.equity, ".csv"))
  colnames(df.pubeq) <- gsub("_World", "", colnames(df.pubeq))
  df.pubeq$Date <- as.Date(df.pubeq$Date)
  
  df.public$RF <- NULL
  df.public <- merge(df.public, df.pubeq[, c("Date", "RF", "MKT")], by = "Date", all.x = TRUE)
  remove(df.pubeq)
  df.public <- df.public[df.public$Date > as.Date("2000-01-31"), ]
}

df.public$Alpha <- 1

if (export.data) {
  write.csv2(df.public, paste0("empirical/data_private_public/public_", public.filename, ".csv"))
}

# Load private data
if(!use.simulation) {
  
  if (private.source == "preqin") {
    year.tag <- ""
    df.private.cfs <- read.csv(paste0("empirical/", data.prepared.folder, "/preqin_cashflows_", weighting, year.tag, "_NAV.csv"))
    colnames(df.private.cfs)
  }
  
  if (private.source == "pitchbook") {
    year.tag <- "_2023"
    df.private.cfs <- read.csv(paste0("empirical/", data.prepared.folder, cutoff, "/pitchbook_cashflows_", weighting, year.tag, ".csv"))
    colnames(df.private.cfs)
    print(max(df.private.cfs$Date))
  }
  
} else {
  # df.private.cfs <- read.csv(paste0("simulation/", data.prepared.folder, "/simulated_cashflows_", weighting, ".csv"))
  df.private.cfs <- read.csv(paste0("simulation/data_prepared_sim/", simulation.filename))
  
}

# split too large data.frames into partitions
if (no.partitions > 1) {
  data.table::setDT(df.private.cfs)
  lev <- unique(df.private.cfs$Fund.ID)
  part <- data.table::data.table(Fund.ID = lev, partition = as.integer(dplyr::ntile(seq_along(lev), no.partitions)))
  df.private.cfs <- as.data.frame(part[df.private.cfs, on = "Fund.ID"][partition == part.to.keep])
  cache.folder.tag <- paste0(cache.folder.tag, "_part", part.to.keep)
  rm(lev) ; rm(part)
  df.private.cfs$partition <- NULL
}

df.private.cfs$Date <- as.Date(df.private.cfs$Date)
df.private.cfs$Fund.ID <- as.factor(paste(df.private.cfs$Fund.ID, df.private.cfs$type, sep = "_"))


to.monthly <- function(df.ss) {
  # fill zero cash flows (necessary for estimation)
  if (use.simulation) {
    max.date <- as.Date("2020-01-01")
  } else {
    max.date <- max(df.ss$Date) + 1
  }
  df.m <- data.frame(Date = seq(min(df.ss$Date) + 1, max.date, by = "month") - 1)
  df.ss <- merge(df.ss, df.m, by = "Date", all = TRUE)
  df.ss$CF[is.na(df.ss$CF)] <- 0
  for(col in c("type", "Vintage", "Fund.ID")) df.ss[, col] <- df.ss[1, col]
  return(df.ss)
}

if (FALSE) {
  list.private <- list()
  i <- 1
  for (fuid in levels(df.private.cfs$Fund.ID)) {
    print(i)
    i <- i + 1
    df.ss <- df.private.cfs[df.private.cfs$Fund.ID == fuid, ]
    list.private[[as.character(fuid)]] <- to.monthly(df.ss)
  }
}

#split.private <- split(df.private.cfs, df.private.cfs$Fund.ID)
#system.time(list.private <- lapply(split.private, to.monthly))
#system.time(df.private.cfs <- as.data.frame(data.table::rbindlist(list.private)))
system.time(df.private.cfs <- as.data.frame(data.table::rbindlist(lapply(split(df.private.cfs, df.private.cfs$Fund.ID), to.monthly))))


df.private.cfs$type <- as.factor(as.character(df.private.cfs$type))
df.private.cfs$Fund.ID <- as.factor(as.character(df.private.cfs$Fund.ID))
length(levels(df.private.cfs$type))



# merge private and public data
df0 <- merge(df.private.cfs, df.public, by="Date")
df0$Fund.ID <- as.factor(df0$Fund.ID)
rm(df.private.cfs)

# check if we have enough public data
min(df0$Vintage)
min(df0$Date)

if (public.filename %in% c("DebtFactorsUSD", "DebtFactorsEUR", "DebtFactorsEURUSD") ) {
  df0 <- df0[df0$Vintage > 2000, ]  
}
#df0 <- df0[df0$Vintage <= 2011, ]


# Export df0
if (export.data) {
  write.csv2(df0, paste0("data_private_public/df0_",
                         private.source, "_",
                         public.filename,"_",
                         weighting,
                         cutoff, ".csv"), 
             row.names = FALSE)
}


# Summary statistics for df0
aggregate(Vintage ~ type , df0, min)
number.of <- function(x) length(table(x))
aggregate(Vintage ~ type , df0, number.of)
aggregate(as.character(Fund.ID) ~ type , df0, number.of)

if (!use.simulation) {
  # print summary: funds per vintage
  df1 <- df0
  df1 <- df1[!duplicated(df1$Fund.ID), ]
  df1 <- as.data.frame.matrix(table(df1$Vintage, df1$type))
  # df1 <- df1[, colnames(df1) %in% c("BO", "DD", "INF", "MEZZ", "NATRES", "PD", "RE", "VC")]
  df1["Total", ] <- as.integer(colSums(df1))
  print(xtable::xtable(df1, caption = "Number of funds per vintage year.", 
                       label =paste0("tab:", private.source, "_data")), include.rownames = TRUE)
  rm(df1)
}

# 2.1) getNPVs function ----
source("helper/getNPVs.R")

# 2.2) err.sqr.calc function ----

if(sdf.model == "linear") {
  f1 <- function(df.ss, max.month, par0) {
    return(getNPVs(df.ss$CF, 
                   exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0)))), 
                   max.month))
  }
}

if(sdf.model == "linear.duration") {
  f1 <- function(df.ss, max.month, par0) {
    return(getNPVsDuration(df.ss$CF, 
                   exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0)))), 
                   max.month))
  }
}

if(sdf.model == "linear.single.date") {
  f1 <- function(df.ss, max.month, par0) {
    return(getNPVsSingleDate(df.ss$CF, 
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


# 2.3) Asymptotic gradient hessian -----
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
    if (use.simulation) types <- levels(df0$type)
  }
  if (public.filename %in% bond.files) {
    factors <- c("TERM", "CORP", "HY", "LIQ")
    types <- levels(df0$type)
    # types <- c("PD", "DD", "MEZZ")
  }
  
  if (factors.to.use != "") {
    factors <- factors.to.use
  }
  
  # 2) RUNNER
  l <- list()
  for (lambda in lambdas) {
    for(max.month in max.months) {
      print(paste("Max.M", max.month, "Lambda", lambda, " - ", Sys.time()))
      
      if(do.parallel) {
        cl <- parallel::makeForkCluster(3)
        doParallel::registerDoParallel(cl)    
      }

      output <- foreach (type = types, .combine = "rbind", .verbose = FALSE) %do% {
        df.optim.in <- df.estimate[df.estimate$type == type, ]
        df.optim.in$Fund.ID <- (as.character(df.optim.in$Fund.ID))
        df.val <- df.validate[df.validate$type == type, ]
        df.val$Fund.ID <- (as.character(df.val$Fund.ID))

        factor.list <- list()
        for (factor in factors) {
          if(nrow(df.optim.in) == 0 | nrow(df.val) == 0) next
          
          if(!include.alpha.term) {
            # two factor model
            par <- c("MKT" = 1)
          } else {
            # three factor model
            par <- c("MKT" = 1, "Alpha"=0)
          }
          
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
          # print(type)
          # print(par)

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
          # print(cov)
          for(v in names(cov$SE)) res[, v] <- cov$SE[v]
          res$Wald.p.value.MKT_1 <- cov$Wald.p.value.MKT_1
          
          # independet with other funds
          cov.indep <- estim.CoVM(dfx, par=par, indep = TRUE)
          for(v in names(cov.indep$SE)) res[, paste0(v, ".indep")] <- cov.indep$SE[v]
          
          # print(res)
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
        cache.path <- paste0(data.out.folder, "/cache_", paste0(public.filename, "_", cache.folder.tag), "/")
        if(!dir.exists(cache.path)) dir.create(cache.path)
        path.cache <- paste0(cache.path, 
                             format(Sys.time() + 2 * 60 * 60, "%Y-%m-%d_%H%M%S"), 
                             "_cached_res.csv")
        write.csv(output, path.cache)
      }
      
      # print(output)
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


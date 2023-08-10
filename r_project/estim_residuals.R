###### estimate residuals
# 0) Prologue -----
if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


# 1) load data -------
use.vintage.year.pfs <- TRUE
public.filename <- "q_factors"
private.source <- "preqin"
weighting <- "FW"
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
sub.folder <- paste(public.filename, private.source, weighting, sep = "#")

# load df0
file.name <- paste("df0", private.source, public.filename, weighting, sep ="_")
file.name <- paste(file.name, "csv", sep=".")
file.name
df0 <- read.csv2(paste0("data_private_public/", file.name))
df0$Date <- as.Date(df0$Date)

rm(file.name)

# load df.public
df.public <- read.csv2(paste0("data_private_public/public_", public.filename, ".csv"))
df.public$Date <- as.Date(df.public$Date)

# set parameters
max.month <- 180
lambda <- 0


# source getNPVs.R
source("getNPVs.R")

# load cached SDF factor estimates
sdf.model <- paste0("cache_q_factors_", weighting, "_SL")
filenames <- list.files(paste0("data_out/", sdf.model), pattern="*cached_res.csv", full.names=TRUE)
df.sdf <- data.frame(Reduce(rbind.all.columns, lapply(filenames, read.csv)))
df.sdf[is.na(df.sdf)] <- 0
types <- levels(as.factor(df.sdf$Type))
type <- types[1]
type <- "PE"

sdf.factors <- colnames(df.sdf)[grep(".indep", colnames(df.sdf))]
sdf.factors <- sub(".indep", "", sdf.factors)
sdf.factors <- sub("SE.", "", sdf.factors)
sdf.factors

avergae.sdf.par <- function(type) {
  print(type)
  df <- df.sdf[, c("Type", sdf.factors)]
  df <- aggregate(df[, sdf.factors], list(Type = df[, "Type"]), mean)
  df <- df[df$Type == type, sdf.factors]
  par <- as.numeric(as.character(df[1, ]))
  names(par) <- colnames(df[1, ])
  return(par)
}
par <- avergae.sdf.par(type)
par

pars <- list("AvgPar" = par)
pars

ensemble.of.par <- function(type) {
  df <- df.sdf[df.sdf$Type == type, ]
  df <- df[, sdf.factors]
  y <- list()
  sdf.names <- colnames(df)
  for (i in 1:nrow(df)) {
    x <- as.numeric(df[i, ])
    names(x) <- sdf.names
    y[[paste0("Ensemble", i)]] <- x
  }
  return(y)
}
pars <- ensemble.of.par(type)



# 2) alpha return component -------------

# IMPORTANT: this is an alternative way to determine the alpha term
# this method here is conditional a factor model that has no intercept (= alpha term)
# compare to: Excess-IRR method described by Phalippou and Gottschalg (2009)

# make df.idi
df.idi <- data.frame(Date = unique(df0$Date),
                     idi.return = 0)
df.idi$date.chr <- as.character(df.idi$Date)

# linear SDF
f1.alpha <- function(df.ss, max.month, par0, par.alpha, df.idi) {
  
  # where to place small return?
  df.ss <- merge(df.ss, df.idi, by = "Date")
  df.ss$idi.return = par.alpha
  
  return(getNPVs(df.ss$CF, 
                 exp(cumsum(log(
                   1 + (as.matrix(df.ss[, names(par0)]) %*% par0)
                   + df.ss$idi.return
                 ))), 
                 max.month))
}
f1.alpha(df0, 120, par, 0, df.idi)

# L2 lasso
err.sqr.calc.alpha <- function(par, max.month, lambda, df, df.idi, par.factor) {
  dfx <- split(df, df$Fund.ID)
  npvs <- sapply(dfx, f1.alpha, max.month=max.month, par0=c("RF" = 1, par.factor), 
                 par.alpha=par, df.idi)
  p <- ifelse(length(par) == 1, par, par[-1])
  return(
    sqrt(sum(npvs^2)) / length(npvs) + lambda / length(p) * sum(abs(p))
  )
}
err.sqr.calc.alpha(0, 120, 0, df0, df.idi, par)

#par <- c(MKT = 1.33, HML = -0.15, SMB = 0.2, HDY = 0.3, QLT = 0.21) # from SDF paper for BO (2020 Preqin)
#par <- c(MKT = 1.4, HML = -0.1, SMB = -0.09, HDY = -0.15, QLT = -0.1) # from SDF paper for BO (2023 Pitchbook)

# optimize alpha
res.alpha <- optimize(err.sqr.calc.alpha, 
                      interval = c(-1, 1),
                      lambda = lambda,
                      max.month = max.month,
                      df = df0[df0$type == type, ],
                      par.factor = par,
                      df.idi = df.idi)
(1 + res.alpha$minimum)^12 - 1

# 3) idiosyncratic return component ------

map.types = list(
  INF = "INF_FUNDS",
  NATRES = "NATRES_FUNDS",
  PD = "PE_PD",
  PE = "ALL_ALL",
  RE ="RE_FUNDS",
  VC = "PE_VC"
)

# TODO: align monthly and quarterly returns !!!
make.df.idi <- function(use.nav.returns=FALSE, par0=NA, df.public=NA, col.return="PE_BO") {
  
  # make df.idi
  df.idi <- data.frame(Date = unique(df0$Date),
                       idi.return = 0)
  df.idi$date.chr <- as.character(df.idi$Date)
  
  if ( use.nav.returns ) {
    
    if (FALSE) {
      # load nav returns (QUARTERLY NAV RETURNS)
      file.nav.returns <- "nav_returns_pitchbook_acs"
      df.nav.ret <- read.csv2(paste0("nav_returns/", file.nav.returns, ".csv"))

      quarter2date <- function(x) {
        year <- substring(x, 1, 4)
        quarter <- as.numeric(substring(x, 6, 6))
        date <- paste(year, ifelse(quarter==1, "03-31",
                                   ifelse(quarter==2, "06-30",
                                          ifelse(quarter==3, "09-30", "12-31"))), sep ="-")
        date <- as.Date(date)
        return(date)
      }
      df.nav.ret$Date <- quarter2date(df.nav.ret$AsOf)
      df.nav.ret <- df.nav.ret[df.nav.ret$asset_class_segment == "PE_BO", ]
      df.idi <- merge(df.idi, df.nav.ret[, c("Date", "StandardReturn_filledQ")], by="Date", all.x = TRUE)
      df.idi$StandardReturn_filledQ[is.na(df.idi$StandardReturn_filledQ)] <- 0
      df.idi$idi.return <- df.idi$StandardReturn_filledQ
      df.idi$StandardReturn_filledQ <- NULL
    } else if (FALSE) {
      file.nav.returns <- "Pitchbook NAV Returns Truncated"
      # newly calculated Pitchbook NAV returns from AM data_prep (also QUARTERLY !!!)
      df.nav.ret <- read.csv(paste0("nav_returns/", file.nav.returns, ".csv"))
      df.nav.ret$Date <- as.Date(df.nav.ret$date)
      df.idi <- merge(df.idi, df.nav.ret[, c("Date", col.return)], by="Date", all.x = TRUE)
      df.idi$idi.return <- as.numeric(gsub(",", ".", df.idi[, col.return]))
      df.idi$idi.return[is.na(df.idi$idi.return)] <- 0
      df.idi[, col.return] <- NULL
    } else {
      # CA NAV 100 return time series
      file.nav.returns <- "ca_index_100"
      df.nav.ret <- read.csv(paste0("nav_returns/", file.nav.returns, ".csv"))
      df.nav.ret$Date <- as.Date(df.nav.ret$Date)
      df.idi <- merge(df.idi, df.nav.ret[, c("Date", col.return)], by="Date", all.x = TRUE)
      df.idi$idi.return <- as.numeric(gsub(",", ".", df.idi[, col.return]))
      df.idi$idi.return[is.na(df.idi$idi.return)] <- 0
    }
    
    
    # calculate residuals by subtracting the factor return
    df.helper <- merge(df.idi, df.public, by="Date")
    par0 <- c("RF" = 1, par0)
    factor.return <- as.matrix(df.helper[, names(par0)]) %*% par0
    df.idi$total.return <- df.idi$idi.return
    df.idi$factor.return <- factor.return
    df.idi$idi.return <- df.idi$idi.return - factor.return
  }
  
  return(df.idi)
}
df.idi <- make.df.idi(TRUE, par, df.public, map.types[[type]])

plot(df.idi$Date, df.idi$idi.return, type="l", xlab="Date", ylab="residuals", main="NAV Return Residuals (Monthly/Quarterly)")
lines(df.idi$Date, df.idi$total.return, col="red")
legend("topright", bty="n", legend = c("Total return", "Residual"), col = c("red", "black"), lty=1)
abline(h=0, col="blue")

# Compare NAV Returns to Public Returns
if ( FALSE ) {
  df.msci <- df.public
  df.msci <- df.msci[df.msci$Date > as.Date("1996-12-31"), ]
  df.msci$Return <- df.msci$RF + df.msci$MKT
  df.msci$Total.Return <- cumprod(1 + df.msci$Return)
  par0 <- c("RF" = 1, par)
  df.msci$Factor.Return <- as.matrix(df.msci[, names(par0)]) %*% par0
  df.msci$Total.Factor.Return <- cumprod(1 + df.msci$Factor.Return)
  rownames(df.msci) <- NULL
  just.quarter.ends <- seq(1, 92) * 3
  df.msci <- df.msci[just.quarter.ends, ]
  df.msci$Total.Return.lag <- c(1, df.msci$Total.Return[-nrow(df.msci)])
  df.msci$Total.Return.Q <- df.msci$Total.Return / df.msci$Total.Return.lag - 1
  df.msci$Total.Factor.Return.lag <- c(1, df.msci$Total.Factor.Return[-nrow(df.msci)])
  df.msci$TotalFactor.Return.Q <- df.msci$Total.Factor.Return / df.msci$Total.Factor.Return.lag - 1
  
  returns.bo.quarterly <- df.idi$total.return[df.idi$total.return != 0]
  plot(acf(returns.bo.quarterly, 5), main="Quarterly Buyout NAV Return Autocorrelation (1997-2019)")
  legend("topright", bty="n", cex=1.2, legend=c(paste("Return Mean:", round(mean(returns.bo.quarterly), 3)),
                                                paste("Return StDev:", round(sd(returns.bo.quarterly), 3))))
  plot(acf(df.msci$Total.Return.Q, 5), main="Quarterly MSCI World Return Autocorrelation (1997-2019)")
  legend("topright", bty="n", cex=1.2, legend=c(paste("Return Mean:", round(mean(df.msci$Total.Return.Q), 3)),
                                                paste("Return StDev:", round(sd(df.msci$Total.Return.Q), 3))))
  plot(acf(df.msci$TotalFactor.Return.Q, 5), main="Quarterly Factor Return Autocorrelation (1997-2019)")
  legend("topright", bty="n", cex=1.2, legend=c(paste("Return Mean:", round(mean(df.msci$TotalFactor.Return.Q), 3)),
                                                paste("Return StDev:", round(sd(df.msci$TotalFactor.Return.Q), 3))))
  summary(returns.bo.quarterly)
  sd(returns.bo.quarterly)
  summary(df.msci$Total.Return.Q)
  sd(df.msci$Total.Return.Q)
  summary(df.msci$TotalFactor.Return.Q)
}


# linear SDF
f1.idi <- function(df.ss, max.month, par0, par.idi.date, par.idi, df.idi) {
  
  # where to place small return?
  df.ss <- merge(df.ss, df.idi, by = "Date")
  df.ss$idi.return[df.ss$Date == par.idi.date] = df.ss$idi.return[df.ss$Date == par.idi.date] + par.idi
  
  return(getNPVs(df.ss$CF, 
                 exp(cumsum(log(
                   1 + (as.matrix(df.ss[, names(par0)]) %*% par0)
                   + df.ss$idi.return
                 ))), 
                 max.month))
}

# L2 lasso
err.sqr.calc.idi <- function(par, max.month, lambda, df, par.idi.date, df.idi, par.factor) {
  dfx <- split(df, df$Fund.ID)
  npvs <- sapply(dfx, f1.idi, max.month=max.month, par0=c("RF" = 1, par.factor), 
                 par.idi.date=par.idi.date, par.idi=par, df.idi)
  p <- ifelse(length(par) == 1, par, par[-1])
  return(
    sqrt(sum(npvs^2)) / length(npvs) + lambda / length(p) * sum(abs(p))
  )
}

RETURN.BOUNDARY <- 1

# one boosting step
one.cwb.iteration <- function(df.idi, lambda, max.month, df.in, par) {
  
  # initialize
  best.res <- list(objective = Inf)
  
  # find best date
  for (date in df.idi$date.chr) {
    # print(date)
    res <- optimize(err.sqr.calc.idi, 
                    interval = c(-RETURN.BOUNDARY, RETURN.BOUNDARY),
                    lambda = lambda,
                    max.month = max.month,
                    df = df.in,
                    par.factor = par,
                    df.idi = df.idi,
                    par.idi.date = date
    )
    res
    
    if (res$objective < best.res$objective) {
      best.res <- res
      best.res$date <- date
    }
  }
  
  return(data.frame(best.res))
}

# system.time(res <- one.cwb.iteration(df.idi, lambda, max.month, df.in, par))
# res

one.cwb.iteration2 <- function(df.idi, lambda, max.month, df.in, par) {
  
  # find best date
  foo <- function(date) {
    res <- optimize(err.sqr.calc.idi, 
                    interval = c(-RETURN.BOUNDARY, RETURN.BOUNDARY),
                    lambda = lambda,
                    max.month = max.month,
                    df = df.in,
                    par.factor = par,
                    df.idi = df.idi,
                    par.idi.date = date
    )
    res$date <- date
    return(data.frame(res))
  }
  
  method <- "parLapply" # Windows
  method <- "mclapply" # Mac
  no.cores <- parallel::detectCores() - 1
  
  if (method == "mclapply") {
    # Unfortunately, 'mc.cores' > 1 not supported in Windows
    mc.cores <- ifelse(.Platform$OS.type == "unix", no.cores, 1)
    list.res <- parallel::mclapply(df.idi$date.chr, foo, mc.cores = mc.cores)
  } else if (method == "parLapply") {
    # does not work with Rccp functions
    cl <- parallel::makeCluster(no.cores)
    export2cluster <- c("err.sqr.calc.idi", "RETURN.BOUNDARY", 
                        "df.in", "f1.idi", "df.idi", "getNPVs",
                        "lambda", "max.month", "par")
    parallel::clusterExport(cl, export2cluster, envir=.GlobalEnv)
    #parallel::clusterExport(cl, export2cluster, envir=environment())
    
    list.res <- parallel::parLapply(cl, df.idi$date.chr, foo)
    parallel::stopCluster(cl)
  } else {
    list.res <- lapply(df.idi$date.chr, foo)
  }
  
  df.res <- data.frame(do.call(rbind, list.res))
  best.res <- df.res[df.res$objective == min(df.res$objective), ]
  
  return(best.res)
}

#system.time(res <- one.cwb.iteration2(df.idi, lambda, max.month, df.in, par))
#res

df.in <- df0[df0$type == type, ]

# multiple boosting steps
multi.cwb.iterations <- function(n.boost, 
                                 lambda, max.month, df0, par, 
                                 type,
                                 alpha=0, 
                                 old.out = NA, 
                                 load.cached.idi = FALSE,
                                 use.nav.returns = FALSE,
                                 df.public = NA, 
                                 ensemble="") {
  
  df.in <- df0[df0$type == type, ]
  
  # load cached csv
  if (load.cached.idi) {
    df.idi <- read.csv2("data_idi/df.idi.csv")
    df.res.before <- read.csv2("data_idi/df.res.csv")
    
    df.idi$X <- NULL
    df.idi$Date <- as.Date(df.idi$Date)
    
    df.res.before$X <- NULL
  }
  
  if ( is.list(old.out)) {
    df.idi = old.out$df.idi
    df.res.before = old.out$df.res
  } else if (!load.cached.idi) {
    df.res.before <- NA
    df.idi <- NA
  }
  
  # make df.idi
  if (!is.data.frame(df.idi)) {
    df.idi <- make.df.idi(use.nav.returns, par, df.public)
    
    if (!use.nav.returns) df.idi$idi.return <- alpha
  }

  # calc factor return
  df.fac <- df.in[!duplicated(df.in$Date), c("Date", names(par))]
  df.idi <- merge(df.idi, df.fac, by = "Date")
  df.idi$factor.return <- (as.matrix(df.idi[, names(par)]) %*% par)
  df.idi <- df.idi[, !(colnames(df.idi) %in% names(par))]
  
  # set hyper-parameters
  damper <- 0.33
  
  list.res <- list()
  # loop over boosting steps
  for (i in 1:n.boost) {
    print(paste(ensemble, i))
    proc.time <- system.time(res <- one.cwb.iteration2(df.idi, lambda, max.month, df.in, par))
    res$alpha <- alpha
    res$damper <- damper
    res$elapsed <- proc.time["elapsed"]
    
    # update idiosyncratic return time-series
    df.idi$idi.return[df.idi$date.chr == res$date] <- df.idi$idi.return[df.idi$date.chr == res$date] + res$minimum * damper
    
    # NICE TO KNOW
    print(paste("filled returns:", nrow(df.idi[df.idi$idi.return != alpha, ])))
    
    # store optimization result
    list.res[[as.character(i)]] <- res
    print(res)
    
  }
  df.res <- data.frame(do.call(rbind, list.res))
  
  if (is.data.frame(df.res.before)) {
    df.res <- rbind(df.res.before, df.res)
    rownames(df.res) <- NULL
  }
  
  df.par <- data.frame(t(par))
  
  # add type col
  df.par$type <- type
  df.res$type <- type
  df.idi$type <- type
  
  # add ensemble col
  df.idi$ensemble <- ensemble
  df.res$ensemble <- ensemble
  df.par$ensemble <- ensemble
  
  # prepare output
  out <- list(
    df.res = df.res,
    df.idi = df.idi,
    df.par = df.par
  )
  return(out)
}

write.out <- function(out, type, tag="") {
  # write
  write.csv2(out$df.idi, paste0("data_idi/",sub.folder, "/df_idi_", type, tag, ".csv"))
  write.csv2(out$df.res, paste0("data_idi/",sub.folder, "/df_res_", type, tag, ".csv"))
  write.csv2(out$df.par, paste0("data_idi/",sub.folder, "/df_par_", type, tag, ".csv"))
}

read.out <- function(type, tag="") {
  out <- list()
  for (name in c("idi", "par", "res")) {
    file <- paste0("data_idi/",sub.folder, "/df_",name, "_", type, tag,".csv")
    if(! file.exists(file)) return(NA)
    
    df <- read.csv2(file)
    if ("Date" %in% colnames(df)) df$Date <- as.Date(df$Date)
    df$X <- NULL
    out[[paste0("df.", name)]] <- df
  }
  return(out)
}
read.out(type)

res.alpha$minimum

no.iterations <- 200
print(paste("approximated run time:", length(pars) * no.iterations * 0.5 / 60, "hours for type", type))
for (ensemble in names(pars)) {
  out <- multi.cwb.iterations(no.iterations, lambda, max.month, df0, pars[[ensemble]], type, 
                              use.nav.returns=FALSE, df.public=df.public, ensemble=ensemble)
  # out <- multi.cwb.iterations(2, lambda, max.month, df.in, par, res.alpha$minimum, load.cached.idi=TRUE)
  # out <- multi.cwb.iterations(50, lambda, max.month, df.in, par, old.out = out)
  
  # combine old.out and new.out
  old.out <- read.out(type)
  if(is.list(old.out)) {
    for(name in names(old.out)) {
      out[[name]] <- rbind(old.out[[name]], out[[name]])
    }
  }
  write.out(out, type)
}


# 4. Analyze output -----
count.dates <- function(df.res) {
  y <- c()
  for (i in 1:nrow(df.res)) {
    df <- df.res[1:i, ]
    y <- c(y, length(unique(df$date)))
  }
  return(y)
}


plot(out$df.idi$Date, out$df.idi$idi.return, type="l", xlab = "Date", ylab = "Return", main=type)
plot(acf(out$df.idi$idi.return))
plot(out$df.res$objective, type = "l", xlab="iterations", ylab="objective", main=type)
legend("topright", bty="n", legend = paste("# different dates:", length(unique(out$df.res$date))))
plot(count.dates(out$df.res), type="s", col="blue", xlab="iterations", ylab="# dates", main=type)



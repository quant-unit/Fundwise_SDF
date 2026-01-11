#####################
### Ang et al. (2018) Estimating Private Equity Returns from Limited Partner Cash Flows
#####################
# 0) load data ------
if(sys.nframe() == 0L) rm(list = ls())

source("load_data.R")
source("getNPVs.R")

VERBOSE <- FALSE

upper.incomplete.gamma <- function(x, a, target.value=0) {
  return(pracma::incgam(x, a) - target.value)
}

inverse.upper.incomplete.gamma <- function(x, a, lower = 0.001, upper = 1000) {
  result <- uniroot(
    upper.incomplete.gamma, lower = lower, upper = upper, a = a, target.value = x
    )
  return(result$root)
}

update.sigma.delta.h <- function(v, N, var.min, var.max) {
  a <- 0.5 * v + 0.5 * var.min
  b <- 0.5 * N + 0.5 * var.max

  u <- try(
    runif(1, 
          min = upper.incomplete.gamma(b / var.min, a),
          max = upper.incomplete.gamma(b / var.max, a)
    ))
  
  while ( is.na(u) ) {
    v <- v / 2
    a <- 0.5 * v + 0.5 * var.min
    b <- 0.5 * N + 0.5 * var.max
    
    u <- try(
      runif(1, 
            min = upper.incomplete.gamma(b / var.min, a),
            max = upper.incomplete.gamma(b / var.max, a)
      ))

    print(paste("u", u, "a", a, "b", b))
  }

  sigma <- sqrt( b / inverse.upper.incomplete.gamma(u, a) )
  return(sigma)
}

# make df.idi
df.idi <- data.frame(Date = unique(df0$Date),
                     idi.return = 0)

# define Present Value Ratio according to Ang et al. 2018
if(sdf.model == "linear") {
  
  calc.linear.return <- function(df.ss, par0, only.errors = FALSE) {
    ret <- df.ss$idi.return # residual error term
    if(!only.errors) {
      ret <- ret + (as.matrix(df.ss[, names(par0)]) %*% par0) # beta * public.factors
    }
    return(ret)
  }
  
  PVR <- function(df.ss, par0, df.idi) {
    
    # add risk-free return (if not done before)
    # par0=c("RF" = 1, par0)
    
    # merge error terms
    df.ss$idi.return <- NULL
    df.ss <- merge(df.ss, df.idi, by = "Date")

    cf <- df.ss$CF
    ret <- exp(cumsum(log(1 + calc.linear.return(df.ss, par0, only.errors=TRUE)))) # residual term
    max_month <- 1
    
    con <- ifelse(cf < 0, cf, 0)
    dis <- ifelse(cf > 0, cf, 0)
    
    PV_con <- getNPVs(con, ret, max_month)
    PV_dis <- getNPVs(dis, ret, max_month)
    
    present.value.ratio <- log(PV_dis / - PV_con)
    
    return(present.value.ratio)
  }
  
}

# init
par <- c(MKT = 0.5, HML = -0.2)
sigma.par.init <- c(MKT = 0.1, HML = 0.05)   # beta factor (a.k.a. par) stdv
sigma.delta.init <- 0.01 # residual stdv
sigma.h.init <- 0.5     # PVR stdv

var.min <- 0.01
var.max <- 0.1

# update
sigma.par <- sigma.par.init
sigma.delta <- sigma.delta.init
sigma.h <- sigma.h.init

# calc likelihood function
dfx <- df0[df0$type == "PE", ]
dfx$Fund.ID <- as.factor(as.character(dfx$Fund.ID))
dfx <- split(dfx, dfx$Fund.ID)

calc.log.likelihood <- function(dfx, par, df.idi, sigma.h) {
  pvrs <- sapply(dfx, PVR, df.idi=df.idi, par0=c("RF" = 1, par))
  log.likelihood <- (sum(log(dnorm(pvrs, mean = - 0.5 * sigma.h^2, sd = sigma.h))))
  return(log.likelihood)
}
calc.log.likelihood(dfx, par, df.idi, sigma.h)

calc.log.prior.delta <- function(idi.return, sigma.delta) {
    log.prior <- sum(log(dnorm(idi.return, mean=0, sd = sigma.delta)))
    return(log.prior)
}
calc.log.prior.delta(df.idi$idi.return, sigma.delta)

calc.log.proposal.delta <- function(one.idi.return, sigma.delta) {
    return(calc.log.prior.delta(one.idi.return, sigma.delta))
}
calc.log.proposal.delta(0, 2)

MH_propose <- function(one.idi.return, sigma.delta) {
  proposed.delta <- -10
  while( proposed.delta < -0.5) { # avoid stark negative returns
    proposed.delta <- rnorm(1, mean = one.idi.return, sd = sigma.delta)
  }
  return(proposed.delta)
}

Acceptance.Fun <- function(dfx, par, df.idi, sigma.h, sigma.delta, delta.proposal, j) {
  df.idi.new <- df.idi
  df.idi.new$idi.return[j] <- delta.proposal
  delta.old <- df.idi$idi.return[j]
  
  # assuming symmetric transition kernel q(a,b) = q(b,a)
  ratio <- exp(0
    # new log.posterior
    + calc.log.likelihood(dfx, par, df.idi.new, sigma.h) # new (proposal)
    + calc.log.prior.delta(df.idi.new$idi.return, sigma.delta) # new (proposal)
    # + calc.log.proposal.delta(delta.old, sigma.delta) # old (previous)
    # old log.posterior
    - calc.log.likelihood(dfx, par, df.idi, sigma.h) # old (previous)
    - calc.log.prior.delta(df.idi$idi.return, sigma.delta) # old (previous)
    # - calc.log.proposal.delta(delta.proposal, sigma.delta) # new (proposal)
  )
  
  if (VERBOSE) {
    print( df.idi.new$idi.return[j] )
    print( calc.log.likelihood(dfx, par, df.idi.new, sigma.h) )
    print( calc.log.prior.delta(df.idi.new$idi.return, sigma.delta) )
    print( calc.log.likelihood(dfx, par, df.idi, sigma.h) )
    print( calc.log.prior.delta(df.idi$idi.return, sigma.delta) )
  }
  
  return(min(1, ratio))
}

type <- "PE"
par0 <- par
sigma.par0 <- sigma.par
iterations <- 100

mcmc.one.iteration <- function(
    type="PE", par0=par, sigma.par0=sigma.par,
    sigma.h0=sigma.h, sigma.delta0=sigma.delta,
    iterations=100) {
  
  # 0. Init
  beta.prior <- par0
  par0 <- c("RF" = 1, par0) # add risk-free rate
  
  # 0.1 make df.idi0
  df.idi0 <- data.frame(Date = unique(df0$Date), idi.return = 0)

  # 0.2 make df.in and dfx and df.fac
  df.in <- df0[df0$type == type, ]
  df.in$Fund.ID <- as.factor(as.character(df.in$Fund.ID))
  df.fac <- df.in[!duplicated(df.in$Date), c("Date", names(par0))]
  
  # 0.3 inital candidates
  df.idi0.candidat <- df.idi0
  sigma.h0.candidat <- sigma.h0
  sigma.delta0.candidat <- sigma.delta0
  par0.candidate <- par0
  sigma.par0.candidate <- sigma.par0
  
  list.of.results <- list()
  for (i in 1:iterations) {
    # 0.4 init in loop
    dfx <- split(df.in, df.in$Fund.ID)
    
    # 1. Update df.idi0$idi.return: Metropolis Hastings for delta
    N.rows <- nrow(df.idi0.candidat)
    # N.rows <- 3
    for (j in 1:N.rows) { # iterate over 1:nrow(df.idi0.candidat)
      one.idi.return <- df.idi0.candidat$idi.return[j]
      delta.proposal <- MH_propose(one.idi.return, sigma.delta0)
      
      AC <- Acceptance.Fun(dfx, par, df.idi0.candidat, sigma.h0.candidat, sigma.delta0, delta.proposal, j)
      
      if (VERBOSE) print(paste(j, one.idi.return, sigma.h0.candidat, sigma.delta0, delta.proposal, AC))
      
      if(is.na(AC)) {
        print("delta.proposal")
        print(delta.proposal)
        stop("Isch over") # next
      } 
      
      if (VERBOSE & (AC < 0.999)) print(paste(j, round(AC, 2)))
      
      if (runif(1) < AC) {
        df.idi0.candidat$idi.return[j] <- delta.proposal
      } # else nothing
    }
    
    # 2. Update par: Regression for beta
    # 2.1 calc factor returns
    df.idi0.candidat <- merge(df.idi0.candidat, df.fac, by="Date", all.x = TRUE)
    # 2.2 calc total returns
    df.idi0.candidat$return <- calc.linear.return(df.idi0.candidat, par0, only.errors=TRUE)
    df.idi0.candidat$return <- df.idi0.candidat$return - df.idi0.candidat$RF # excess return
    # 2.3 run regression
    indep.variables <- setdiff(names(par0), "RF")
    form <- as.formula(paste("return", "~", paste(indep.variables, collapse="+"), "+ 0"))
    reg <- lm(form, data=df.idi0.candidat)
    # 2.4 update beta
    sd.beta <- summary(reg)$coefficients[, "Std. Error"]
    cov_matrix <- vcov(reg) # Extract covariance matrix of regression coefficients
    cholesky_decomp <- data.frame(chol(cov_matrix)) # Perform Cholesky decomposition

    for(par.name in names(beta.prior)) { # add noise sigma_beta
      var.reg <- (sd.beta[par.name])^2
      var.prior <- (sigma.par0[par.name])^2
      w.prior <- var.reg / ( var.reg +  var.prior)
      beta <- (1 - w.prior) * coef(reg)[par.name] + w.prior * beta.prior[par.name]
      beta.sd <- sqrt(1 - w.prior) * cholesky_decomp[par.name, par.name]
      par0.candidate[par.name] <- rnorm(1, mean=beta, sd=beta.sd)
    }
    
    # 2.5 clean df.idi0.candidat
    df.idi0.candidat <- df.idi0.candidat[, !(colnames(df.idi0.candidat) %in% names(par0.candidate))]
    df.idi0.candidat$return <- NULL

    # 3. Update sigma.delta
    v <- t(reg$residuals) %*% reg$residuals
    sigma.delta0.candidat <- update.sigma.delta.h(v, N=N.rows, var.min=var.min, var.max=var.max)
    
    # 4. Update sigma.h0.candidat
    dfx <- split(df.in, df.in$Fund.ID)
    pvrs <- sapply(dfx, PVR, df.idi=df.idi0.candidat, par0=par0.candidate)
    v <- t(pvrs) %*% pvrs
    if ( sum(is.na(pvrs)) ) {
      print(paste("SKIP RUN", i))
      # RESET
      par0.candidate <- par0
      sigma.par0.candidate <- sigma.par0
      df.idi0.candidat <- df.idi0
      sigma.delta0.candidat <- sigma.delta0
      # sigma.h0.candidat <- sigma.h0 # not yet updated
      # SKIP
      next
    }
    sigma.h0.candidat <- update.sigma.delta.h(v, N=length(pvrs), var.min=var.min, var.max=var.max)
    sigma.h0.candidat
    
    # 5. use candidates for update
    par0 <- par0.candidate
    sigma.par0 <- sigma.par0.candidate
    df.idi0 <- df.idi0.candidat
    sigma.h0 <- sigma.h0.candidat
    sigma.delta0 <- sigma.delta0.candidat
    
    # X. RESULT
    result <- data.frame(t(par0), 
                         "sigma.delta0" = sigma.delta0,
                         "sigma.h0" = sigma.h0, 
                         "iter" = i)
    list.of.results[[paste(i)]] <- result
    print(result)
  }
  
  df.res <- data.frame(do.call(rbind, list.of.results))
  
  final.results <- list(
    "df.idi0" = df.idi0,
    "df.res" = df.res
  )
  
  return(final.results)
}

system.time(final.results <- mcmc.one.iteration(type="PE", iterations = 1000))
final.results$df.res
summary(final.results$df.res)

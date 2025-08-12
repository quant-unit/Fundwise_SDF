#####################
### Ang et al. (2018) Estimating Private Equity Returns from Limited Partner Cash Flows
### translated from original MatLab files by GPT-o1
### https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fjofi.12688&file=jofi12688-sup-0001-SupMat.zip
####################
# Fun: Mnpv_f ---------

Mnpv_f <- function(alpha, bta, Factors, RF, f, Out, In) {
  T <- nrow(RF)
  disc <- rep(1, T + 1)
  
  for (c in 2:(T + 1)) {
    disc[c] <- disc[c - 1] / (1 + RF[c - 1] + sum(Factors[c - 1, ] * bta) + alpha + f[c - 1])
  }
  
  pv <- numeric(ncol(Out))
  for (i in 1:ncol(Out)) {
    pv[i] <- sum(Out[, i] * disc) - sum(In[, i] * disc)
  }
  
  mnpv <- mean(pv)^2
  return(mnpv)
}

# Fun: MPME_full ------

MPME_full <- function(alpha, bta, Factors, RF, f, Out, In, weight) {
  T <- nrow(RF)
  disc <- rep(1, T + 1)
  
  for (c in 2:(T + 1)) {
    disc[c] <- disc[c - 1] / (1 + RF[c - 1] + sum(Factors[c - 1, ] * bta) + alpha + f[c - 1])
    # Alternative without 'f' term:
    # disc[c] <- disc[c - 1] / (1 + RF[c - 1] + sum(Factors[c - 1, ] * bta) + alpha)
  }
  
  pme <- numeric(ncol(Out))
  for (i in 1:ncol(Out)) {
    pme[i] <- sum(Out[, i] * disc) / sum(In[, i] * disc)
  }
  
  mPME <- (sum(pme * weight) / sum(weight) - 1)^2
  return(mPME)
}


# Fun: MPME_simple ------

MPME_simple <- function(F, Out, In, weight) {
  T <- length(F)
  disc <- rep(1, T + 1)
  
  for (c in 2:(T + 1)) {
    disc[c] <- disc[c - 1] / (1 + F[c - 1])
  }
  
  pme <- numeric(ncol(Out))
  for (i in 1:ncol(Out)) {
    pme[i] <- sum(Out[, i] * disc) / sum(In[, i] * disc)
  }
  
  mPME <- (sum(pme * weight) / sum(weight) - 1)^2
  return(mPME)
}

# Fun: PME ------

PME <- function(Factors, Out, In, RF, E_f_ret, E_g_ret, DivAll, InvAll, theta_mat, Num_Factors) {
  T <- nrow(Factors)
  
  # Number of parameters
  nl <- 2 + Num_Factors
  
  # Compute means of estimated parameters
  btamean <- colMeans(theta_mat[, 3:nl])
  Amean <- colMeans(theta_mat[, 1:2])
  
  # Define the objective function to minimize
  objective_function <- function(alpha) {
    MPME_full(alpha, btamean, Factors, RF, E_f_ret, DivAll, InvAll, rep(1, ncol(DivAll)))
  }
  
  # Find alpha_full that minimizes MPME_full
  result <- optim(par = 0, fn = objective_function, method = "Nelder-Mead")
  alpha_full <- result$par
  eval <- result$value
  
  deltaF <- Amean[1] - alpha_full
  F <- RF + E_g_ret - deltaF
  
  # Compute discount factors
  disc <- rep(1, T + 1)
  for (ci in 2:(T + 1)) {
    disc[ci] <- disc[ci - 1] / (1 + F[ci - 1])
  }
  
  # Compute PME for each fund
  P <- numeric(ncol(DivAll))
  for (ci in 1:ncol(DivAll)) {
    P[ci] <- sum(DivAll[, ci] * disc) / sum(InvAll[, ci] * disc)
  }
  
  # Display results
  cat('PME for full sample (mean, median)\n')
  cat('Mean PME should be one by construction\n')
  mean_P <- mean(P)
  median_P <- median(P)
  print(c(mean_P, median_P))
  
  # The following lines store results in predefined variables
  # Assuming variables 'storg', 'storf', 'storg_std', 'storf_std', 'store', 'wc', 'std_g_ret', 'std_f_ret', 'LL_new', 'Num_Funds' are defined
  # storg[1:T, wc] <- F
  # storf[1:T, wc] <- E_f_ret
  # storg_std[1:T, wc] <- std_g_ret
  # storf_std[1:T, wc] <- std_f_ret
  # store[wc, 1:nl] <- c((1 + alpha_full)^4 - 1, (1 + Amean[1])^4 - 1, btamean)
  # store[wc, 22:25] <- c(mean_P, median_P, LL_new, Num_Funds)
  # store[wc, 26:27] <- Amean
  # a <- c(Num_Factors, sd(theta_mat[, 1:2]), 0, sd(theta_mat[, 3:nl]))
  # store[wc, 10:(9 + length(a))] <- a
  
  return(P)
}


# Fun: LL_allfunds_Norm.m ------

# Log Likelihood for all funds
#
# Moment equation:
# PV(dividends) / PV(investments) ~ lognormal(mu, sigma) with mean 1.
#
# Outputs:
# LL   - Log likelihood
# PME  - Sample PV ratios
# sig  - Standard deviation of the log of PV ratios
# mpme - Mean of the log of PV ratios

LL_allfunds_Norm <- function(div, invest, g_ret, max_sig, mpme = 0) {
  
  # Number of funds and time periods
  Num_Funds <- ncol(div)
  T <- nrow(div) - 1  # CF_mat has one extra row than the time vector
  
  # Replicate g_ret to match dimensions
  r_mat <- matrix(rep(g_ret, Num_Funds), ncol = Num_Funds)
  
  # Compute cumulative returns from time zero
  discount_mat <- rbind(rep(1, Num_Funds), 1 / (1 + r_mat))
  discount_mat <- apply(discount_mat, 2, cumprod)
  
  # Compute Present Values (PVs) of dividends and investments
  PVdiv <- colSums(discount_mat * div)
  PVinvest <- colSums(discount_mat * invest)
  
  # Compute PV ratios (assumed to be log-normally distributed)
  x <- PVdiv / PVinvest  # x is a vector of length Num_Funds
  
  # Demean the log of PV ratios
  # Assume the distribution is centered at 1.
  # The mean of the log-normal distribution is exp(0.5 * sig^2), so mu = 0.5 * sig^2
  sig <- min(sd(log(x)), max_sig + rnorm(1) * 0.01)
  mu <- 0.5 * sig^2
  
  # Construct log likelihood of the log-normal distribution
  LL <- log(1 / (x * sig)) - 0.5 / sig^2 * (log(x) - mu)^2
  LL <- sum(LL)
  
  # Output PV ratios and mean of the log of PV ratios
  PME <- x       # Sample PV ratios (PME)
  mpme <- mean(log(PME))  # Mean of the log of PV ratios
  
  # Return results as a list
  return(list(LL = LL, PME = PME, sig = sig, mpme = mpme))
}



# Step: Draw_g_ret.m -------

# This code assumes that all necessary variables (g_ret, Factors, bta, alpha, RF, max_sig, max_sig_g, max_ret, min_ret, delta, mpme, div, Inv, T, Burn, j, acceptg_ret, std_g_ret, E_g_ret, PME_arr, mpme_arr) are defined elsewhere in your R environment. Also, it assumes that the function LL_allfunds_Norm is properly defined in R and returns a list containing the elements LL, PME, sig, and mpme.

Draw_g_ret <- function() {
  # Store old g_ret in case of no update
  g_ret_old <- g_ret
  
  # Initialize accept vector and step lengths
  accept_vec <- rep(0, T)
  stepLength <- rep(0.05, T)
  sig_g <- max_sig_g / 2  # Could be set to min
  
  # Loop over time periods in random order
  for (t in sample(1:T)) {
    # Compute temporary mean
    temp_mu <- alpha + Factors[t, ] %*% bta
    org_std <- sig_g
    
    # Proposal standard deviation
    prop_std <- stepLength[t]
    
    # Propose new qg
    qg <- temp_mu + prop_std * rnorm(1)
    
    # Update g_ret_new
    g_ret_new <- g_ret
    g_ret_new[t] <- qg
    
    # Adjusted log-likelihood for using prop_std instead of org_std to propose
    LL_Adjust <- -0.5 * (qg - temp_mu)^2 / org_std^2 +
      0.5 * (qg - temp_mu)^2 / prop_std^2 +
      0.5 * (g_ret[t] - temp_mu)^2 / org_std^2 -
      0.5 * (g_ret[t] - temp_mu)^2 / prop_std^2
    
    # Compute log-likelihoods
    LL_old_results <- LL_allfunds_Norm(div, Inv, g_ret + RF, max_sig, mpme)
    LL_old <- LL_old_results$LL
    PME_old <- LL_old_results$PME
    
    LL_new_results <- LL_allfunds_Norm(div, Inv, g_ret_new + RF, max_sig, mpme)
    LL_new <- LL_new_results$LL
    PME_new <- LL_new_results$PME
    
    # Accept/reject
    Likratio <- exp(LL_new - LL_old + LL_Adjust)
    accept_prob <- min(Likratio, 1)
    
    accept <- 0
    if (runif(1) < accept_prob &&
        qg < max_ret && qg > min_ret &&
        abs(mean(PME_new) - 1) < max(delta, abs(mean(PME_old) - 1))) {
      # Accept with probability accept_prob
      accept <- 1
      g_ret[t] <- qg
    }
    
    # Record acceptance
    accept_vec[t] <- accept
  }
  
  # Update acceptance rate after burn-in period
  if (j > Burn) {
    acceptg_ret <- (acceptg_ret * (j - Burn - 1) + accept_vec) / (j - Burn)
  }
  
  # Adjust step lengths periodically
  if (j > Burn && (j %% Burn) == 0) {
    stepLength <- pmin(pmax(std_g_ret, 0.02), 0.5)
  }
  
  # Compute LL, PME, sig, mpme
  LL_results <- LL_allfunds_Norm(div, Inv, g_ret + RF, max_sig, mpme)
  LL <- LL_results$LL
  PME <- LL_results$PME
  sig <- LL_results$sig
  mpme <- LL_results$mpme
  
  # Append PME and mpme to arrays
  PME_arr <- c(PME_arr, PME)
  mpme_arr <- c(mpme_arr, mpme)
  
  # Update E_g_ret and std_g_ret
  if (j <= Burn) {
    iter <- j  # Number of iterations
    prev_x <- E_g_ret  # Previous mean
    prev_x2 <- std_g_ret^2 + E_g_ret^2  # Previous sum of squares
    
    # Update mean
    E_g_ret <- ((iter - 1) * prev_x + g_ret) / iter
    
    # Update sum of squares
    tmp <- ((iter - 1) * prev_x2 + g_ret^2) / iter
    std_g_ret <- sqrt(pmax(tmp - E_g_ret^2, 0))
  }
  
  if (j > Burn) {
    iter <- j - Burn  # Number of iterations after burn-in
    prev_x <- E_g_ret  # Previous mean
    prev_x2 <- std_g_ret^2 + E_g_ret^2  # Previous sum of squares
    
    # Update mean
    E_g_ret <- ((iter - 1) * prev_x + g_ret) / iter
    
    # Update sum of squares
    tmp <- ((iter - 1) * prev_x2 + g_ret^2) / iter
    std_g_ret <- sqrt(pmax(tmp - E_g_ret^2, 0))
  }
}


# Step: Draw_bta.m ------

Draw_bta <- function() {
  # Initialize variables
  TrimPercent <- 1 / T
  
  # Subset Y and X
  Y <- g_ret[5:length(g_ret)]
  X <- cbind(rep(1, T - 4), Factors[5:nrow(Factors), ])
  
  # Trimming extreme values from Y
  indtemp <- (Y <= quantile(Y, probs = 1 - TrimPercent / 2)) & (Y >= quantile(Y, probs = TrimPercent / 2))
  X <- X[indtemp, ]
  Y <- Y[indtemp]
  
  # Prior mean and variance
  MU <- c(p_m_alpha, p_m_bta)
  LAM <- diag(c(p_std_alpha^2, p_std_bta^2))
  
  # Least Absolute Deviations (LAD) estimation loop
  LADLoop <- 3
  E <- diag(nrow(X))
  for (k in 1:LADLoop) {
    BTA <- solve(t(X) %*% E %*% X) %*% (t(X) %*% E %*% Y)
    residuals <- Y - X %*% BTA
    E <- diag(1 / pmax(abs(residuals), 1e-10))
  }
  
  # Updating BTA with prior information
  gp <- cstd / 4
  BTA <- (gp * BTA + MU) / (gp + 1)
  std2 <- solve(t(X) %*% X) * var(Y) * gp / (gp + 1)
  
  # Drawing BTA from a normal distribution
  BTAtmp <- as.numeric(BTA) + t(chol(std2)) %*% rnorm(length(BTA))
  STDtmp <- sqrt(diag(std2))
  
  # Checking bounds for bta
  if (all(BTAtmp[2:(Num_Factors + 1)] < max_bta) && all(BTAtmp[2:(Num_Factors + 1)] > min_bta)) {
    bta <- BTAtmp[2:(Num_Factors + 1)]
  }
  
  # Calculating alpha and f_ret
  alpha <- mean(g_ret[5:length(g_ret)] - Factors[5:nrow(Factors), ] %*% bta)
  f_ret <- g_ret - Factors %*% bta - alpha
  
  # Updating expectations and standard deviations after burn-in period
  if (j > Burn) {
    iter <- j - Burn
    prev_x <- E_f_ret
    prev_x2 <- std_f_ret^2 + E_f_ret^2
    E_f_ret <- ((iter - 1) * prev_x + f_ret) / iter
    tmp <- ((iter - 1) * prev_x2 + f_ret^2) / iter
    std_f_ret <- sqrt(pmax(tmp - E_f_ret^2, 0))
  }
  
  phi <- 0  # Set phi to zero
  
  # Drawing sig_g from a truncated inverse gamma distribution
  u <- g_ret[5:length(g_ret)] - Factors[5:nrow(Factors), ] %*% bta - alpha
  d1 <- p_v5 + sum(u^2)
  p_m5 <- d1 / max_sig_g^2 * set_sig_g
  nu1 <- (T - 4) + p_m5
  
  # Inverse gamma CDF and inverse functions
  invgamcdf <- function(x, shape, rate) {
    pgamma(1 / x, shape = shape, rate = rate, lower.tail = FALSE)
  }
  
  invgaminv <- function(p, shape, rate) {
    1 / qgamma(1 - p, shape = shape, rate = rate)
  }
  
  # Calculating probabilities for truncation
  P_max <- invgamcdf(max_sig_g^2, nu1 / 2, d1 / 2)
  P_min <- invgamcdf(min_sig_g^2, nu1 / 2, d1 / 2)
  P <- runif(1) * (P_max - P_min) + P_min
  
  # Drawing sig_g based on calculated probability
  if (P == 1) {
    sig_g <- max_sig_g
  } else if (P == 0) {
    sig_g <- min_sig_g
  } else {
    sig_g <- sqrt(invgaminv(P, nu1 / 2, d1 / 2))
  }
  
}


# Main: MainCode.m --------

# Code for estimation
# g_t = alpha + phi * g_{t-1} + beta' * F_t + sigma_g * eps_t     Latent
# Observed CFs, unobserved returns
# Important: Data is created in Format.R (assumed equivalent to Format.m)

# Clear workspace
rm(list = ls())

# Load required packages
library(R.matlab)  # For reading .mat files
library(readxl)    # For reading Excel files

# Load data from .mat file (make sure the path is correct)
data <- readMat('C:/Users/lphalippou/Dropbox/PE_Bayesian/Code/NewData.mat')

# Extract variables from the loaded data
CashOut <- data$CashOut
NAV <- data$NAV
CashIn <- data$CashIn
DB <- data$DB

# Add NAVs to CashOut
CashOut <- CashOut + NAV

# Select only vintage 1994-2008
a <- which(DB[3, ] >= 1994 & DB[3, ] <= 2008)
DB <- DB[, a]
CashOut <- CashOut[, a]
CashIn <- CashIn[, a]

# Compute Ncashflows
Ncashflows <- numeric(nrow(CashOut))
for (i in 1:nrow(CashOut)) {
  Ncashflows[i] <- length(which(abs(CashOut[i, ] - CashIn[i, ]) > 0))
}

# Adjust CashIn and CashOut
CashIn[94, ] <- CashIn[94, ] + CashIn[95, ]
CashIn <- CashIn[9:94, ]

CashOut[94, ] <- CashOut[94, ] + CashOut[95, ]
CashOut <- CashOut[9:94, ]

# Add factors from Excel file
FactorM <- read_excel('Factors.xlsx', sheet = 'Sheet1')  # Ensure the file path is correct

# Adjust factors
FactBeg <- (1994 - 1986) * 4 + 2   # Because cash flows are considered at the end of the quarter
FactEnd <- (2015 - 1986) * 4 + 2   # The end is always June 2015

# Subtract Risk-free rate (Rf) from specific factors
for (h in c(8, 20, 28, 29, 38, 39, 40)) {
  FactorM[[h]] <- FactorM[[h]] - FactorM[[12]]
}

# Subtract Vanguard index from specific factors
for (h in c(5, 6, 46)) {
  FactorM[[h]] <- FactorM[[h]] - FactorM[[7]]
}

# Extract Risk-free rate and set T
RF <- FactorM[FactBeg:FactEnd, 12]
T <- (2015 - 1994) * 4 + 1

# Initialize variables
wc <- 1
Burn <- 500
Num_Sim <- 5000

# Main loops over typ, wfactor, and er
# For typ = 1:6 and 10 (6=NatRes, 1=vc, 2=bo, 3=re, 4=credit, 5=FoF, 10=pool)
for (typ in c(1:6, 10)) {
  # wfactor: 1=1F, 3=3F, 4=4F, 5=1F, 8=T1, 9=T3, 10=T4, 11=5FF
  for (wfactor in c(1, 3:5, 8:11)) {
    for (er in 1:10) {
      
      set.seed(er)
      # Call external functions (to be implemented in R)
      # typa()       # Provides Iall set of funds
      # FactorModel()  # Sets up the factor model
      
      # Assuming typa() defines Iall
      # Remove outliers for estimation but keep them for alpha computation
      DBAll <- DB[, Iall]
      DivAll <- CashOut[, Iall]
      InvAll <- CashIn[, Iall]
      p1 <- quantile(DBAll[10, ], probs = 0.05)
      p99 <- quantile(DBAll[10, ], probs = 0.95)
      Irestrict <- which(DBAll[10, ] >= p1 & DBAll[10, ] <= p99)
      div <- DivAll[, Irestrict]          # Dividends
      Inv <- InvAll[, Irestrict]          # Investments
      CF_mat <- div - Inv
      DBin <- DBAll[, Irestrict]
      Num_Funds <- ncol(div)
      
      # Call parameter setup function (to be implemented)
      # Parameter()
      
      # Initialize parameter matrix
      theta_mat <- NULL
      
      # Start burn-in iterations
      for (j in 1:Burn) {
        Draw_g_ret()  # Function to draw g_ret (to be implemented)
        Draw_bta()    # Function to draw beta parameters (to be implemented)
        
        # Update parameter matrix
        thistheta <- c(alpha, phi, t(bta), sig_g)
        theta_mat <- rbind(theta_mat, thistheta)
      }
      
      # Start actual simulations
      theta_mat <- NULL
      PME_arr <- NULL
      
      for (j in (Burn + 1):(Burn + Num_Sim)) {
        Draw_g_ret()  # Function to draw g_ret (to be implemented)
        Draw_bta()    # Function to draw beta parameters (to be implemented)
        
        # Update parameter matrix
        r2_numerator <- var((Factors[2:nrow(Factors), ] - phi * Factors[1:(nrow(Factors) - 1), ]) %*% bta)
        r2_denominator <- var(g_ret[2:length(g_ret)] - phi * g_ret[1:(length(g_ret) - 1)])
        r2 <- r2_numerator / r2_denominator
        
        thistheta <- c(alpha, phi, t(bta), sig_g)
        theta_mat <- rbind(theta_mat, thistheta)
      }
      
      # Display results
      cat('     alpha    phi    bta1    bta2     bta3    sig_g   A   B   sig\n')
      cat('Mean:\n')
      print(colMeans(theta_mat))
      cat('Standard Deviation:\n')
      print(apply(theta_mat, 2, sd))
      
      # Collect results
      kip <- colMeans(theta_mat)
      if (!exists("Beta")) {
        Beta <- matrix(NA, nrow = 1000, ncol = length(kip))  # Adjust the number of rows as needed
      }
      Beta[wc, 1:length(kip)] <- kip
      
      # Call output function to process or save results (to be implemented)
      # output()
      
      wc <- wc + 1
    }
  }
}

# Save results (adjust the object to save as needed)
# save(Beta, file = "Full.RData")

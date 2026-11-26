# ==============================================================================
# Regime 3 Monte Carlo Extension Specification (Drop-In Module for r_project)
# ==============================================================================
# Purpose:
#   Extends the existing mc_step4 simulator (Regime 1 core c_v = 0 and
#   Regime 2/3a vintage shocks M4 W_v) with Regime 3b: an overlapping
#   calendar-time latent PE shock e_u orthogonal to F^pub, plus the
#   leave-band-out cross-vintage pair-product diagnostic U_cross_{0,m}.
#
# Economic link:
#   Each deal invested at month a and exiting at month b accumulates
#     N_{a,b}^{(e)} = exp( sum_{u = a+1}^b e_u - 0.5 * Var(sum_{u = a+1}^b e_u) ),
#   where e_u is stationary zero-mean Gaussian (i.i.d. or monthly AR(1) with
#   persistence phi_e) independent of public factors F^pub, so Condition (X)
#   holds exactly: E[N_{a,b}^{(e)} | F^pub] = 1.
#
# Calibration rule (variance-matched to M4 static vintage shock sigma_W):
#   Choose marginal SD sigma_e so that the within-vintage cross-fund PV
#   covariance matches the static vintage shock variance sigma_W^2
#   (e.g., sigma_W = 0.306 for Buyout M4 2c; sigma_W = 0.50 for VC-scale).
# ==============================================================================

#' Calibrate monthly calendar-time shock SD sigma_e to match target cross-fund
#' log-payoff covariance sigma_W^2 under a deal-timing distribution (a, h).
#'
#' @param sigma_W Target within-vintage shared shock SD (e.g., 0.306 or 0.50).
#' @param I_months Investment period length in months (default 60 = 5 yr).
#' @param h_min_months Minimum holding period in months (default 3 = 1 qtr).
#' @param h_max_months Maximum holding period in months (default 138 = 11.5 yr).
#' @param phi_m Monthly AR(1) persistence of e_u (0 for i.i.d.; e.g., 0.947 -> 0.52/yr).
#' @param n_mc Monte Carlo pairs for computing expected overlap covariance.
#' @return Named list with sigma_e and lookup vector var_sum_h[1..h_max_months].
calibrate_calendar_shock <- function(sigma_W,
                                     I_months = 60L,
                                     h_min_months = 3L,
                                     h_max_months = 138L,
                                     phi_m = 0.0,
                                     n_mc = 5000L) {
  a1 <- sample.int(I_months + 1L, n_mc, replace = TRUE) - 1L
  a2 <- sample.int(I_months + 1L, n_mc, replace = TRUE) - 1L
  h1 <- sample(h_min_months:h_max_months, n_mc, replace = TRUE)
  h2 <- sample(h_min_months:h_max_months, n_mc, replace = TRUE)
  b1 <- a1 + h1
  b2 <- a2 + h2

  if (abs(phi_m) < 1e-12) {
    ov <- pmax(0, pmin(b1, b2) - pmax(a1, a2))
    mean_ov <- mean(ov)
    var_unit <- seq_len(h_max_months)
  } else {
    max_t <- I_months + h_max_months
    idx <- seq_len(max_t)
    K_mat <- phi_m ^ abs(outer(idx, idx, "-"))
    ov <- numeric(n_mc)
    for (k in seq_len(n_mc)) {
      ov[k] <- sum(K_mat[(a1[k] + 1L):b1[k], (a2[k] + 1L):b2[k]])
    }
    mean_ov <- mean(ov)
    var_unit <- vapply(seq_len(h_max_months), function(h) {
      if (h == 1L) return(1.0)
      ks <- seq_len(h - 1L)
      h + 2.0 * sum((h - ks) * (phi_m ^ ks))
    }, numeric(1))
  }

  sigma_e <- sigma_W / sqrt(mean_ov)
  list(
    sigma_e = sigma_e,
    phi_m = phi_m,
    mean_overlap = mean_ov,
    var_sum_h = (sigma_e ^ 2) * var_unit
  )
}

#' Draw a monthly calendar-time latent PE shock path e_1..e_T and return
#' exact mean-one multiplicative deal multipliers N_{a,b}^{(e)}.
#'
#' @param T_months Total calendar span in months.
#' @param cal_spec Output of calibrate_calendar_shock().
#' @param deal_entry_idx Integer vector of calendar entry months a_d (1-based).
#' @param deal_exit_idx Integer vector of calendar exit months b_d (> a_d).
#' @return Numeric vector of mean-one multipliers N_{a_d, b_d}^{(e)}.
draw_calendar_shock_multipliers <- function(T_months,
                                            cal_spec,
                                            deal_entry_idx,
                                            deal_exit_idx) {
  sigma_e <- cal_spec$sigma_e
  phi_m <- cal_spec$phi_m

  eps <- rnorm(T_months)
  e <- numeric(T_months)
  if (abs(phi_m) < 1e-12) {
    e <- sigma_e * eps
  } else {
    e[1L] <- sigma_e * eps[1L]
    s_inn <- sigma_e * sqrt(1.0 - phi_m ^ 2)
    for (t in 2L:T_months) {
      e[t] <- phi_m * e[t - 1L] + s_inn * eps[t]
    }
  }

  cum_e <- c(0.0, cumsum(e))
  sum_e <- cum_e[deal_exit_idx] - cum_e[deal_entry_idx]
  h_d <- deal_exit_idx - deal_entry_idx
  var_e <- cal_spec$var_sum_h[h_d]
  exp(sum_e - 0.5 * var_e)
}

#' Leave-band-out cross-vintage pair-product diagnostic objective U_cross_{H,m}(theta).
#' Included to document Proposition 1(iv): rank-1 collapse and finite-V indefiniteness.
#'
#' @param eps_v Numeric vector of length V containing vintage-level bar_epsilon_{v,H}(theta).
#' @param vintage_years Integer vector of length V of vintage years y_v.
#' @param m Non-negative exclusion band in years (e.g., m = 0, 5, 7, 10).
#' @return Scalar leave-band-out cross-vintage objective value.
u_cross_leave_band_out <- function(eps_v, vintage_years, m = 7L) {
  dist_mat <- abs(outer(vintage_years, vintage_years, "-"))
  mask <- (dist_mat > m)
  n_pairs <- sum(mask)
  if (n_pairs == 0L) stop("Exclusion band m removes all vintage pairs.")
  as.numeric(crossprod(eps_v, (mask * 1.0) %*% eps_v)) / n_pairs
}

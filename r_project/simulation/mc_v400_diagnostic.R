# WP1.3 Diagnostic: Root cause of V = 400 studentization distortion and low N_eff
# Run from r_project/:  LANG=en_US.UTF-8 Rscript simulation/mc_v400_diagnostic.R
#
# Dissects the finite-sample studentization behavior of equal-vintage U_{v,H}
# (q_v = 1, single-factor MKT, H in {0, 15} yr) across:
#   Part 1: Mechanism ablations at V = 400, n = 20 (500 reps each):
#           resampled historical vs. Gaussian vs. low-volatility factor paths;
#           lognormal vs. variance-matched symmetric vs. light symmetric deal
#           shocks; lower sigma (0.05) and shorter deal life (2/2 yr).
#   Part 2: Asymptotic scaling across V in {25, 100, 400, 800, 1600} (n = 20).
#   Part 3: Asymptotic scaling at n = 5 across V in {400, 1600}.
# Also compares the plug-in fixed-b Wald test against the null-evaluated
# centered fixed-b score test (scr0), verifying bit-for-bit equivalence of the
# vectorized sandwich computation with mc.sandwich() before running.
# Writes results/mc_step4/v400_studentization_diagnostic.csv.
suppressMessages({
  source("simulation/mc_driver.R")
  if (!exists("fixedb.crit.t")) source("helper/fixedb.R")
})

D <- mc.load.designs()
pub <- mc.public.factors()
workers <- if (exists("n.parallel.workers")) n.parallel.workers(500, verbose = FALSE) else min(8L, parallel::detectCores(logical = FALSE))

build_fast_direct <- function(vint, n_funds, deals_vint, deals_fund, a, b, mult, cal) {
  X <- as.matrix(cal[, c("RF", "MKT"), drop = FALSE])
  rownames(X) <- NULL
  Tlen <- nrow(X)
  V <- length(vint)
  n_units <- V * n_funds
  v_idx <- match(deals_vint, vint)
  u_idx <- (v_idx - 1L) * n_funds + deals_fund
  k0.v <- as.integer(tapply(a, v_idx, min))
  grids <- data.frame(k0 = k0.v, n = Tlen - k0.v + 1L)
  unit.grid <- rep(seq_len(V), each = n_funds)
  structure(list(
    X = X, grids = grids, unit.grid = unit.grid,
    cf.unit = c(u_idx, u_idx), cf.row = c(a, b), cf.val = c(rep(-1.0, length(a)), mult),
    n.units = n_units, fund.id = as.character(seq_len(n_units)),
    group = unit.grid, n.groups = V,
    group.n = setNames(rep(as.integer(n_funds), V), as.character(vint)),
    vintage = vint
  ), class = "fast.units")
}

fast_u_vec <- function(par, fp, max.month) {
  e <- fast.unit.errors(par, fp, max.month)
  s1 <- rowsum(e, fp$group, reorder = TRUE)[, 1]
  s2 <- rowsum(e^2, fp$group, reorder = TRUE)[, 1]
  n  <- fp$group.n
  (s1^2 - s2) / (n * (n - 1))
}

fast_sandwich_diag <- function(fp, bhat, max.month, delta = 1e-4) {
  g_vint <- fp$vintage
  V <- length(g_vint)
  S <- hac.bandwidth.llsw(V)
  b.rat <- S / V
  Wmat <- outer(g_vint, g_vint, bartlett.weight, bw = S)

  eval_at <- function(b) {
    l0 <- fast_u_vec(c(MKT = b), fp, max.month)
    lp <- fast_u_vec(c(MKT = b + delta), fp, max.month)
    lm <- fast_u_vec(c(MKT = b - delta), fp, max.month)
    g <- (lp - lm) / (2 * delta)
    h <- (lp + lm - 2 * l0) / (delta^2)
    psi_hac <- as.numeric(t(g) %*% Wmat %*% g) / V
    psi_iid <- mean(g^2)
    Hbar <- mean(h)
    se_hac <- sqrt(psi_hac / (V * Hbar^2))
    se_iid <- sqrt(psi_iid / (V * Hbar^2))
    neff <- sum(g^2)^2 / sum(g^4)
    kurt_g <- mean(g^4) / (mean(g^2)^2)
    skew_g <- mean(g^3) / (mean(g^2)^(3/2))
    max_g2 <- max(g^2) / sum(g^2)
    lev <- h / sum(h)
    g_hc2 <- g / pmax(1 - lev, 0.1)
    psi_hc2 <- as.numeric(t(g_hc2) %*% Wmat %*% g_hc2) / V
    se_hc2 <- sqrt(psi_hc2 / (V * Hbar^2))
    gc <- g - mean(g)
    psi_c <- as.numeric(t(gc) %*% Wmat %*% gc) / V
    t_score <- abs(mean(g)) / sqrt(psi_c / V)
    c(se_hac = se_hac, se_iid = se_iid, se_hc2 = se_hc2, psi_hac = psi_hac, psi_iid = psi_iid,
      Hbar = Hbar, neff = neff, kurt_g = kurt_g, skew_g = skew_g, max_g2 = max_g2,
      max_lev = max(lev), g_top_neg = as.numeric(g[which.max(abs(g))] < 0), t_score = t_score)
  }
  a0 <- eval_at(1.0)
  ah <- eval_at(bhat)
  c(bhat = bhat, err = bhat - 1, hac_b = b.rat,
    se_h = ah[["se_hac"]], se_0 = a0[["se_hac"]], se_iid_h = ah[["se_iid"]], se_hc2_h = ah[["se_hc2"]],
    psi_h = ah[["psi_hac"]], psi_0 = a0[["psi_hac"]], psi_iid_0 = a0[["psi_iid"]],
    H_h = ah[["Hbar"]], H_0 = a0[["Hbar"]],
    neff_h = ah[["neff"]], neff_0 = a0[["neff"]],
    kurt_0 = a0[["kurt_g"]], skew_0 = a0[["skew_g"]], max_g2_0 = a0[["max_g2"]],
    max_lev_0 = a0[["max_lev"]], g_top_neg_0 = a0[["g_top_neg"]], t_score0 = a0[["t_score"]])
}

# Verify equivalence against mc.fit / mc.sandwich on v_V400_n5 rep 1
dr1 <- mc.draw(D[["v_V400_n5"]]$design, 1, pub)
fp_orig <- mc.as.fast(dr1, pub, "fund")
fit_orig <- mc.fit(MC.ESTIMATORS$U.eq, fp_orig, "MKT", 180)
v1 <- as.integer(names(dr1$design$n.funds))
fp_fast <- build_fast_direct(v1, 5L, dr1$deals$vintage, dr1$deals$fund, dr1$deals$start, dr1$deals$end, dr1$deals$mult, dr1$public)
chk <- fast_sandwich_diag(fp_fast, fit_orig$MKT, 180)
stopifnot(abs(chk["se_h"] - fit_orig$SE.MKT) < 1e-7, abs(chk["neff_h"] - fit_orig$neff.MKT) < 1e-7)

run_cell <- function(V, n_funds, reps = 500, horizons = c(0, 180),
                     factor_mode = "resample", shock_mode = "lognormal",
                     sigma = 0.130, inv = 5, hold = 11.5, label = "") {
  src <- pub[pub$Date >= as.Date("1967-01-31") & pub$Date <= as.Date("2021-12-31"), ]
  mu_rf <- mean(src$RF); mu_mkt <- mean(src$MKT); sd_mkt <- sd(src$MKT)
  vint <- seq(2005L - V + 1L, 2005L)
  first <- seq(as.Date(sprintf("%d-01-01", min(vint))), as.Date("2021-12-01"), by = "month")
  dates <- as.Date(format(first + 32, "%Y-%m-01")) - 1
  Tlen <- length(dates)
  first.row <- match(as.Date(paste0(vint, "-01-31")), dates)
  q_shift <- (3L - as.integer(format(dates, "%m")) %% 3L) %% 3L
  fund.vint <- rep(vint, each = n_funds)
  fund.first <- rep(first.row, each = n_funds)
  fund.id <- rep(seq_len(n_funds), times = V)
  no.deals <- 15L
  nd <- length(fund.vint) * no.deals
  deal.fund <- rep(seq_along(fund.vint), each = no.deals)
  deals_vint <- fund.vint[deal.fund]
  deals_fund <- fund.id[deal.fund]
  inv_m <- as.integer(12 * inv); hold_m <- as.integer(12 * hold)
  s2 <- log(1 + sigma^2)

  one <- function(s) {
    set.seed(s)
    if (factor_mode == "resample") {
      cal <- src[sample.int(nrow(src), Tlen, replace = TRUE), c("RF", "MKT"), drop = FALSE]
    } else if (factor_mode == "gaussian") {
      cal <- data.frame(RF = rep(mu_rf, Tlen), MKT = pmax(rnorm(Tlen, mu_mkt, sd_mkt), -0.5))
    } else if (factor_mode == "lowvol") {
      cal <- data.frame(RF = rep(mu_rf, Tlen), MKT = rnorm(Tlen, mu_mkt, 0.015))
    }
    lG <- c(0, cumsum(log(1 + cal$RF + cal$MKT)))
    s.month <- fund.first[deal.fund] - 1L + sample.int(inv_m, nd, replace = TRUE)
    e.month <- pmin(s.month + sample.int(hold_m, nd, replace = TRUE), Tlen)
    a <- pmin(pmin(s.month, Tlen) + q_shift[pmin(s.month, Tlen)], Tlen)
    b <- pmin(e.month + q_shift[e.month], Tlen)
    k <- b - a
    if (shock_mode == "lognormal") {
      mult <- exp(lG[b + 1L] - lG[a + 1L] + rnorm(nd, -k * s2 / 2, sqrt(k * s2)))
    } else if (shock_mode == "sym_matched") {
      mult <- exp(lG[b + 1L] - lG[a + 1L]) * (1 + rnorm(nd, 0, sqrt(exp(k * s2) - 1)))
    } else if (shock_mode == "sym_light") {
      mult <- exp(lG[b + 1L] - lG[a + 1L]) * (1 + rnorm(nd, 0, 0.5))
    }
    fp <- build_fast_direct(vint, n_funds, deals_vint, deals_fund, a, b, mult, cal)
    r0 <- min(fp$grids$k0)
    rf_u <- fp$X[r0:Tlen, "RF"]; mk_u <- fp$X[r0:Tlen, "MKT"]
    rhs <- 0.08 - 1 - rf_u
    br <- c(max(c(-1, rhs[mk_u > 0] / mk_u[mk_u > 0])), min(c(4, rhs[mk_u < 0] / mk_u[mk_u < 0])))

    do.call(rbind, lapply(horizons, function(m) {
      f <- function(p) mean(fast_u_vec(c(MKT = p), fp, m))
      ans <- stats::optim(1.0, f, method = "Brent", lower = br[1], upper = br[2])
      bhat <- unname(ans$par)
      if (abs(bhat - br[1]) < 1e-6 || abs(bhat - br[2]) < 1e-6) {
        return(data.frame(rep = s, H = m / 12, at_bound = TRUE, bhat = bhat, err = bhat - 1,
                          hac_b = NA, se_h = NA, se_0 = NA, se_iid_h = NA, se_hc2_h = NA,
                          psi_h = NA, psi_0 = NA, psi_iid_0 = NA, H_h = NA, H_0 = NA,
                          neff_h = NA, neff_0 = NA, kurt_0 = NA, skew_0 = NA, max_g2_0 = NA,
                          max_lev_0 = NA, g_top_neg_0 = NA, t_score0 = NA))
      }
      dg <- fast_sandwich_diag(fp, bhat, m)
      data.frame(rep = s, H = m / 12, at_bound = FALSE, as.list(dg))
    }))
  }

  t0 <- proc.time()[["elapsed"]]
  df <- do.call(rbind, parallel::mclapply(seq_len(reps), one, mc.cores = workers))
  dt <- proc.time()[["elapsed"]] - t0
  do.call(rbind, lapply(unique(df$H), function(h) {
    x <- df[df$H == h & !df$at_bound & is.finite(df$se_h), ]
    err <- x$err; b.rat <- unique(x$hac_b); cr <- fixedb.crit.t(0.05, b.rat)
    th <- err / x$se_h; t0s <- err / x$se_0; tmc <- err / sd(err); tsc <- x$t_score0
    data.frame(
      label = label, V = V, n = n_funds, H = h, reps = nrow(x), b = b.rat, crit = cr,
      sd_err = sd(err), mean_se = mean(x$se_h), se_sd = mean(x$se_h)/sd(err),
      med_se_sd = median(x$se_h)/sd(err), se0_sd = mean(x$se_0)/sd(err),
      rej_h = mean(abs(th) > cr), left_h = sum(th < -cr), right_h = sum(th > cr),
      rej_0 = mean(abs(t0s) > cr), rej_score0 = mean(tsc > cr), rej_mcsd = mean(abs(tmc) > cr),
      corr_err_seh = cor(err, x$se_h), corr_err_se0 = cor(err, x$se_0),
      med_neff_h = median(x$neff_h), med_neff_0 = median(x$neff_0), neff_over_V = median(x$neff_h) / V,
      med_kurt_0 = median(x$kurt_0), med_skew_0 = median(x$skew_0),
      med_max_g2 = median(x$max_g2_0), med_max_lev = median(x$max_lev_0),
      frac_gtop_neg = mean(x$g_top_neg_0), secs = dt, stringsAsFactors = FALSE
    )
  }))
}

all_res <- list(
  run_cell(400, 20, 500, c(0, 180), "resample", "lognormal",   0.130, 5, 11.5, "1. Baseline (resamp+LN)"),
  run_cell(400, 20, 500, c(0, 180), "gaussian", "lognormal",   0.130, 5, 11.5, "2. GaussFac + LN deal"),
  run_cell(400, 20, 500, c(0, 180), "lowvol",   "lognormal",   0.130, 5, 11.5, "3. LowVolFac + LN deal"),
  run_cell(400, 20, 500, c(0, 180), "resample", "sym_matched", 0.130, 5, 11.5, "4. ResampFac + SymMatch"),
  run_cell(400, 20, 500, c(0, 180), "gaussian", "sym_matched", 0.130, 5, 11.5, "5. GaussFac + SymMatch"),
  run_cell(400, 20, 500, c(0, 180), "resample", "sym_light",   0.130, 5, 11.5, "6. ResampFac + LightSym"),
  run_cell(400, 20, 500, c(0, 180), "lowvol",   "sym_light",   0.130, 5, 11.5, "7. LowVolFac + LightSym"),
  run_cell(400, 20, 500, c(0, 180), "resample", "lognormal",   0.050, 5, 11.5, "8. ResampFac + LN sig05"),
  run_cell(400, 20, 500, c(0, 180), "resample", "lognormal",   0.130, 2,  2.0, "9. ResampFac + LN 2/2yr"),
  run_cell(25,   20, 500, c(0, 180), "resample", "lognormal", 0.130, 5, 11.5, "Scale V=25  n=20"),
  run_cell(100,  20, 500, c(0, 180), "resample", "lognormal", 0.130, 5, 11.5, "Scale V=100 n=20"),
  run_cell(800,  20, 500, c(0, 180), "resample", "lognormal", 0.130, 5, 11.5, "Scale V=800 n=20"),
  run_cell(1600, 20, 500, c(0, 180), "resample", "lognormal", 0.130, 5, 11.5, "Scale V=1600 n=20"),
  run_cell(400,  5,  500, c(0, 180), "resample", "lognormal", 0.130, 5, 11.5, "Scale V=400  n=5"),
  run_cell(1600, 5,  500, c(0, 180), "resample", "lognormal", 0.130, 5, 11.5, "Scale V=1600 n=5")
)
out <- do.call(rbind, all_res)
write.csv(out, "results/mc_step4/v400_studentization_diagnostic.csv", row.names = FALSE)
cat("wrote results/mc_step4/v400_studentization_diagnostic.csv\n")

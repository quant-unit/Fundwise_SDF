# Standalone Synthetic Preqin-Calibrated Verification Script
# Run from jbes_replication/:  LANG=en_US.UTF-8 Rscript run_synthetic_verification.R
#
# Verifies end-to-end without proprietary Preqin data:
#   1. Core algebraic identities (horizon rescaling T1, U_{v,H} identity, fixed-b quantiles).
#   2. Synthetic Preqin-calibrated Buyout (25 vintages, 383 funds) and VC (27 vintages, 432 funds)
#      cash-flow panels with calibrated deal volatility and lognormal fund sizes on historical q^5 factors.
#   3. Full estimation & inference pipeline: single-factor U_{v,H} vs. L2 VYP (EW & FW) across
#      H in {0, 5, 15, 30} yr, joint (alpha, beta) at H in {0, 15} yr, calendar-lag vs. index-lag
#      Bartlett fixed-b SEs, effective vintage count N_eff, and centered fixed-b AR 95% sets.
# Writes output/synthetic_verification_summary.csv.

suppressMessages({
  source("simulation/mc_driver.R")
  if (!exists("fixedb.crit.t")) source("helper/fixedb.R")
})

cat("=== JBES Replication Archive: End-to-End Synthetic Verification ===\n")
cat("R version:", R.version.string, "\n\n")

pub <- mc.public.factors()
dir.create("output", showWarnings = FALSE)

# 1. Verify exact algebraic rescaling (T1) and U_{v,H} O(n_v) identity on synthetic draw
n_bo <- mc.data.shape("BO")
n_vc <- mc.data.shape("VC")
cat(sprintf("Loaded public vintage shapes: BO (%d vintages, %d funds), VC (%d vintages, %d funds; %d with n_v>=2)\n",
            length(n_bo), sum(n_bo), length(n_vc), sum(n_vc), sum(n_vc >= 2)))

des_bo <- mc.design(n_bo, theta0 = c(MKT = 1), sigma = 0.130, no.deals = 15,
                    investment.period = 5, max.holding.period = 11.5, sizes = list(sdlog = 1.19))
des_vc <- mc.design(n_vc, theta0 = c(MKT = 1), sigma = 0.154, no.deals = 15,
                    investment.period = 5, max.holding.period = 11.5, sizes = list(sdlog = 0.90))

dr_bo <- mc.draw(des_bo, seed = 202601L, public = pub)
dr_vc <- mc.draw(des_vc, seed = 202602L, public = pub)

# Identity check T1: eps_H(beta) == A_H(beta) * PV(beta)
fp_chk <- mc.as.fast(dr_bo, pub, "vyp")
e0 <- fast.unit.errors(c(MKT = 1.25), fp_chk, 0)
e15 <- fast.unit.errors(c(MKT = 1.25), fp_chk, 180)
sc15 <- fast.horizon.scale(c(MKT = 1.25), fp_chk, 180)
err_t1 <- max(abs(e15 - sc15 * e0))
cat(sprintf("Check T1 (exact horizon rescaling identity): max |eps_H - A_H * PV| = %.3e\n", err_t1))
stopifnot(err_t1 < 1e-12)

# Identity check U_{v,H} (Eq. (6)): the O(n_v) formula used in estimation equals the
# direct sum over distinct fund pairs, for EW and FW fund weights, at H = 15 yr
u_direct <- function(e, fp, s, vintage.weights = "equal") {
  u <- vapply(seq_len(fp$n.groups), function(k) {
    i <- which(fp$group == k)
    if (length(i) < 2) return(NA_real_)
    P <- outer(e[i], e[i]); W <- outer(s[i], s[i])
    (sum(P) - sum(diag(P))) / (sum(W) - sum(diag(W)))
  }, 0)
  w <- u.unit.weights(fp$group.n, vintage.weights)
  ok <- !is.na(w)
  mean(w[ok] * u[ok])
}
err_u <- 0
for (lvl in c("fund", "fund.fw")) {
  fp_u <- mc.as.fast(dr_bo, pub, lvl)
  s_u <- if (lvl == "fund.fw") fp_u$size else rep(1, fp_u$n.units)
  e_u <- fast.unit.errors(c(MKT = 1.25), fp_u, 180)
  q_fast <- fast.u.criterion(e_u, fp_u, s_u)
  q_pair <- u_direct(e_u, fp_u, s_u)
  err_u <- max(err_u, abs(q_fast - q_pair) / max(1, abs(q_pair)))
}
cat(sprintf("Check U_{v,H} identity (O(n_v) formula vs. direct distinct-pair sum, EW and FW): max rel. diff = %.3e\n", err_u))
stopifnot(err_u < 1e-10)

# Identity check fixed-b critical values at empirical V in {23, 25, 27}
c25 <- fixedb.crit.t(0.05, hac.bandwidth.llsw(25) / 25)
c23 <- fixedb.crit.t(0.05, hac.bandwidth.llsw(23) / 23)
c27 <- fixedb.crit.t(0.05, hac.bandwidth.llsw(27) / 27)
cat(sprintf("Check fixed-b 5%% critical values: V=25 -> %.3f, V=23 -> %.3f, V=27 -> %.3f\n\n", c25, c23, c27))
stopifnot(abs(c25 - 2.790) < 0.005, abs(c23 - 2.867) < 0.005, abs(c27 - 2.722) < 0.005)

# Helper: centered fixed-b AR 95% confidence set (Eq. (14))
ar_set_1d <- function(fp, m = 0, alpha = 0.05) {
  br <- feasible.interval(mc.domain(fp, "MKT"))
  z <- fast.horizon.scale(c(MKT = 1), fp, m)
  V <- fp$n.units; S <- hac.bandwidth.llsw(V); b.rat <- S / V; crit <- fixedb.crit.t(alpha, b.rat)
  K <- outer(fp$vintage, fp$vintage, bartlett.weight, bw = S)
  margin <- function(beta) {
    e <- fast.unit.errors(c(MKT = beta), fp, 0)
    if (any(!is.finite(e))) return(NA_real_)
    s <- z * e; gbar <- mean(s); sc <- s - gbar
    om <- sum(K * outer(sc, sc)) / V
    if (!is.finite(om) || om <= 0) return(if (abs(gbar) < 1e-12) -crit else Inf)
    abs(gbar) / sqrt(om / V) - crit
  }
  grid <- seq(br[1], br[2], length.out = 401L); mv <- vapply(grid, margin, 0)
  in.set <- is.finite(mv) & mv <= 0
  runs <- rle(in.set); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1L
  ivs <- lapply(which(runs$values), function(r) {
    i1 <- starts[r]; i2 <- ends[r]
    lo <- if (i1 == 1L) br[1] else uniroot(margin, c(grid[i1 - 1L], grid[i1]), tol = 1e-8)$root
    hi <- if (i2 == length(grid)) br[2] else uniroot(margin, c(grid[i2], grid[i2 + 1L]), tol = 1e-8)$root
    c(lo, hi)
  })
  if (!length(ivs)) return(c(lo = NA_real_, hi = NA_real_))
  c(lo = ivs[[1]][1], hi = ivs[[length(ivs)]][2])
}

# 2. Run full estimation & inference pipeline on synthetic BO & VC panels (EW and FW)
rows <- list()
for (nm in c("BO", "VC")) {
  dr <- if (nm == "BO") dr_bo else dr_vc
  for (w in c("EW", "FW")) {
    fp_fund <- mc.as.fast(dr, pub, if (w == "FW") "fund.fw" else "fund")
    fp_vyp  <- mc.as.fast(dr, pub, if (w == "FW") "vyp.fw" else "vyp")
    est_u   <- if (w == "FW") MC.ESTIMATORS$U.eq.fw else MC.ESTIMATORS$U.eq
    est_l2  <- if (w == "FW") MC.ESTIMATORS$L2.vyp.fw else MC.ESTIMATORS$L2.vyp

    ar0 <- ar_set_1d(fp_vyp, 0)
    for (m in c(0, 60, 180, 360)) {
      fu <- mc.fit(est_u, fp_fund, "MKT", m)
      fl <- mc.fit(est_l2, fp_vyp, "MKT", m)
      p1_u <- if (!fu$at.bound && is.finite(fu$SE.MKT)) fixedb.q1.pvalue(((fu$MKT - 1) / fu$SE.MKT)^2, fu$hac.b) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        panel = nm, weighting = w, model = "MKT", H_yr = m / 12,
        U_beta = fu$MKT, U_se = fu$SE.MKT, U_p1 = p1_u, U_neff = fu$neff.MKT,
        L2_beta = fl$MKT, L2_se = fl$SE.MKT,
        AR_lo = if (m == 0) ar0[["lo"]] else NA_real_,
        AR_hi = if (m == 0) ar0[["hi"]] else NA_real_,
        U_alpha_pct = NA_real_, L2_alpha_pct = NA_real_, stringsAsFactors = FALSE
      )
    }
    for (m in c(0, 180)) {
      fu_ab <- mc.fit(est_u, fp_fund, "Alpha", m)
      fl_ab <- mc.fit(est_l2, fp_vyp, "Alpha", m)
      rows[[length(rows) + 1L]] <- data.frame(
        panel = nm, weighting = w, model = "Alpha", H_yr = m / 12,
        U_beta = fu_ab$MKT, U_se = fu_ab$SE.MKT, U_p1 = NA_real_, U_neff = fu_ab$neff.MKT,
        L2_beta = fl_ab$MKT, L2_se = fl_ab$SE.MKT,
        AR_lo = NA_real_, AR_hi = NA_real_,
        U_alpha_pct = 100 * fu_ab$Coef, L2_alpha_pct = 100 * fl_ab$Coef, stringsAsFactors = FALSE
      )
    }
  }
}

out <- do.call(rbind, rows)
write.csv(out, "output/synthetic_verification_summary.csv", row.names = FALSE)
cat("=== Synthetic Panel Verification Results (true beta = 1, alpha = 0) ===\n")
print(out[out$model == "MKT" & out$H_yr %in% c(0, 15, 30),
          c("panel", "weighting", "H_yr", "U_beta", "U_se", "U_p1", "U_neff", "L2_beta", "AR_lo", "AR_hi")],
      row.names = FALSE, digits = 3)
cat("\nJoint (alpha, beta) synthetic check:\n")
print(out[out$model == "Alpha",
          c("panel", "weighting", "H_yr", "U_beta", "U_alpha_pct", "L2_beta", "L2_alpha_pct")],
      row.names = FALSE, digits = 3)
cat("\nALL VERIFICATION CHECKS PASSED. Output saved to output/synthetic_verification_summary.csv\n")

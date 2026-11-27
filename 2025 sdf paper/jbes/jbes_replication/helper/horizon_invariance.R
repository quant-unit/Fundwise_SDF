# Specification tests for the identifying conditions of U_{v,H} (Theorem T7: factor-path
# exogeneity (X), conditional independence of funds within a vintage, no shared vintage
# shock c_v = 0). Used by results/horizon_invariance_test.R (Preqin data) and
# simulation/mc_horizon_invariance.R (size and power). Needs simulation/mc_driver.R.
#
# 1. Horizon-invariance test (main text). Under T7, U_{v,H} targets theta_0 at every
#    horizon, so beta_{H0} = beta_{H1}. The two single-factor U estimates are the stacked
#    M-estimator of the separable criterion Q_{H0}(b0) + Q_{H1}(b1) over the same vintage
#    units; the contrast b0 - b1 is studentized with the Bartlett sandwich of the stacked
#    per-vintage scores (calendar-distance kernel, S = ceil(1.3 sqrt V), b = S / V, the
#    paper's convention) and compared with the fixed-b distribution (q = 1). Under c_v > 0
#    the target tilts with H (T8, Proposition 2), so the contrast is nonzero.
# 2. Two-moment shared-shock test at H = 0 (Online Supplement). Per vintage with n_v >= 2:
#    g1_v(beta) = vintage mean pricing error bar_eps_v (zero mean at beta_0 under (P)) and
#    g2_v(beta) = U_v(beta) (zero mean at beta_0 under T7, positive under c_v > 0). The
#    centered fixed-b Wald statistic of (g1, g2) (q = 2) is inverted over beta; an empty
#    acceptance set rejects the T7 conditions jointly.

# Vintage units and unit losses of a U estimator (the pipeline's: l_v = w_v U_v)
u.units <- function(est, fp) {
  w <- u.unit.weights(fp$group.n, est$vintage.weights)
  keep <- names(w)[!is.na(w)]
  units <- lapply(keep, function(k) fast.subset(fp, fp$key == k))
  names(units) <- keep
  vint <- vapply(units, function(u) u$vintage[1], 0)
  ord <- order(vint)
  list(units = units[ord], vint = vint[ord], w = w)
}

# Influence contributions H^{-1} g_v (one per vintage) of a single-parameter U fit at
# horizon m, evaluated at par; same derivatives as sandwich.covariance()
u.influence <- function(est, uu, par, m, delta = MC.SETTINGS$delta) {
  loss <- function(p, u) uu$w[[u$key[1]]] * fast.u.criterion(fast.unit.errors(p, u, m), u, mc.u.sizes(est, u), vintage.weights = est$vintage.weights)
  g <- vapply(uu$units, function(u) num.grad(function(p) loss(p, u), par, delta), 0)
  h <- mean(vapply(uu$units, function(u) num.hess(function(p) loss(p, u), par, delta), 0))
  g / h
}

# Horizon-invariance test of beta_{m0} = beta_{m1} (single-factor MKT model).
# f0, f1: mc.fit() rows at m0 and m1 (computed if NULL). Returns NA statistics when
# either fit is on a bound of the admissible domain (no interior asymptotics).
u.horizon.contrast <- function(est, fp, m0 = 0, m1 = 180, f0 = NULL, f1 = NULL) {
  if (is.null(f0)) f0 <- mc.fit(est, fp, "MKT", m0)
  if (is.null(f1)) f1 <- mc.fit(est, fp, "MKT", m1)
  uu <- u.units(est, fp)
  V <- length(uu$units); S <- hac.bandwidth.llsw(V); b <- S / V
  out <- data.frame(m0 = m0, m1 = m1, beta0 = f0$MKT, beta1 = f1$MKT, diff = f0$MKT - f1$MKT,
    se0 = NA_real_, se1 = NA_real_, corr = NA_real_, se.diff = NA_real_, t = NA_real_, p = NA_real_,
    V = V, S = S, b = b, bound = isTRUE(f0$at.bound) || isTRUE(f1$at.bound))
  if (out$bound) return(out)
  i0 <- u.influence(est, uu, c(MKT = f0$MKT), m0)
  i1 <- u.influence(est, uu, c(MKT = f1$MKT), m1)
  K <- outer(uu$vint, uu$vint, bartlett.weight, bw = S)
  q <- function(a, c) as.numeric(t(a) %*% K %*% c) / V^2
  v0 <- q(i0, i0); v1 <- q(i1, i1); v01 <- q(i0, i1); vd <- v0 + v1 - 2 * v01
  out$se0 <- sqrt(v0); out$se1 <- sqrt(v1); out$corr <- v01 / sqrt(v0 * v1)
  if (is.finite(vd) && vd > 0) {
    out$se.diff <- sqrt(vd); out$t <- out$diff / out$se.diff
    out$p <- fixedb.q1.pvalue(out$t^2, b)
  }
  out
}

# Per-vintage moments (g1, g2) at H = 0 for the U-eligible vintages (n_v >= 2)
u.two.moments <- function(est, fp, beta) {
  e <- fast.unit.errors(c(MKT = beta), fp, 0)
  if (any(!is.finite(e))) return(NULL)
  s <- mc.u.sizes(est, fp)
  sum.by <- function(x) { out <- numeric(fp$n.groups); r <- rowsum(x, fp$group, reorder = TRUE); out[as.integer(rownames(r))] <- r[, 1]; out }
  se <- sum.by(e); ss <- sum.by(s)
  g1 <- se / ss
  g2 <- (se^2 - sum.by(e^2)) / (ss^2 - sum.by(s^2))
  ok <- fp$group.n >= 2
  vint <- vapply(seq_len(fp$n.groups), function(k) fp$vintage[fp$group == k][1], 0)
  list(g = cbind(g1 = g1[ok], g2 = g2[ok]), vint = vint[ok])
}

# Centered fixed-b Wald statistic (q = 2) of the two moments at beta, and its p-value
u.two.moment.stat <- function(est, fp, beta) {
  mm <- u.two.moments(est, fp, beta)
  if (is.null(mm)) return(c(W = NA_real_, p = NA_real_))
  g <- mm$g; V <- nrow(g); S <- hac.bandwidth.llsw(V); b <- S / V
  gbar <- colMeans(g); gc <- sweep(g, 2, gbar)
  K <- outer(mm$vint, mm$vint, bartlett.weight, bw = S)
  omega <- t(gc) %*% K %*% gc / V
  W <- tryCatch(as.numeric(V * t(gbar) %*% solve(omega, gbar)), error = function(e) NA_real_)
  c(W = W, p = if (is.finite(W)) fixedb.pvalue(W, b, 2) else NA_real_)
}

# Inversion over a beta grid on the admissible domain: acceptance set at level alpha and
# the largest p-value (the test rejects the T7 conditions when it is below alpha)
u.two.moment.test <- function(est, fp, alpha = 0.05, n.grid = 401L) {
  br <- feasible.interval(mc.domain(fp, "MKT"))
  grid <- seq(br[1], br[2], length.out = n.grid)
  st <- t(vapply(grid, function(bt) u.two.moment.stat(est, fp, bt), c(W = 0, p = 0)))
  acc <- is.finite(st[, "p"]) & st[, "p"] >= alpha
  i.max <- which.max(replace(st[, "p"], !is.finite(st[, "p"]), -1))
  data.frame(p.max = st[i.max, "p"], beta.at.p.max = grid[i.max],
    set.lo = if (any(acc)) min(grid[acc]) else NA_real_, set.hi = if (any(acc)) max(grid[acc]) else NA_real_,
    set.share = mean(acc), set.connected = if (any(acc)) all(diff(which(acc)) == 1L) else NA,
    lo.bound = br[1], hi.bound = br[2])
}

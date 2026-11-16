######################
### Asymptotic inference for the LMD estimators (revision plan B3)
### Sourced by estim_model_optimized.R and helper/test_code_repairs.R
######################
#
# Criterion Q_n(theta) = n^{-1} sum_i l_i(theta) over cross-sectional units i
# (vintage pools or funds; vintages for U_pair), with unit loss
#   L2:      l_i = eps_i(theta)^2
#   U_pair:  l_v = U_{v,H}(theta)
# (the reported objective sqrt(sum eps^2)/n is a monotone transform with the
# same minimizer). Sandwich covariance (Appendix C, prop:clt_H):
#   H   = n^{-1} sum_i nabla^2 l_i(theta_hat)
#   Psi = n^{-1} sum_i sum_j w(v_i, v_j) nabla l_i nabla l_j'   (Bartlett in vintage distance)
#   V   = H^{-1} Psi H^{-1} / n
# Wald test of H0: theta = theta_null:  W = d' V^{-1} d,  d = theta_hat - theta_null.
#
# Repairs relative to the former code in estim_model_optimized.R:
#   (1) derivatives of the unit loss l_i = eps_i^2, not of |eps_i| (err.sqr.calc
#       on a single unit returned sqrt(eps^2) / 1 = |eps|);
#   (2) gradient perturbs one parameter and keeps the others
#       (was p.plus <- par[x] + delta, which dropped all other parameters);
#   (3) off-diagonal Hessian perturbs only parameters i and j
#       (was par + delta, which moved all parameters at once);
#   (4) Wald statistic uses V^{-1} (was V).

num.grad <- function(fun, par, delta = 1e-4) {
  g <- vapply(seq_along(par), function(k) {
    p.plus <- p.minus <- par
    p.plus[k] <- par[k] + delta
    p.minus[k] <- par[k] - delta
    (fun(p.plus) - fun(p.minus)) / (2 * delta)
  }, 0)
  names(g) <- names(par)
  g
}

num.hess <- function(fun, par, delta = 1e-4) {
  k <- length(par)
  H <- matrix(0, k, k, dimnames = list(names(par), names(par)))
  f0 <- fun(par)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (j < i) next
      if (i == j) {
        p.p <- p.m <- par
        p.p[i] <- par[i] + delta
        p.m[i] <- par[i] - delta
        H[i, i] <- (fun(p.p) + fun(p.m) - 2 * f0) / delta^2
      } else {
        pp <- pm <- mp <- mm <- par
        pp[i] <- par[i] + delta; pp[j] <- par[j] + delta
        pm[i] <- par[i] + delta; pm[j] <- par[j] - delta
        mp[i] <- par[i] - delta; mp[j] <- par[j] + delta
        mm[i] <- par[i] - delta; mm[j] <- par[j] - delta
        H[i, j] <- H[j, i] <- (fun(pp) - fun(pm) - fun(mp) + fun(mm)) / (4 * delta^2)
      }
    }
  }
  H
}

bartlett.weight <- function(x, y, bw) {
  d <- abs((x - y) / bw)
  ifelse(d <= 1, 1 - d, 0)
}

# units:     list of unit data (anything unit.loss accepts)
# unit.loss: function(par, unit) -> l_i(par)
# vintages:  vintage year of each unit (for the kernel)
# indep:     TRUE -> Psi uses own-unit products only
sandwich.covariance <- function(units, par, unit.loss, vintages, bw = 12, indep = FALSE,
                                theta.null = NULL, delta = 1e-4) {
  n <- length(units)
  G <- t(vapply(units, function(u) num.grad(function(p) unit.loss(p, u), par, delta), numeric(length(par))))
  if (length(par) == 1) G <- matrix(G, ncol = 1)
  colnames(G) <- names(par)
  W <- if (indep) diag(n) else outer(vintages, vintages, bartlett.weight, bw = bw)
  psi <- t(G) %*% W %*% G / n
  hessian <- Reduce("+", lapply(units, function(u) num.hess(function(p) unit.loss(p, u), par, delta))) / n
  # A singular Hessian (parameters not identified in this sample, e.g. a CV
  # fold with few vintages and several factors) gives NA SEs and Wald test
  # instead of stopping the whole run
  h.inv <- tryCatch(solve(hessian), error = function(e) NULL)
  if (is.null(h.inv)) {
    se <- setNames(rep(NA_real_, length(par)), paste0("SE.", names(par)))
    na.m <- matrix(NA_real_, length(par), length(par), dimnames = list(names(par), names(par)))
    return(list(Psi = psi, Hessian = hessian, CoVM = na.m, par = par, SE = se,
                Wald.stat = NA_real_, Wald.p.value.MKT_1 = NA_real_, se.singular = TRUE))
  }
  covm <- h.inv %*% psi %*% t(h.inv) / n
  dimnames(covm) <- list(names(par), names(par))
  std.err <- sqrt(diag(covm))
  names(std.err) <- paste0("SE.", names(par))
  if (is.null(theta.null)) theta.null <- c(1, rep(0, length(par) - 1)) # H0: beta_MKT = 1, others 0
  d <- par - theta.null
  wald.stat <- tryCatch(as.numeric(t(d) %*% solve(covm) %*% d), error = function(e) NA_real_) # singular V: no Wald test
  list(
    Psi = psi,
    Hessian = hessian,
    CoVM = covm,
    par = par,
    SE = std.err,
    Wald.stat = wald.stat,
    Wald.p.value.MKT_1 = 1 - pchisq(wald.stat, df = length(par)),
    se.singular = FALSE
  )
}

# Zero-bias horizon H0*: continuous-time root of
#   tau_bar(H; r) = H - 1/r + H / (exp(r H) - 1) = T_eff        (main text, Proposition 2(i))
# against the root on the monthly grid tau = 0, 1/12, ..., H of the estimator,
#   sum_k (k/12) G^k / sum_k G^k = T_eff,  G = exp(r / 12),  k = 0, ..., 12 H
# (Online Supplement, proof of Proposition 2(i); H between whole months by linear
# interpolation of the grid average). Grid of Online Supplement Table C.1 (tab:H_star_lookup),
# deterministic columns, r = log(1 + g).
# Run from r_project/:  Rscript simulation/optimal_H_theory/discrete_vs_continuous_H0.R
tau.cont <- function(H, r) if (r == 0) H / 2 else H - 1 / r + H / expm1(r * H)
tau.month <- function(H, r) {
  avg <- function(K) { k <- 0:K; w <- exp(r / 12 * k); sum(k / 12 * w) / sum(w) }
  K <- 12 * H; K0 <- floor(K)
  avg(K0) + (avg(K0 + 1) - avg(K0)) * (K - K0)
}
root <- function(f, T.eff) uniroot(function(H) f(H) - T.eff, c(1e-6, 200), tol = 1e-10)$root
grid <- expand.grid(T.eff = c(3, 5, 7, 8.25, 10), g = c(0, 0.08, 0.10, 0.12, 0.133, 0.15, 0.20))
grid$H.cont <- mapply(function(T.eff, g) root(function(H) tau.cont(H, log1p(g)), T.eff), grid$T.eff, grid$g)
grid$H.month <- mapply(function(T.eff, g) root(function(H) tau.month(H, log1p(g)), T.eff), grid$T.eff, grid$g)
grid$diff.months <- 12 * (grid$H.month - grid$H.cont)
print(format(grid, digits = 4), row.names = FALSE)
cat(sprintf("max |monthly - continuous| = %.2f months\n", max(abs(grid$diff.months))))

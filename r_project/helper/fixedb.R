######################
### Fixed-b inference for the Bartlett (Newey-West) sandwich
### Decision of 17 November 2026 (revision plan, execution log):
###   truncation S = ceil(1.3 * sqrt(V)) (Lazarus, Lewis, Stock and Watson 2018,
###   JBES 36(4)), V = number of vintages in the sample, b = S / V, and
###   Kiefer-Vogelsang (2005, Econometric Theory 21) fixed-b critical values.
######################
#
# Our sandwich weights pairs of units by 1 - |v - w| / bw over vintage years
# (bartlett.weight() in helper/inference.R). This is the Bartlett kernel
# k(j / M) of both papers with M = S = bw, so b = bw / V.
#
# Fixed-b distribution: the limit of the Wald statistic computed with the
# Bartlett HAC at bandwidth M = bT, evaluated at the estimate (residuals sum to
# zero). It is simulated as Kiefer-Vogelsang did: Gaussian location model with
# T = 1000 and the HAC of the demeaned series. For the Bartlett kernel
#   Omega = (2 / (M T)) sum_t S_t S_t' - (1 / (M T)) sum_{t <= T-M} (S_t S_{t+M}' + S_{t+M} S_t'),
# with S_t the partial sums of the demeaned series; this equals the direct
# kernel sum exactly (checked in helper/test_code_repairs.R).
# The t statistic for one restriction has W = t^2.
#
# Two calibrations (decision of 17 November 2026, after the external reviews):
#   q = 1 (all t-tests in the tables): conditional ("Rao-Blackwell") method. In
#     the Gaussian location model the mean is independent of the demeaned
#     series, so with Omega the HAC estimate, P(W > w | Omega) = 2 Phi(-sqrt(w Omega))
#     exactly, and p(w) = E[2 Phi(-sqrt(w Omega))]. Averaging this over simulated
#     Omega has far smaller Monte Carlo error than counting draws of W
#     (critical-value SE about 0.001 with 1e6 draws). Exact b: the distribution
#     depends on M = round(b T); Omega quantiles are stored per M in
#     FIXEDB.Q1.FILE and simulated on first use (fixed seed per M).
#   q >= 2 (joint Wald tests, not printed): quantile table FIXEDB.FILE.
# Estimation rows store the Wald statistic, q and b; p-values are computed at
# analysis time, so a recalibration never invalidates estimation results.

FIXEDB.FILE <- "helper/fixedb_quantiles.rds"
FIXEDB.Q1.FILE <- "helper/fixedb_q1_omega.rds"
FIXEDB.T <- 1000
FIXEDB.Q1.REPS <- 1e6
FIXEDB.Q1.N <- 4000 # stored quantiles of Omega per M (midpoints (i - 0.5) / N)

hac.bandwidth.llsw <- function(V) ceiling(1.3 * sqrt(V))

# Wald statistics of the q-dimensional location model, one per replication
fixedb.simulate <- function(b, q = 1, reps = 1e5, T = 1000, chunk = 5000, seed = 1) {
  set.seed(seed)
  M <- max(1, round(b * T))
  idx1 <- seq_len(T - M)
  idx2 <- M + idx1
  out <- numeric(0)
  for (k in seq_len(ceiling(reps / chunk))) {
    ybar <- matrix(0, q, chunk)
    S <- vector("list", q)
    for (i in seq_len(q)) {
      y <- matrix(rnorm(T * chunk), T, chunk)
      ybar[i, ] <- colMeans(y)
      S[[i]] <- apply(sweep(y, 2, ybar[i, ]), 2, cumsum)
    }
    om <- array(0, c(q, q, chunk))
    for (i in seq_len(q)) for (j in i:q) {
      cross <- colSums(S[[i]][idx1, , drop = FALSE] * S[[j]][idx2, , drop = FALSE]) +
        colSums(S[[j]][idx1, , drop = FALSE] * S[[i]][idx2, , drop = FALSE])
      v <- (2 / (M * T)) * colSums(S[[i]] * S[[j]]) - (1 / (M * T)) * cross
      om[i, j, ] <- v
      om[j, i, ] <- v
    }
    w <- if (q == 1) T * ybar[1, ]^2 / om[1, 1, ] else vapply(seq_len(chunk), function(r) {
      d <- sqrt(T) * ybar[, r]
      as.numeric(t(d) %*% solve(matrix(om[, , r], q, q), d))
    }, 0)
    out <- c(out, w)
  }
  out[seq_len(reps)]
}

# Upper-tail quantiles of W on a grid of tail probabilities, for a grid of b and q
FIXEDB.TAIL <- sort(unique(round(c(0.999, seq(0.99, 0.1, by = -0.01), seq(0.1, 0.01, by = -0.001), 10^seq(-2, -4, length.out = 41)), 10)), decreasing = TRUE)

fixedb.build <- function(b.grid = seq(0.02, 1, by = 0.02), q.grid = 1:4, reps = 1e5, file = FIXEDB.FILE) {
  tab <- array(NA_real_, c(length(b.grid), length(q.grid), length(FIXEDB.TAIL)),
    dimnames = list(b = b.grid, q = q.grid, tail = FIXEDB.TAIL))
  cells <- expand.grid(i = seq_along(b.grid), j = seq_along(q.grid))
  n.cores <- if (exists("n.parallel.workers")) n.parallel.workers(nrow(cells), gb.per.worker = 1.5) else 1
  qs <- parallel::mclapply(seq_len(nrow(cells)), function(k) {
    i <- cells$i[k]; j <- cells$j[k]
    w <- fixedb.simulate(b.grid[i], q.grid[j], reps = reps, seed = 1000 * i + j)
    quantile(w, 1 - FIXEDB.TAIL, names = FALSE, type = 8)
  }, mc.cores = n.cores)
  for (k in seq_len(nrow(cells))) tab[cells$i[k], cells$j[k], ] <- qs[[k]]
  out <- list(tab = tab, b = b.grid, q = q.grid, tail = FIXEDB.TAIL, reps = reps, T = 1000,
    note = "Upper-tail quantiles of the fixed-b Wald statistic, Bartlett kernel; see helper/fixedb.R")
  saveRDS(out, file, compress = "xz")
  invisible(out)
}

.fixedb.env <- new.env()
fixedb.table <- function(file = FIXEDB.FILE) {
  if (is.null(.fixedb.env$tab)) {
    path <- if (file.exists(file)) file else file.path(here::here(), "r_project", file)
    if (!file.exists(path)) path <- here::here(file)
    .fixedb.env$tab <- readRDS(path)
  }
  .fixedb.env$tab
}

# quantiles of W for one q at b, linear interpolation between grid values of b
.fixedb.quantiles <- function(b, q) {
  ft <- fixedb.table()
  if (!is.finite(b) || b <= 0 || b > max(ft$b) + 1e-9 || !(q %in% ft$q)) return(NULL)
  j <- match(q, ft$q)
  if (b <= min(ft$b)) return(ft$tab[1, j, ]) # below the grid: smallest b tabulated (conservative)
  hi <- which(ft$b >= b - 1e-12)[1]
  lo <- max(1, hi - (ft$b[hi] > b + 1e-12))
  if (lo == hi) return(ft$tab[hi, j, ])
  wgt <- (b - ft$b[lo]) / (ft$b[hi] - ft$b[lo])
  (1 - wgt) * ft$tab[lo, j, ] + wgt * ft$tab[hi, j, ]
}

# Upper-tail p-value of a Wald statistic W (q restrictions) at fixed b.
# Two-sided p-value of a t statistic: fixedb.pvalue(t^2, b, 1).
# Returns the smallest tabulated tail probability (1e-4) as an upper bound
# when W lies beyond the table, and NA when b or q is outside the table.
fixedb.pvalue <- function(W, b, q = 1) {
  n <- max(length(W), length(b), length(q))
  W <- rep_len(W, n); b <- rep_len(b, n); q <- rep_len(q, n)
  ft <- fixedb.table()
  vapply(seq_len(n), function(k) {
    if (is.na(W[k]) || is.na(q[k])) return(NA_real_)
    if (q[k] == 1) return(fixedb.q1.pvalue(W[k], b[k])) # conditional method
    qs <- .fixedb.quantiles(b[k], q[k])
    if (is.null(qs)) return(NA_real_)
    if (W[k] <= qs[1]) return(max(ft$tail)) # p above 0.999
    if (W[k] >= qs[length(qs)]) return(min(ft$tail))
    exp(approx(qs, log(ft$tail), xout = W[k])$y)
  }, 0)
}

# Critical value of a t statistic, two-sided level alpha (e.g. 0.05), at fixed b
# (conditional method; fixedb.crit.t.table gives the quantile-table value)
fixedb.crit.t <- function(alpha, b) {
  if (is.null(fixedb.q1.omega(b))) return(NA_real_)
  stats::uniroot(function(c) fixedb.q1.pvalue(c^2, b) - alpha, c(0.5, 50), tol = 1e-10)$root
}
fixedb.crit.t.table <- function(alpha, b) {
  ft <- fixedb.table()
  qs <- .fixedb.quantiles(b, 1)
  if (is.null(qs)) return(NA_real_)
  sqrt(exp(approx(log(ft$tail), log(qs), xout = log(alpha))$y))
}

# ---- q = 1: conditional p-values ----------------------------------------------

# Omega draws (Bartlett HAC of the demeaned Gaussian series, bandwidth M, length T)
fixedb.omega.simulate <- function(M, reps = FIXEDB.Q1.REPS, T = FIXEDB.T, chunk = 5000, seed = 100000 + M) {
  set.seed(seed)
  idx1 <- seq_len(T - M); idx2 <- M + idx1
  out <- numeric(0)
  for (k in seq_len(ceiling(reps / chunk))) {
    y <- matrix(rnorm(T * chunk), T, chunk)
    S <- apply(sweep(y, 2, colMeans(y)), 2, cumsum)
    cross <- if (M < T) 2 * colSums(S[idx1, , drop = FALSE] * S[idx2, , drop = FALSE]) else 0
    out <- c(out, (2 / (M * T)) * colSums(S^2) - (1 / (M * T)) * cross)
  }
  out[seq_len(reps)]
}

fixedb.q1.M <- function(b, T = FIXEDB.T) max(1L, as.integer(round(b * T)))

.fixedb.q1.env <- new.env()
.fixedb.q1.path <- function(file = FIXEDB.Q1.FILE) {
  if (file.exists(file)) return(file)
  alt <- file.path(here::here(), "r_project", file)
  if (file.exists(dirname(alt))) alt else here::here(file)
}
fixedb.q1.store <- function(file = FIXEDB.Q1.FILE) {
  if (is.null(.fixedb.q1.env$store)) {
    path <- .fixedb.q1.path(file)
    .fixedb.q1.env$store <- if (file.exists(path)) readRDS(path) else
      list(T = FIXEDB.T, reps = FIXEDB.Q1.REPS, N = FIXEDB.Q1.N, omega = list(),
           note = "Quantiles of the fixed-b Omega (Bartlett, q = 1) per bandwidth M; see helper/fixedb.R")
  }
  .fixedb.q1.env$store
}

# Simulate and store the Omega quantiles of the given M values (parallel); saves the store
fixedb.q1.precompute <- function(M.values, file = FIXEDB.Q1.FILE, cores = 1) {
  st <- fixedb.q1.store(file)
  todo <- setdiff(as.integer(M.values), as.integer(names(st$omega)))
  if (!length(todo)) return(invisible(st))
  qs <- parallel::mclapply(todo, function(M) {
    stats::quantile(fixedb.omega.simulate(M, st$reps, st$T), (seq_len(st$N) - 0.5) / st$N, names = FALSE, type = 7)
  }, mc.cores = cores)
  for (k in seq_along(todo)) st$omega[[as.character(todo[k])]] <- qs[[k]]
  st$omega <- st$omega[order(as.integer(names(st$omega)))]
  .fixedb.q1.env$store <- st
  saveRDS(st, .fixedb.q1.path(file), compress = "xz")
  invisible(st)
}

fixedb.q1.omega <- function(b) {
  if (!is.finite(b) || b <= 0 || b > 1 + 1e-9) return(NULL)
  M <- fixedb.q1.M(b); st <- fixedb.q1.store()
  if (is.null(st$omega[[as.character(M)]])) st <- fixedb.q1.precompute(M)
  st$omega[[as.character(M)]]
}

# two-sided p-value of t (W = t^2) at fixed b, conditional method
fixedb.q1.pvalue <- function(W, b) {
  n <- max(length(W), length(b)); W <- rep_len(W, n); b <- rep_len(b, n)
  vapply(seq_len(n), function(k) {
    if (is.na(W[k])) return(NA_real_)
    om <- fixedb.q1.omega(b[k]); if (is.null(om)) return(NA_real_)
    mean(2 * stats::pnorm(-sqrt(max(W[k], 0) * om)))
  }, 0)
}

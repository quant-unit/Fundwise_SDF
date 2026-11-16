######################
### Automated verification tests for the revision code repairs (revision plan, Section "Code and data repairs")
### Run from r_project/:  LANG=en_US.UTF-8 Rscript helper/test_code_repairs.R
######################
# Implemented so far:
#   B1  simulator: compounding over (start+1):end, mean-one multiplicative shocks
#   B2  H = 0 convention: origin date not discounted (helper/objective.R)
#   B3  gradient, Hessian, sandwich SEs and Wald test (helper/inference.R)
#   B4  constant-composition VYP pooling  (helper/pool_vintage.R, empirical/prep_preqin.R)
#   B5  common vintage grid (origin t0, data end T.max) and U_{v,H} criterion
#   B6  FW pools normalized by total vintage commitment
#   B7  SR_CV = full-sample estimate / SD_CV in the LaTeX tables
#   B8  common H = 0 held-out validation loss
#   B9  provenance in results and tables; Brent for one-parameter models
#   Sample rule: main sample = complete-history funds; all funds = robustness
#
# Each test prints PASS/FAIL; the script exits non-zero if any test fails.

suppressPackageStartupMessages(library(here))
source(here("helper/pool_vintage.R"))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R")))))
source(here("helper/objective.R"))
source(here("helper/inference.R"))

n.fail <- 0
check <- function(name, ok, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0("  (", detail, ")") else ""))
  if (!isTRUE(ok)) n.fail <<- n.fail + 1
}

# NPV of a dated cash-flow vector under a per-date discount function
npv <- function(dates, cf, disc) sum(cf * disc(dates))

# A family of discount curves: constant annual rates (incl. negative) and
# benchmark curves R^b(beta) = prod(1 + RF + beta * MKT) from the q5 factor file.
make.discount.curves <- function(origin) {
  curves <- lapply(c(-0.05, 0, 0.08, 0.25), function(r) {
    function(d) (1 + r)^(-as.numeric(as.Date(d) - origin) / 365.25)
  })
  q.file <- here("empirical/data_prepared_2026/q_factors.csv")
  if (file.exists(q.file)) {
    q <- read.csv(q.file)
    q$Date <- as.Date(q$Date)
    q <- q[q$Date >= origin, ]
    for (beta in c(0.5, 1, 1.5)) {
      cum <- cumprod(1 + q$RF + beta * q$MKT)
      ym <- format(q$Date, "%Y-%m")
      curves[[length(curves) + 1]] <- local({
        cum <- cum; ym <- ym
        function(d) {
          k <- match(format(as.Date(d), "%Y-%m"), ym)
          ifelse(is.na(k), 1, 1 / cum[k])
        }
      })
    }
  }
  curves
}

# ---------------------------------------------------------------------------
# B4/B6 unit tests on a synthetic vintage (no data needed)
# ---------------------------------------------------------------------------
cat("\n== B4/B6: synthetic vintage ==\n")
d <- as.Date(c("2000-03-31", "2000-06-30", "2000-09-30", "2000-12-31", "2001-03-31"))
toy <- rbind(
  # fund A reports on every date
  data.frame(Fund.ID = "A", Date = d, CF = c(-0.5, -0.3, 0.1, 0.4, 0.9), NAV = c(0.5, 0.8, 0.8, 0.5, 0.2), Size = 100),
  # fund B reports on three dates with a gap and liquidates early
  data.frame(Fund.ID = "B", Date = d[c(1, 3, 4)], CF = c(-0.8, 0.2, 1.1), NAV = c(0.8, 0.7, 0.3), Size = 300),
  # fund C starts late
  data.frame(Fund.ID = "C", Date = d[4:5], CF = c(-1.0, 1.6), NAV = c(1.0, 0.1), Size = 50)
)
toy$CON <- pmin(toy$CF, 0)
toy$DIS <- pmax(toy$CF, 0)
toy$Vintage <- 2000
toy$type <- "BO"
ids <- c("A", "B", "C")
sizes <- c(A = 100, B = 300, C = 50)

for (fsw in c(FALSE, TRUE)) {
  lab <- if (fsw) "FW" else "EW"
  f <- toy
  if (fsw) for (col in c("CF", "NAV", "CON", "DIS")) f[[col]] <- f[[col]] * f$Size
  p <- pool.vintage(f, fund.size.weighting = fsw)
  w <- if (fsw) sizes / sum(sizes) else rep(1 / 3, 3)
  names(w) <- ids

  # pooled CF on each date = sum over all funds / time-invariant denominator
  denom <- if (fsw) sum(sizes) else 3
  expect.cf <- sapply(d, function(t) sum(f$CF[f$Date == t])) / denom
  check(paste(lab, "pooled CF uses fixed denominator"), max(abs(p$CF - expect.cf)) < 1e-15)
  check(paste(lab, "N.funds and Commitment columns"),
        all(p$N.funds == 3) && all(p$Commitment == sum(sizes)))

  # exact NPV identity for every discount curve
  curves <- make.discount.curves(as.Date("1999-12-31"))
  err <- max(sapply(curves, function(disc) {
    pv.pool <- npv(p$Date, p$CF, disc)
    pv.funds <- sapply(ids, function(i) {
      fi <- f[f$Fund.ID == i, ]
      npv(fi$Date, fi$CF / if (fsw) sizes[[i]] else 1, disc)
    })
    abs(pv.pool - sum(w * pv.funds))
  }))
  check(paste(lab, "pooled NPV = weighted avg of fund NPVs (", length(curves), "curves)"),
        err < 1e-14, sprintf("max abs err %.1e", err))

  # NAV: carried forward within each fund's life, zero after liquidation
  nav.raw <- toy$NAV * if (fsw) toy$Size else 1
  nav.B.prev <- nav.raw[toy$Fund.ID == "B" & toy$Date == d[1]] # B has no 2000-06-30 report
  exp.nav.q2 <- (nav.raw[toy$Fund.ID == "A" & toy$Date == d[2]] + nav.B.prev) / denom
  exp.nav.q5 <- (nav.raw[toy$Fund.ID == "A" & toy$Date == d[5]] +
                 nav.raw[toy$Fund.ID == "C" & toy$Date == d[5]]) / denom # B liquidated
  check(paste(lab, "NAV carried forward / zero after liquidation"),
        abs(p$NAV[2] - exp.nav.q2) < 1e-15 && abs(p$NAV[5] - exp.nav.q5) < 1e-15)
}

# the old reporter-mean pooling violates the identity on the same data
old <- aggregate(CF ~ Date, data = toy, mean)
disc <- make.discount.curves(as.Date("1999-12-31"))[[3]]
err.old <- abs(npv(old$Date, old$CF, disc) -
               mean(sapply(ids, function(i) npv(toy$Date[toy$Fund.ID == i], toy$CF[toy$Fund.ID == i], disc))))
check("reporter-mean pooling fails the identity (sanity)", err.old > 1e-3, sprintf("err %.3f", err.old))

# ---------------------------------------------------------------------------
# B4/B6 tests on the rebuilt Preqin files (skipped if not yet generated)
# ---------------------------------------------------------------------------
main.dir <- here("empirical/data_prepared_2026_c")
regions <- c("North America", NA, "Europe")
tol.file <- 1e-12 # CSV round trip keeps 15 significant digits

for (data.dir in c(main.dir, paste0(main.dir, "_allfunds"))) {
is.main <- data.dir == main.dir
cat("\n== B4/B6: Preqin files in", data.dir, "==\n")
abs.sums <- list()
for (region in regions) {
  for (lab in c("EW", "FW")) {
    reg.tag <- if (is.na(region)) "" else paste0("_", region)
    f.fund <- file.path(data.dir, paste0("preqin_cashflows_", lab, reg.tag, "_NAV.csv"))
    f.vyp <- file.path(data.dir, paste0("preqin_cashflows_", lab, "_VYP", reg.tag, "_NAV.csv"))
    rlab <- paste(lab, if (is.na(region)) "ALL" else region)
    if (!file.exists(f.fund) || !file.exists(f.vyp)) {
      # a missing file means prep_preqin.R did not finish: fail, don't skip
      check(paste(rlab, "prepared files exist"), FALSE, basename(f.vyp))
      next
    }
    fund <- read.csv(f.fund)
    vyp <- read.csv(f.vyp)
    fund$Date <- as.Date(fund$Date)
    vyp$Date <- as.Date(vyp$Date)
    fsw <- lab == "FW"

    # sample rule: the main sample holds complete-history funds only, and its
    # pools match the funnel table; the all-funds sample adds truncated funds
    if (is.main) {
      check(paste(rlab, "main sample holds complete-history funds only"),
            all(fund$Complete) &&
              # live cutoff = data end (2022-03-31) minus live.window.months (12)
              all(fund$Age.Last > 5 | as.Date(fund$Last.Date) >= as.Date("2021-03-31")))
    }
    f.funnel <- file.path(main.dir, paste0("preqin_sample_funnel", reg.tag, ".csv"))
    if (file.exists(f.funnel)) {
      fn <- read.csv(f.funnel)
      nv <- vyp[!duplicated(paste(vyp$type, vyp$Vintage)), ]
      k <- match(paste(nv$type, nv$Vintage), paste(fn$type, fn$Vintage))
      expect.n <- if (is.main) fn$N.complete[k] else fn$N.complete[k] + fn$N.truncated[k]
      check(paste(rlab, "pool fund counts match the sample funnel"),
            !anyNA(k) && all(nv$N.funds == expect.n) &&
              all(fn$N.raw == fn$N.complete + fn$N.truncated + fn$N.late.start + fn$N.missing.nav))
    }

    max.err <- 0
    n.checked <- 0
    abs.sum <- c()
    for (key in unique(paste(vyp$type, vyp$Vintage))) {
      tp <- sub(" .*", "", key)
      vin <- as.integer(sub(".* ", "", key))
      pv <- vyp[vyp$type == tp & vyp$Vintage == vin, ]
      fv <- fund[fund$type == tp & fund$Vintage == vin, ]
      sz <- fv$Size[!duplicated(fv$Fund.ID)]
      names(sz) <- fv$Fund.ID[!duplicated(fv$Fund.ID)]
      w <- if (fsw) sz / sum(sz) else rep(1 / length(sz), length(sz))
      names(w) <- names(sz)
      abs.sum[key] <- sum(abs(pv$CF))
      curves <- make.discount.curves(as.Date(paste0(vin - 1, "-12-31")))
      for (disc in curves) {
        pv.pool <- npv(pv$Date, pv$CF, disc)
        pv.funds <- sapply(names(sz), function(i) {
          fi <- fv[fv$Fund.ID == i, ]
          npv(fi$Date, fi$CF / if (fsw) sz[[i]] else 1, disc)
        })
        scale <- max(1, abs(pv.pool))
        max.err <- max(max.err, abs(pv.pool - sum(w * pv.funds)) / scale)
        n.checked <- n.checked + 1
      }
      if (pv$N.funds[1] != length(sz)) max.err <- Inf
    }
    check(paste(rlab, "pooled NPV = weighted avg of fund NPVs"),
          max.err < tol.file,
          sprintf("%d vintage x curve checks, max rel err %.1e", n.checked, max.err))

    # B6: sum of |pooled CF| per vintage is O(1) for both EW and FW, on the
    # estimation sample (BO and VC, seasoned vintages <= 2010, n_v >= 2).
    # Unseasoned 2021 vintages and single-fund niche pools legitimately have
    # tiny totals and are excluded.
    nf <- vyp[!duplicated(paste(vyp$type, vyp$Vintage)), ]
    est <- paste(nf$type, nf$Vintage)[nf$type %in% c("BO", "VC") & nf$Vintage <= 2010 & nf$N.funds >= 2]
    est <- intersect(est, names(abs.sum))
    if (length(est) > 0 && is.main) {
      rng <- range(abs.sum[est])
      check(paste(rlab, "sum |pooled CF| per vintage is O(1) (BO/VC, V<=2010, n_v>=2)"),
            rng[1] > 0.3 && rng[2] < 20,
            sprintf("%d pools, range %.2f to %.2f", length(est), rng[1], rng[2]))
    }
    abs.sums[[rlab]] <- abs.sum[est]
  }
  # B6: FW and EW pools are on the same per-unit-of-commitment scale. Checked
  # on the main sample only: in the all-funds sample, near-empty truncated
  # funds make some small pools legitimately tiny, so the ratio is unstable.
  ew <- abs.sums[[paste("EW", if (is.na(region)) "ALL" else region)]]
  fw <- abs.sums[[paste("FW", if (is.na(region)) "ALL" else region)]]
  if (is.main && length(ew) > 0 && length(fw) > 0) {
    k <- intersect(names(ew), names(fw))
    ratio <- fw[k] / ew[k]
    check(paste(if (is.na(region)) "ALL" else region, "FW/EW scale ratio per vintage in [0.33, 3]"),
          all(ratio > 1 / 3 & ratio < 3),
          sprintf("range %.2f to %.2f", min(ratio), max(ratio)))
  }
}
}

# ---------------------------------------------------------------------------
# B2: horizon convention (synthetic)
# ---------------------------------------------------------------------------
cat("\n== B2: H = 0 convention ==\n")
set.seed(1)
nm <- 120
rf <- rep(0.003, nm)
mkt <- rnorm(nm, 0.006, 0.045)
cf <- c(-0.4, rep(0, 5), -0.3, rep(0, 29), 0.5, rep(0, 23), 0.9, rep(0, nm - 61))
ret.of <- function(beta, alpha = 0) cumprod(1 + rf + beta * mkt + alpha)
r <- ret.of(1.2)
pv.direct <- cf[1] + sum(cf[-1] * r[1] / r[-1])
check("max.month 0 is an alias of max.month 1", horizon.error(cf, r, 0) == horizon.error(cf, r, 1))
check("H = 0 error equals CF_t0 + sum_{t > t0} CF_t / R_{t0,t}",
      abs(horizon.error(cf, r, 1) - pv.direct) < 1e-12, sprintf("diff %.1e", abs(horizon.error(cf, r, 1) - pv.direct)))
lim <- horizon.error(cf, ret.of(1.2, alpha = 1e3), 1)
check("H = 0 error -> CF_t0 as alpha -> infinity (no extra discounting of t0)",
      abs(lim - cf[1]) < 1e-9, sprintf("value at alpha = 1000/month: %.6f", lim))

# ---------------------------------------------------------------------------
# B1: simulator timing and conditionally mean-one deal multipliers
# ---------------------------------------------------------------------------
cat("\n== B1: simulator ==\n")
invisible(capture.output(suppressMessages(source(here("simulation/simulation_multifactor.R")))))
q.sim <- df.q5[as.integer(format(df.q5$Date, "%Y")) >= 1990, ]
G.path <- 1 + 0.001 + q.sim$RF + 1.3 * q.sim$MKT
check("deal multiple at stdv = 0 equals R_{start,end} = prod_{start+1..end} G (and 1 if end = start)",
      abs(deal.multiplier(G.path, 10, 70, 0) - prod(G.path[11:70])) < 1e-15 && deal.multiplier(G.path, 5, 5, 0.2) == 1,
      sprintf("former start:end compounding adds G_start - 1 = %.4f", G.path[10] - 1))
b1.dir <- file.path(tempdir(), "b1_sim")
invisible(capture.output(gen <- create.simulation.multifactor(
  no.deals = 15, investment.period = 5, max.holding.period = 10, alpha = 0.001, beta_MKT = 1.3,
  no.samples = 1, no.funds = 3, min.vin = 1990, max.vin = 1993, stdvs = 0,
  scenario_id = "b1test", output_folder = b1.dir
)))
sim.f <- read.csv(gen$fund_path)
sim.f$Date <- as.Date(sim.f$Date)
pv.sim <- sapply(split(sim.f, sim.f$Fund.ID), function(d) {
  qq <- df.q5[as.integer(format(df.q5$Date, "%Y")) >= d$Vintage[1], ]
  R <- cumprod(1 + 0.001 + qq$RF + 1.3 * qq$MKT)
  sum(d$CF * R[1] / R[match(d$Date, qq$Date)])
})
check("simulated funds at stdv = 0 have PV(theta_0) = 0 (linear SDF)", max(abs(pv.sim)) < 1e-13,
      sprintf("%d funds, max |PV| %.1e", length(pv.sim), max(abs(pv.sim))))
gen.files <- function(nc) {
  d <- file.path(tempdir(), paste0("b1_par_", nc))
  invisible(capture.output(g <- create.simulation.multifactor(
    no.deals = 15, investment.period = 5, max.holding.period = 10, alpha = 0.001, beta_MKT = 1.2,
    no.samples = 6, no.funds = 3, min.vin = 1990, max.vin = 1993, stdvs = 0.3, simulate.mkt = TRUE,
    scenario_id = "par", output_folder = d, n.cores = nc
  )))
  unname(tools::md5sum(c(g$fund_path, g$vyp_path, g$factors_path)))
}
check("parallel sample generation gives byte-identical files to serial generation",
      identical(gen.files(1), gen.files(4)))
set.seed(4)
M <- 20000
R.true <- prod(G.path[11:70])
ratio <- replicate(M, deal.multiplier(G.path, 10, 70, 0.2)) / R.true
z <- (mean(ratio) - 1) / (sd(ratio) / sqrt(M))
check("multiplicative shocks: E[m | factor path] = R_{start,end} at stdv = 0.2 (5-year deal)",
      abs(z) < 3, sprintf("mean m / R = %.4f, z = %.2f, sd = %.2f", mean(ratio), z, sd(ratio)))
add <- replicate(M, deal.multiplier(1, 0, 1, 2, shock.model = "additive"))
check("additive Gaussian shocks with default truncation are biased (stress case only)",
      mean(add) > 1.3, sprintf("E[(G + e)_+] for G = 1, sd = 2: %.3f (theory 1.396)", mean(add)))

# ---------------------------------------------------------------------------
# B3: gradient, Hessian, sandwich covariance and Wald test on a toy NLS model
# ---------------------------------------------------------------------------
cat("\n== B3: inference on a 3-parameter toy model ==\n")
set.seed(3)
nu <- 40
toy.units <- lapply(seq_len(nu), function(i) list(x = rnorm(1), z = runif(1, 0.2, 1.5), y = rnorm(1, 1, 0.5), v = 1980 + i))
eps.toy <- function(p, u) u$y - p[1] - p[2] * u$x - exp(p[3] * u$z)
loss.toy <- function(p, u) eps.toy(p, u)^2
th <- c(MKT = 0.3, b = 0.8, c = 0.4)
grad.an <- function(p, u) { de <- -c(1, u$x, u$z * exp(p[3] * u$z)); 2 * eps.toy(p, u) * de }
hess.an <- function(p, u) {
  de <- -c(1, u$x, u$z * exp(p[3] * u$z))
  d2 <- matrix(0, 3, 3); d2[3, 3] <- -u$z^2 * exp(p[3] * u$z)
  2 * (outer(de, de) + eps.toy(p, u) * d2)
}
# relative error: central differences with delta = 1e-4 have O(delta^2) truncation error
g.err <- max(sapply(toy.units, function(u) { a <- grad.an(th, u); max(abs(num.grad(function(p) loss.toy(p, u), th) - a)) / max(abs(a)) }))
h.err <- max(sapply(toy.units, function(u) { a <- hess.an(th, u); max(abs(num.hess(function(p) loss.toy(p, u), th) - a)) / max(abs(a)) }))
check("numerical gradient of the unit loss matches analytic (all 3 parameters)", g.err < 1e-7, sprintf("max rel err %.1e", g.err))
check("numerical Hessian incl. off-diagonals matches analytic", h.err < 1e-7, sprintf("max rel err %.1e", h.err))
vins <- sapply(toy.units, `[[`, "v")
sw <- sandwich.covariance(toy.units, th, loss.toy, vins, bw = 12, indep = FALSE)
G <- t(sapply(toy.units, grad.an, p = th))
W <- outer(vins, vins, function(a, b) pmax(1 - abs(a - b) / 12, 0))
psi.an <- t(G) %*% W %*% G / nu
H.an <- Reduce("+", lapply(toy.units, hess.an, p = th)) / nu
V.an <- solve(H.an) %*% psi.an %*% solve(H.an) / nu
d <- th - c(1, 0, 0)
W.an <- as.numeric(t(d) %*% solve(V.an) %*% d)
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
check("sandwich covariance (Bartlett, bw = 12) matches analytic", rel(sw$CoVM, V.an) < 1e-7, sprintf("rel err %.1e", rel(sw$CoVM, V.an)))
check("Wald statistic uses V^{-1} and matches analytic", abs(sw$Wald.stat - W.an) / W.an < 1e-7,
      sprintf("W = %.4f vs %.4f; old formula d'Vd = %.2e", sw$Wald.stat, W.an, as.numeric(t(d) %*% V.an %*% d)))
sw.i <- sandwich.covariance(toy.units, th, loss.toy, vins, indep = TRUE)
psi.i <- t(G) %*% G / nu
check("independent-units sandwich matches analytic", rel(sw.i$CoVM, solve(H.an) %*% psi.i %*% solve(H.an) / nu) < 1e-7)
# SE must be evaluated at the estimated horizon, not a fixed one
pv.u <- lapply(1:25, function(i) { set.seed(100 + i); list(cf = c(-1, rep(0, 59), exp(rnorm(1, 0.6, 0.4)), rep(0, 299)), mkt = rnorm(360, 0.007, 0.045), v = 1985 + i) })
l.h <- function(m) function(p, u) horizon.error(u$cf, cumprod(1 + 0.003 + p[1] * u$mkt), m)^2
se.1 <- sandwich.covariance(pv.u, c(MKT = 1.1), l.h(1), sapply(pv.u, `[[`, "v"))$SE
se.180 <- sandwich.covariance(pv.u, c(MKT = 1.1), l.h(180), sapply(pv.u, `[[`, "v"))$SE
check("SE depends on the horizon at which it is evaluated (H = 0 vs 15y)", abs(se.1 / se.180 - 1) > 0.05,
      sprintf("SE(H=0) = %.3f, SE(H=15y) = %.3f", se.1, se.180))

# ---------------------------------------------------------------------------
# B5: rescaling identity eps_H = A_H * PV, A_H from the factor path alone
# ---------------------------------------------------------------------------
cat("\n== B5: horizon rescaling identity (synthetic) ==\n")
err <- max(sapply(c(1, 2, 30, 60, 120, 180), function(m) {
  for (b in c(0.5, 1, 1.8)) {
    rr <- ret.of(b, 0.001)
    a <- mean(rr[1:min(m, nm)] / rr[1])
    e <- abs(horizon.error(cf, rr, m) - a * horizon.error(cf, rr, 1))
    if (e > 1e-12) return(e)
  }
  0
}))
check("eps_H = A_H * PV_0 with A_H = mean_{tau < t0 + m} R_tau / R_t0 (window truncated at T.max)", err < 1e-12,
      sprintf("max abs err %.1e", err))

# ---------------------------------------------------------------------------
# B5: U_{v,H} criterion (synthetic)
# ---------------------------------------------------------------------------
cat("\n== B5: cross-fund pair-product criterion ==\n")
e <- c(0.31, -0.12, 0.05, 0.22)
pairs.ew <- mean(outer(e, e)[upper.tri(diag(4)) | lower.tri(diag(4))])
check("EW U equals mean of cross-fund pair products", abs(u.pair(e) - pairs.ew) < 1e-15)
check("EW U equals ebar^2 - s^2/n", abs(u.pair(e) - (mean(e)^2 - var(e) / length(e))) < 1e-15)
sz <- c(100, 300, 50, 250)
eu <- e # per-unit errors; FW fund-level errors are size * per-unit
off <- upper.tri(diag(4)) | lower.tri(diag(4))
pairs.fw <- sum((outer(sz, sz) * outer(eu, eu))[off]) / sum(outer(sz, sz)[off])
check("FW U equals size-weighted mean of cross-fund pair products", abs(u.pair(sz * eu, sz) - pairs.fw) < 1e-15)
check("vintages with one fund are dropped; n.funds weights q_v = n_v",
      abs(u.criterion(c(e, 0.7, 0.1, 0.2), c(rep("a", 4), "b", "c", "c"), vintage.weights = "n.funds") -
            (4 * u.pair(e) + 2 * u.pair(c(0.1, 0.2))) / 6) < 1e-15)
set.seed(2)
mu <- 0.05; sig <- 0.3; nv <- 3; V <- 40000
sim <- matrix(rnorm(V * nv, mu, sig), V, nv)
u.mc <- apply(sim, 1, u.pair)
l2.mc <- rowMeans(sim)^2
z <- (mean(u.mc) - mu^2) / (sd(u.mc) / sqrt(V))
check("U is unbiased for mu^2 with n_v = 3 (squared pooled error is not)",
      abs(z) < 4 && mean(l2.mc) - mu^2 > 0.9 * sig^2 / nv,
      sprintf("E[U] = %.4f vs mu^2 = %.4f (z = %.1f); E[ebar^2] = %.4f", mean(u.mc), mu^2, z, mean(l2.mc)))

# ---------------------------------------------------------------------------
# B5 on Preqin data: common grid and pooled = mean of fund horizon errors
# ---------------------------------------------------------------------------
cat("\n== B5: common vintage grid on Preqin North America BO (V <= 2010) ==\n")
f.fund <- file.path(main.dir, "preqin_cashflows_EW_North America_NAV.csv")
f.vyp <- file.path(main.dir, "preqin_cashflows_EW_VYP_North America_NAV.csv")
q.file <- here("empirical/data_prepared_2026/q_factors.csv")
if (!file.exists(f.fund) || !file.exists(f.vyp) || !file.exists(q.file)) {
  check("B5 data files exist", FALSE)
} else {
  q <- read.csv(q.file)
  q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
  fund <- read.csv(f.fund); vyp <- read.csv(f.vyp)
  fund <- fund[fund$type == "BO" & fund$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage", "Size")]
  vyp <- vyp[vyp$type == "BO" & vyp$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]
  fund$Date <- as.Date(fund$Date); vyp$Date <- as.Date(vyp$Date)
  t.max <- max(fund$Date)
  for (oc in c("first.cf", "predetermined")) {
    pf <- pad.all.units(fund, oc, t.max)
    pv <- pad.all.units(vyp, oc, t.max)
    rng <- aggregate(Date ~ Fund.ID + Vintage, pf, range)
    same.start <- all(tapply(rng$Date[, 1], rng$Vintage, function(x) length(unique(x)) == 1))
    same.end <- all(rng$Date[, 2] == as.numeric(t.max))
    start.ok <- if (oc == "predetermined") {
      all(as.Date(tapply(rng$Date[, 1], rng$Vintage, min), origin = "1970-01-01") <= as.Date(paste0(sort(unique(rng$Vintage)) - 1, "-12-31")))
    } else TRUE
    check(paste0("[", oc, "] all funds of a vintage share t0; all units end at T.max = ", t.max),
          same.start && same.end && start.ok)
    # pooled (VYP) horizon error = mean of fund horizon errors, every H and beta
    max.err <- 0
    for (vin in unique(pv$Vintage)) {
      pv.v <- pv[pv$Vintage == vin, ]
      g <- 1 + q$RF[match(pv.v$Date, q$Date)] + 1.2 * q$MKT[match(pv.v$Date, q$Date)]
      rr <- cumprod(g)
      fv <- split(pf[pf$Vintage == vin, ], pf$Fund.ID[pf$Vintage == vin], drop = TRUE)
      for (m in c(0, 1, 60, 180, 360)) {
        e.f <- sapply(fv, function(d) horizon.error(d$CF, rr, m))
        e.p <- horizon.error(pv.v$CF, rr, m)
        max.err <- max(max.err, abs(e.p - mean(e.f)) / max(1, abs(e.p)))
      }
    }
    check(paste0("[", oc, "] pooled horizon error = mean of fund horizon errors (H = 0..30y, all vintages)"),
          max.err < 1e-12, sprintf("max rel err %.1e", max.err))
    n.early <- sum(tapply(as.numeric(fund$Date), fund$Vintage, min) <
                     as.numeric(as.Date(paste0(sort(unique(fund$Vintage)) - 1, "-12-31"))))
    if (oc == "predetermined") cat(sprintf("       note: %d BO vintages have cash flows before Dec 31 of vintage - 1; their t0 is the first cash-flow date\n", n.early))
  }
}

# ---------------------------------------------------------------------------
# Speed changes: fast grid padding and prepared units give identical results
# ---------------------------------------------------------------------------
cat("\n== Speed: equivalence with the reference code ==\n")
if (exists("fund") && exists("vyp") && exists("q")) {
  same.pad <- TRUE
  for (dd in list(fund, vyp)) for (oc in c("first.cf", "predetermined")) {
    a <- pad.all.units(dd, oc, t.max, fast = FALSE)
    b <- pad.all.units(dd, oc, t.max, fast = TRUE)
    a <- a[order(a$Fund.ID, a$Date), ]; b <- b[order(b$Fund.ID, b$Date), ]
    cols <- intersect(c("Date", "CF", "type", "Vintage", "Fund.ID", "Size"), colnames(a))
    same.pad <- same.pad && nrow(a) == nrow(b) && all(vapply(cols, function(cc) identical(a[[cc]], b[[cc]]), TRUE))
  }
  check("fast grid padding equals merge() padding (fund-level and VYP, both origins)", same.pad)

  pv <- pad.all.units(vyp, "first.cf", t.max)
  qq <- q; qq$Alpha <- 1
  dfm <- merge(pv, qq, by = "Date", all.x = TRUE)
  prep <- prepare.units(dfm, c("RF", "MKT", "ME", "IA", "ROE", "EG", "Alpha"))
  f1.ref <- function(df.ss, m, par0, model) { # reference: data-frame f1() of estim_model_optimized.R
    if (model == "exp.aff") {
      df.ss$Alpha <- exp(1) - 1
      return(horizon.error(df.ss$CF, exp(cumsum(as.matrix(log(1 + df.ss[, names(par0)])) %*% par0)), m))
    }
    horizon.error(df.ss$CF, exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0)))), m)
  }
  same.obj <- TRUE
  for (model in c("linear", "exp.aff")) for (par0 in list(c(RF = 1, MKT = 1.2), c(RF = 1, MKT = 0.9, Alpha = 0.002), c(RF = 1, MKT = 1.1, ME = 0.3))) {
    for (m in c(0, 1, 120, 360)) {
      ref <- sapply(split(dfm, dfm$Fund.ID), f1.ref, m = m, par0 = par0, model = model)
      new <- vapply(prep$units, unit.horizon.error, 0, max.month = m, par0 = par0, sdf.model = model)
      same.obj <- same.obj && identical(unname(ref[names(new)]), unname(new))
    }
  }
  check("prepared-unit horizon errors are identical to the data-frame f1() (linear and exp.aff, 3 models, 4 horizons)", same.obj)
} else {
  check("speed equivalence data available", FALSE)
}

# ---------------------------------------------------------------------------
# B7/B9: table export (SR_CV definition, provenance header)
# ---------------------------------------------------------------------------
cat("\n== B7/B9: LaTeX table export ==\n")
suppressPackageStartupMessages(library(dplyr))
plot.exprs <- parse(here("results/plot_empirical_estimates.R"))
for (ex in plot.exprs) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
      as.character(ex[[2]]) %in% c(".export_latex_tables", ".table_provenance")) eval(ex, globalenv())
}
toy.tab <- rbind(
  data.frame(Type = "BO", weighting = c("EW", "FW"), method = "Asymptotic", Factor = "MKT", max.month = 1,
             MKT = c(1.21, 1.05), SE.MKT = c(0.07, 0.09), t.MKT = c(17.3, 11.7),
             git.commit = "abc1234", git.dirty = FALSE, scenario.id = "toy", data.prepared.folder = "data_prepared_2026_c"),
  data.frame(Type = "BO", weighting = c("EW", "FW"), method = "Cross-Validation", Factor = "MKT", max.month = 1,
             MKT = c(1.40, 0.80), SE.MKT = c(0.20, 0.25), t.MKT = c(7.0, 3.2),
             git.commit = NA, git.dirty = NA, scenario.id = NA, data.prepared.folder = NA)
)
tab.base <- file.path(tempdir(), "b7_toy")
suppressMessages(.export_latex_tables(toy.tab, "BO", "MKT", tab.base))
tex <- readLines(paste0(tab.base, "_tables.tex"))
row <- tex[grepl("^0.1 yr &", tex)][1] # max.month = 1 is labelled 1/12 = 0.1 yr
cells <- trimws(strsplit(sub("\\\\\\\\.*$", "", row), "&")[[1]])
sr.ew <- as.numeric(cells[6]); sr.fw <- as.numeric(cells[11])
check("SR_CV column = full-sample estimate / SD_CV (not fold mean / SD_CV)",
      isTRUE(all.equal(sr.ew, round(1.21 / 0.20, 2))) && isTRUE(all.equal(sr.fw, round(1.05 / 0.25, 2))),
      sprintf("EW %.2f (fold-mean ratio would be 7.00), FW %.2f (3.20)", sr.ew, sr.fw))
check("table header records generation commit and the results' commit, scenario and data folder",
      any(grepl("^% Generated .* at commit", tex)) && any(grepl("commit\\(s\\) abc1234", tex)) &&
        any(grepl("Scenario\\(s\\): toy", tex)) && any(grepl("data_prepared_2026_c", tex)))

cat("\n== B8/B9: estimator source checks ==\n")
est.src <- paste(readLines(here("estim_model_optimized.R")), collapse = "\n")
check("held-out folds also report the common H = 0 loss (validation.error.H0)",
      grepl("res$validation.error.H0 <- err.sqr.calc(par, max.month = 1", est.src, fixed = TRUE))
check("one-parameter models use Brent's method (parameter name restored); results carry provenance and an at.bound flag",
      grepl('method = "Brent"', est.src, fixed = TRUE) && grepl("err.sqr.calc(setNames(p, par.name), ...)", est.src, fixed = TRUE) &&
        grepl("res$git.commit <- provenance$git.commit", est.src, fixed = TRUE) && grepl("res$at.bound <-", est.src, fixed = TRUE))
f.toy <- function(b) (b - 1.2345678)^2 * (1 + 0.1 * b^2)
br <- optim(c(MKT = 1), f.toy, method = "Brent", lower = -1, upper = 4)$par
check("Brent on the [-1, 4] bracket finds the minimum to 1e-6", abs(br - optimize(f.toy, c(-1, 4), tol = 1e-12)$minimum) < 1e-6)

# ---------------------------------------------------------------------------
# Every Preqin scenario in the empirical config resolves to an existing file,
# in the main sample and in the all-funds robustness sample
# ---------------------------------------------------------------------------
cat("\n== Config: Preqin scenarios resolve to prepared files ==\n")
suppressMessages(source(here("helper/config_loader.R")))
cfg <- tryCatch(load_empirical_config(here("config/empirical_scenarios.yaml")), error = function(e) e)
if (inherits(cfg, "error")) {
  check("empirical config parses (run with LANG=en_US.UTF-8)", FALSE, conditionMessage(cfg))
} else {
  for (folder in c(cfg$defaults$data_prepared_folder, paste0(cfg$defaults$data_prepared_folder, "_allfunds"))) {
    missing <- c()
    n.sc <- 0
    for (id in names(cfg$scenarios)) {
      sc <- merge_with_defaults(cfg$scenarios[[id]], cfg$defaults)
      if (sc$private_source != "preqin") next
      if (sc$data_prepared_folder == cfg$defaults$data_prepared_folder) sc$data_prepared_folder <- folder
      p <- empirical_scenario_to_params(sc)
      nc <- if (!is.null(sc$final_nav_discount) && sc$final_nav_discount != 100) paste0("_NC", sc$final_nav_discount) else ""
      f <- here("empirical", p$data_prepared_folder, paste0("preqin_cashflows_", p$weighting, nc, "_NAV.csv"))
      n.sc <- n.sc + 1
      if (!file.exists(f)) missing <- c(missing, id)
    }
    check(paste("all Preqin scenarios resolve in", folder), length(missing) == 0,
          if (length(missing) == 0) sprintf("%d scenarios", n.sc) else paste("missing:", paste(missing, collapse = ", ")))
  }
}

cat(sprintf("\n%d test(s) failed\n", n.fail))
if (n.fail > 0) quit(status = 1)

######################
### Step 4 Monte Carlo: in-memory data generator
### Revision plan, Section I, "Step 4: scope and compute plan" (decisions of
### 18 November 2026); tests in helper/test_mc.R (T2)
######################
#
# One replication is drawn from a design and a seed and kept in memory; it is
# never written to disk in production (principle (a)). The deal model is the
# repaired generator of simulation/simulation_multifactor.R (B1):
#   - fund j of vintage v holds `no.deals` deals of 1 unit each; a deal starts
#     in month s ~ U{1, ..., 12 * investment.period} of the vintage (month 1 =
#     January of year v) and is held h ~ U{1, ..., 12 * max.holding.period}
#     months, exit capped at the data end T.max (the open position is paid out
#     at T.max, as a terminal NAV at fair value);
#   - deal timing is monthly; cash flows are booked at the quarter end of their
#     month (decision (iii)), and the deal compounds between the booked dates
#     a < b: m = R_{a,b}(theta_0) * Z, R_{a,b} = prod_{u = a+1}^{b} G_u,
#     G_u = 1 + alpha + r_f,u + beta' F_u;
#   - Z = exp(N(-k s^2 / 2, k s^2)) * W_v with k = b - a months and
#     s^2 = log(1 + sigma^2): the product of k i.i.d. mean-one monthly
#     lognormal shocks of the old generator, drawn as one normal (the same
#     distribution, one draw per deal instead of k);
#     W_v = exp(N(-c^2 / 2, c^2)) is a shock shared by all deals of vintage v
#     (design margin M4; c = 0 in the core); independent across vintages, or
#     (shared.rho != 0, DGP (b) of the plan) log W_v = -c^2 / 2 + c u_v with
#     u an AR(1) series in calendar years (unit variance, coefficient rho)
#     sampled at the vintage years, so the marginal law is unchanged and gaps
#     shorten the dependence (as in D0).
# Factor path: the historical monthly factors (the default; the paper
# conditions on them), or (path = list(resample = "iid", from, to)) a path
# drawn for each replication by resampling the historical months from..to
# i.i.d. with replacement, all factors and r_f of a month together, on a
# synthetic monthly calendar from January of the first vintage to T.max
# (needed for more vintages than the history holds). The path is drawn first
# from the replication's seed and returned with the draw.
#   E[m | factor path] = R_{a,b}(theta_0): every fund is priced exactly at
#   theta_0 in expectation (PV_i(theta_0) = 0 when sigma = 0 and c = 0).
# Units: funds (EW, every deal 1 unit) or vintage portfolios
# CF_{v,t} = sum_i CF_{v,i,t} / n_v (constant composition, helper/pool_vintage.R).

# Monthly factors on month ends, as estim_model_optimized.R normalizes them
mc.public.factors <- function(file = "empirical/data_prepared_2026/q_factors.csv") {
  p <- read.csv(file)
  p$Date <- as.Date(p$Date)
  p$Date <- as.Date(format(p$Date, "%Y-%m-01"))
  p$Date <- seq.Date(p$Date[1], by = "month", length.out = nrow(p))
  p$Date <- as.Date(format(p$Date + 32, "%Y-%m-01")) - 1
  p$Alpha <- 1
  p
}

# Fund counts per vintage of the empirical sample (aggregate counts from the
# local fund file; the licensed file itself never leaves the machine)
mc.data.shape <- function(type, max.vintage = 2010,
                          file = "empirical/data_prepared_2026_c/preqin_cashflows_EW_North America_NAV.csv") {
  d <- unique(read.csv(file)[, c("Fund.ID", "type", "Vintage")])
  d <- d[d$type == type & d$Vintage <= max.vintage, ]
  n <- table(d$Vintage)
  setNames(as.integer(n), names(n))
}

# A design: fund counts per vintage (named by year) and the DGP
# (shared.rho, path, path.shock and sizes enter the design only when set, so
# the designs that do not use them keep their specification hash and draws)
#   path.shock = list(a, window.years): a vintage shock correlated with the
#     factor path (log hardening, design A): w_v = a log R_v - log E[R_v^a],
#     R_v the gross benchmark return at theta_0 over the first 12 window.years
#     months from January of the vintage year (fewer where the calendar ends).
#     With i.i.d. resampled months E[R_v^a] = mean(G^a)^K_v exactly, so
#     E[exp(w_v)] = 1: pricing holds unconditionally (P), not given the path
#     (X). Needs a resampled path; uses no random numbers.
#     Caveat (log entry "hardening weak links 1 and 3: results"): the vintage
#     PV discounts investments to the vintage origin with the same path, so
#     this shock also breaks the pricing moment at the origin.
#   path.shock = list(a, level = "deal"): the deal-level version (design A'):
#     each deal's log multiple gets a log R_ab - k log mean(G^a), R_ab the
#     benchmark return over the deal's holding period (k months after the
#     booked investment). With i.i.d. resampled months its mean is one given
#     the information at the investment date, so pricing holds at the
#     vintage origin and fails only given the path (the plan's DGP (c)).
#   sizes = list(sdlog): lognormal fund sizes, drawn after every other draw, so
#     the equal-weighted structures equal those of the same design without sizes.
mc.design <- function(n.funds, theta0 = c(MKT = 1), sigma = 0.2, shared.sd = 0,
                      no.deals = 15, investment.period = 5, max.holding.period = 10,
                      t.max = as.Date("2021-12-31"), quarterly = TRUE, shared.rho = 0, path = NULL,
                      path.shock = NULL, sizes = NULL) {
  stopifnot(!is.null(names(n.funds)), all(n.funds >= 1), all(names(theta0) %in% c("MKT", "Alpha", "ME", "IA", "ROE", "EG")),
            abs(shared.rho) < 1, is.null(path) || identical(path$resample, "iid"),
            is.null(path.shock) || (!is.null(path) && is.finite(path.shock$a) &&
              (identical(path.shock$level, "deal") || (is.null(path.shock$level) && isTRUE(path.shock$window.years > 0)))),
            is.null(sizes) || sizes$sdlog >= 0)
  d <- list(n.funds = n.funds, theta0 = theta0, sigma = sigma, shared.sd = shared.sd, no.deals = no.deals,
       investment.period = investment.period, max.holding.period = max.holding.period,
       t.max = t.max, quarterly = quarterly)
  if (shared.rho != 0) d$shared.rho <- shared.rho
  if (!is.null(path)) d$path <- list(resample = "iid", from = as.Date(path$from), to = as.Date(path$to))
  if (!is.null(path.shock)) d$path.shock <- if (identical(path.shock$level, "deal")) list(a = as.numeric(path.shock$a), level = "deal") else
    list(a = as.numeric(path.shock$a), window.years = as.numeric(path.shock$window.years))
  if (!is.null(sizes)) d$sizes <- list(sdlog = as.numeric(sizes$sdlog))
  d
}

# log E[R^a] per month for the path shock: log mean of G^a over the source
# months of the resampled path, G the gross benchmark return at theta_0
mc.path.shock.logmgf <- function(design, public) {
  src <- public[public$Date >= design$path$from & public$Date <= design$path$to, ]
  par0 <- c(RF = 1, design$theta0)
  G <- as.numeric(1 + as.matrix(src[, names(par0)]) %*% par0)
  log(mean(G^design$path.shock$a))
}

# A factor path for one replication: historical months path$from..path$to
# resampled i.i.d. (whole months), on month ends from January of the first
# vintage to T.max
mc.resample.path <- function(design, public) {
  src <- public[public$Date >= design$path$from & public$Date <= design$path$to, ]
  stopifnot(nrow(src) > 0, !anyNA(src))
  first <- seq(as.Date(sprintf("%d-01-01", min(as.integer(names(design$n.funds))))), as.Date(format(design$t.max, "%Y-%m-01")), by = "month")
  dates <- as.Date(format(first + 32, "%Y-%m-01")) - 1
  out <- src[sample.int(nrow(src), length(dates), replace = TRUE), , drop = FALSE]
  out$Date <- dates
  rownames(out) <- NULL
  out
}

# One replication: the deals (vintage, fund, booked start and end rows of the
# factor calendar, multiple). Reproducible from the seed alone.
mc.draw <- function(design, seed, public) {
  set.seed(seed)
  cal <- if (is.null(design$path)) public[public$Date <= design$t.max, ] else mc.resample.path(design, public)
  T <- nrow(cal)
  par0 <- c(RF = 1, design$theta0)
  lG <- c(0, cumsum(log(as.numeric(1 + as.matrix(cal[, names(par0)]) %*% par0)))) # log R up to each row
  vint <- as.integer(names(design$n.funds))
  first.row <- match(as.Date(paste0(vint, "-01-31")), cal$Date)
  stopifnot(!anyNA(first.row))
  fund.vint <- rep(vint, design$n.funds)
  fund.first <- rep(first.row, design$n.funds)
  fund.id <- sequence(design$n.funds)
  D <- design$no.deals
  nd <- length(fund.vint) * D
  deal.fund <- rep(seq_along(fund.vint), each = D)
  # monthly timing (as the old generator: start uniform over the investment
  # period, holding period uniform), then booking at the quarter end
  inv <- max(1L, as.integer(12 * design$investment.period))
  s.month <- fund.first[deal.fund] - 1L + sample.int(inv, nd, replace = TRUE)
  e.month <- pmin(s.month + sample.int(as.integer(12 * design$max.holding.period), nd, replace = TRUE), T)
  book <- function(row) if (design$quarterly) pmin(row + (3L - as.integer(format(cal$Date[row], "%m")) %% 3L) %% 3L, T) else row
  a <- book(pmin(s.month, T)); b <- book(e.month)
  k <- b - a
  s2 <- log(1 + design$sigma^2)
  z <- if (design$sigma > 0) rnorm(nd, -k * s2 / 2, sqrt(k * s2)) else numeric(nd)
  w <- if (design$shared.sd == 0) numeric(length(vint)) else if (is.null(design$shared.rho)) {
    rnorm(length(vint), -design$shared.sd^2 / 2, design$shared.sd)
  } else {
    span <- min(vint):max(vint); u <- rnorm(length(span)); r <- design$shared.rho
    for (t in seq_along(span)[-1]) u[t] <- r * u[t - 1] + sqrt(1 - r^2) * u[t]
    -design$shared.sd^2 / 2 + design$shared.sd * u[match(vint, span)]
  }
  wd <- numeric(nd) # deal-level path shock (design A')
  if (!is.null(design$path.shock) && identical(design$path.shock$level, "deal")) {
    wd <- design$path.shock$a * (lG[b + 1L] - lG[a + 1L]) - k * mc.path.shock.logmgf(design, public)
  } else if (!is.null(design$path.shock)) { # w_v = a log R_v - K_v log mean(G^a)
    K <- pmin(as.integer(12 * design$path.shock$window.years), T - first.row + 1L)
    w <- w + design$path.shock$a * (lG[first.row + K] - lG[first.row]) - K * mc.path.shock.logmgf(design, public)
  }
  mult <- exp(lG[b + 1L] - lG[a + 1L] + z + w[match(fund.vint[deal.fund], vint)] + wd)
  out <- list(design = design, seed = seed, dates = cal$Date,
       deals = data.frame(vintage = fund.vint[deal.fund], fund = fund.id[deal.fund], start = a, end = b, mult = mult))
  if (!is.null(design$path)) out$public <- cal # the replication's own factor path
  if (!is.null(design$sizes)) out$sizes <- data.frame(vintage = fund.vint, fund = fund.id, size = exp(rnorm(length(fund.vint), 0, design$sizes$sdlog)))
  out
}

# Fund sizes of the deals' funds (1 without sizes)
mc.deal.size <- function(draw) {
  if (is.null(draw$sizes)) return(rep(1, nrow(draw$deals)))
  draw$sizes$size[match(paste(draw$deals$vintage, draw$deals$fund), paste(draw$sizes$vintage, draw$sizes$fund))]
}

# Cash flows of a replication as rows (Date, CF, type, Vintage, Fund.ID), for
# the pipeline (tests T1, T3); level "fund" or "vyp" (constant composition)
# fw = TRUE: fund-size weighting as the pipeline's FW files: fund cash flows
# scaled by the fund's size (with a Size column), vintage pools
# sum_i Size_i cf_i / sum_i Size_i (pool.vintage, fund.size.weighting = TRUE)
mc.as.rows <- function(draw, level = c("fund", "vyp"), type = "SIM", fw = FALSE) {
  level <- match.arg(level)
  d <- draw$deals
  sz <- if (fw) mc.deal.size(draw) else rep(1, nrow(d))
  flows <- data.frame(vintage = c(d$vintage, d$vintage), fund = c(d$fund, d$fund),
                      row = c(d$start, d$end), CF = c(-sz, d$mult * sz))
  if (level == "vyp") {
    den <- if (fw) tapply(draw$sizes$size, draw$sizes$vintage, sum) else draw$design$n.funds
    flows$CF <- flows$CF / as.numeric(den[as.character(flows$vintage)])
    flows$fund <- 0L
  }
  agg <- aggregate(CF ~ vintage + fund + row, flows, sum)
  # a zero flow at T.max for every unit: the pipeline ends the grid at the
  # latest date in the file (grid.t.max), which must be the design's T.max
  T <- length(draw$dates)
  units <- unique(agg[, c("vintage", "fund")])
  agg <- rbind(agg, data.frame(units, row = T, CF = 0))
  agg <- aggregate(CF ~ vintage + fund + row, agg, sum)
  out <- data.frame(Date = draw$dates[agg$row], CF = agg$CF, type = type, Vintage = agg$vintage,
             Fund.ID = if (level == "vyp") paste0(agg$vintage, "_", type) else paste0(agg$vintage, "_", agg$fund, "_", type),
             stringsAsFactors = FALSE)
  if (fw && level == "fund") out$Size <- draw$sizes$size[match(paste(agg$vintage, agg$fund), paste(draw$sizes$vintage, draw$sizes$fund))]
  out
}

# The fast-criterion structure of a replication directly (production path):
# the same as prepare.fast() on the padded rows of mc.as.rows() (test T2),
# without building and padding data frames. Origin of a vintage: its first
# booked cash flow (origin.convention "first.cf"); grid end: T.max.
# Levels "fund.fw" and "vyp.fw": fund-size weighting as mc.as.rows(fw = TRUE);
# at "fund.fw" the units carry their sizes ($size), which the FW U criterion
# uses in its pair normalization (as the pipeline, use.size)
mc.as.fast <- function(draw, public, level = c("fund", "vyp", "fund.fw", "vyp.fw"), type = "SIM",
                       factor.cols = unique(c("RF", setdiff(colnames(public), "Date")))) {
  level <- match.arg(level)
  fw <- level %in% c("fund.fw", "vyp.fw"); level <- sub("\\.fw$", "", level)
  if (fw && is.null(draw$sizes)) stop("fund-size weighting needs a design with sizes")
  if (!is.null(draw$public)) public <- draw$public # a resampled path
  d <- draw$deals
  T <- length(draw$dates)
  X <- as.matrix(public[match(draw$dates, public$Date), factor.cols, drop = FALSE])
  rownames(X) <- NULL
  vint <- as.integer(names(draw$design$n.funds))
  k0.v <- vapply(vint, function(v) min(d$start[d$vintage == v]), 0L)
  unit.of <- if (level == "fund") paste0(d$vintage, "_", d$fund, "_", type) else paste0(d$vintage, "_", type)
  ids <- sort(unique(unit.of)) # split() order of prepare.units()
  sz <- if (fw) mc.deal.size(draw) else rep(1, nrow(d))
  scale <- if (level == "vyp") {
    den <- if (fw) tapply(draw$sizes$size, draw$sizes$vintage, sum) else draw$design$n.funds
    sz / as.numeric(den[as.character(d$vintage)])
  } else sz
  flows <- data.frame(unit = match(c(unit.of, unit.of), ids), row = c(d$start, d$end),
                      CF = c(-scale, d$mult * scale), vintage = c(d$vintage, d$vintage))
  agg <- aggregate(CF ~ unit + row + vintage, flows, sum)
  agg <- agg[agg$CF != 0, ]
  unit.vint <- vapply(seq_along(ids), function(j) as.integer(sub("_.*", "", ids[j])), 0L)
  grids <- data.frame(k0 = k0.v, n = T - k0.v + 1L)
  unit.grid <- match(unit.vint, vint)
  o <- order(agg$unit, agg$row)
  agg <- agg[o, ]
  key <- paste(type, unit.vint)
  groups <- sort(unique(key)); group <- match(key, groups)
  structure(list(X = X, grids = grids, unit.grid = unit.grid,
    cf.unit = agg$unit, cf.pos = agg$row - k0.v[unit.grid[agg$unit]] + 1L, cf.row = agg$row, cf.val = agg$CF,
    n.units = length(ids), fund.id = ids, key = setNames(key, ids), group = group, n.groups = length(groups),
    group.n = setNames(tabulate(group, length(groups)), groups),
    vintage = setNames(as.numeric(unit.vint), ids),
    size = setNames(if (fw && level == "fund") vapply(ids, function(u) sz[match(u, unit.of)], 0) else rep(1, length(ids)), ids)),
    class = "fast.units")
}

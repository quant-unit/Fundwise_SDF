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
#     (design margin M4; c = 0 in the core).
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
mc.design <- function(n.funds, theta0 = c(MKT = 1), sigma = 0.2, shared.sd = 0,
                      no.deals = 15, investment.period = 5, max.holding.period = 10,
                      t.max = as.Date("2021-12-31"), quarterly = TRUE) {
  stopifnot(!is.null(names(n.funds)), all(n.funds >= 1), all(names(theta0) %in% c("MKT", "Alpha", "ME", "IA", "ROE", "EG")))
  list(n.funds = n.funds, theta0 = theta0, sigma = sigma, shared.sd = shared.sd, no.deals = no.deals,
       investment.period = investment.period, max.holding.period = max.holding.period,
       t.max = t.max, quarterly = quarterly)
}

# One replication: the deals (vintage, fund, booked start and end rows of the
# factor calendar, multiple). Reproducible from the seed alone.
mc.draw <- function(design, seed, public) {
  set.seed(seed)
  cal <- public[public$Date <= design$t.max, ]
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
  w <- if (design$shared.sd > 0) rnorm(length(vint), -design$shared.sd^2 / 2, design$shared.sd) else numeric(length(vint))
  mult <- exp(lG[b + 1L] - lG[a + 1L] + z + w[match(fund.vint[deal.fund], vint)])
  list(design = design, seed = seed, dates = cal$Date,
       deals = data.frame(vintage = fund.vint[deal.fund], fund = fund.id[deal.fund], start = a, end = b, mult = mult))
}

# Cash flows of a replication as rows (Date, CF, type, Vintage, Fund.ID), for
# the pipeline (tests T1, T3); level "fund" or "vyp" (constant composition)
mc.as.rows <- function(draw, level = c("fund", "vyp"), type = "SIM") {
  level <- match.arg(level)
  d <- draw$deals
  flows <- data.frame(vintage = c(d$vintage, d$vintage), fund = c(d$fund, d$fund),
                      row = c(d$start, d$end), CF = c(rep(-1, nrow(d)), d$mult))
  if (level == "vyp") {
    n <- draw$design$n.funds
    flows$CF <- flows$CF / as.numeric(n[as.character(flows$vintage)])
    flows$fund <- 0L
  }
  agg <- aggregate(CF ~ vintage + fund + row, flows, sum)
  # a zero flow at T.max for every unit: the pipeline ends the grid at the
  # latest date in the file (grid.t.max), which must be the design's T.max
  T <- length(draw$dates)
  units <- unique(agg[, c("vintage", "fund")])
  agg <- rbind(agg, data.frame(units, row = T, CF = 0))
  agg <- aggregate(CF ~ vintage + fund + row, agg, sum)
  data.frame(Date = draw$dates[agg$row], CF = agg$CF, type = type, Vintage = agg$vintage,
             Fund.ID = if (level == "vyp") paste0(agg$vintage, "_", type) else paste0(agg$vintage, "_", agg$fund, "_", type),
             stringsAsFactors = FALSE)
}

# The fast-criterion structure of a replication directly (production path):
# the same as prepare.fast() on the padded rows of mc.as.rows() (test T2),
# without building and padding data frames. Origin of a vintage: its first
# booked cash flow (origin.convention "first.cf"); grid end: T.max.
mc.as.fast <- function(draw, public, level = c("fund", "vyp"), type = "SIM",
                       factor.cols = unique(c("RF", setdiff(colnames(public), "Date")))) {
  level <- match.arg(level)
  d <- draw$deals
  T <- length(draw$dates)
  X <- as.matrix(public[match(draw$dates, public$Date), factor.cols, drop = FALSE])
  rownames(X) <- NULL
  vint <- as.integer(names(draw$design$n.funds))
  k0.v <- vapply(vint, function(v) min(d$start[d$vintage == v]), 0L)
  unit.of <- if (level == "fund") paste0(d$vintage, "_", d$fund, "_", type) else paste0(d$vintage, "_", type)
  ids <- sort(unique(unit.of)) # split() order of prepare.units()
  scale <- if (level == "vyp") 1 / as.numeric(draw$design$n.funds[as.character(d$vintage)]) else rep(1, nrow(d))
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
    vintage = setNames(as.numeric(unit.vint), ids), size = setNames(rep(1, length(ids)), ids)),
    class = "fast.units")
}

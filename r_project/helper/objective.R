######################
### Estimation objective building blocks (revision plan B2, B5)
### Sourced by estim_model_optimized.R and helper/test_code_repairs.R
######################
#
# Horizon convention (B2, assumption A5)
#   Every unit (vintage pool or fund) is placed on a monthly grid that starts at
#   the vintage origin t0 and ends at the common data end T.max. With gross
#   benchmark returns G_u(theta) and R_t = prod_{u <= t} G_u on that grid,
#     PV(theta)        = sum_t CF_t * R_{t0} / R_t         (origin not discounted)
#     eps_H(theta)     = PV(theta) * A_H(theta),
#     A_H(theta)       = mean_{tau = t0, ..., t0 + m - 1} R_tau / R_{t0},
#   where m = max.month is the number of monthly dates averaged (m = 1 is the
#   inception NPV, H = 0; a horizon of H years averages m = 12 H dates starting
#   at t0, truncated at T.max). max.month = 0 is an alias of max.month = 1: the
#   former code discounted the origin date once more (PV / G_{t0}), which tends
#   to 0 as alpha -> infinity.
#
# Grid convention (B5, assumption A5)
#   t0 is identical for all units of a vintage (same type and Vintage):
#     origin.convention = "first.cf"       : first cash-flow date of the vintage
#     origin.convention = "predetermined"  : Dec 31 of year Vintage - 1 (or the
#                                            first cash-flow date if earlier)
#   T.max is the same for all units (the data end), in the empirical and the
#   simulation branch alike. A_H(theta) is therefore a function of the factor
#   path and t0 only and factors out of within-vintage sums and pair products.

month.end.grid <- function(t0, t.max) seq(t0 + 1, t.max + 1, by = "month") - 1

# Vintage origin t0 per unit key (type, Vintage)
vintage.origins <- function(df, origin.convention = "first.cf") {
  key <- paste(df$type, df$Vintage)
  first.cf <- as.Date(tapply(as.numeric(df$Date), key, min), origin = "1970-01-01")
  if (origin.convention == "first.cf") return(first.cf)
  if (origin.convention == "predetermined") {
    vin <- as.integer(sub(".* ", "", names(first.cf)))
    predet <- as.Date(paste0(vin - 1, "-12-31"))
    out <- as.Date(pmin(as.numeric(predet), as.numeric(first.cf)), origin = "1970-01-01")
    names(out) <- names(first.cf)
    return(out)
  }
  stop("unknown origin.convention: ", origin.convention)
}

# Pad one unit's cash flows with zeros on the grid t0, ..., t.max
to.monthly.grid <- function(df.ss, t0, t.max, keep.cols = c("type", "Vintage", "Fund.ID", "Size")) {
  df.m <- data.frame(Date = month.end.grid(t0, t.max))
  df.ss <- merge(df.ss, df.m, by = "Date", all = TRUE)
  df.ss$CF[is.na(df.ss$CF)] <- 0
  for (col in intersect(keep.cols, colnames(df.ss))) {
    df.ss[, col] <- df.ss[which(!is.na(df.ss[, col]))[1], col]
  }
  df.ss[order(df.ss$Date), ]
}

# Same result as to.monthly.grid() without merge(): rows are placed on the
# grid by index (speed). Falls back to merge() if a unit has duplicate dates
# or dates off the month-end grid. Equivalence: helper/test_code_repairs.R.
to.monthly.grid.fast <- function(df.ss, t0, t.max, keep.cols = c("type", "Vintage", "Fund.ID", "Size")) {
  grid <- month.end.grid(t0, t.max)
  idx <- match(as.numeric(df.ss$Date), as.numeric(grid)) # numeric match: Date match() goes via strings
  if (anyNA(idx) || anyDuplicated(df.ss$Date)) return(to.monthly.grid(df.ss, t0, t.max, keep.cols))
  out <- df.ss[match(seq_along(grid), idx), , drop = FALSE]
  out$Date <- grid
  out$CF[is.na(out$CF)] <- 0
  for (col in intersect(keep.cols, colnames(out))) out[[col]] <- df.ss[[col]][1]
  rownames(out) <- NULL
  out
}

pad.all.units <- function(df, origin.convention = "first.cf", t.max = max(df$Date), fast = TRUE) {
  t0 <- vintage.origins(df, origin.convention)
  units <- split(df, df$Fund.ID, drop = TRUE)
  pad <- if (fast) to.monthly.grid.fast else to.monthly.grid
  out <- lapply(units, function(d) pad(d, t0[[paste(d$type[1], d$Vintage[1])]], t.max))
  as.data.frame(data.table::rbindlist(out, use.names = TRUE, fill = TRUE))
}

# ---------------------------------------------------------------------------
# Prepared units (speed): the objective is evaluated hundreds of times per
# estimate, so each unit's cash flows and factor columns are extracted once
# instead of re-splitting the data frame on every evaluation. Results are
# identical to the data-frame path (equivalence: helper/test_code_repairs.R).
# ---------------------------------------------------------------------------
prepare.units <- function(df, factor.cols, split.key = df$Fund.ID) {
  if (inherits(df, "prepared.units")) return(df)
  factor.cols <- intersect(factor.cols, colnames(df))
  # split row indices, not the data frame (split.data.frame is slow)
  X.all <- as.matrix(df[, factor.cols, drop = FALSE])
  cf.all <- df$CF
  date.all <- as.numeric(df$Date)
  has.size <- "Size" %in% colnames(df)
  rows <- split(seq_len(nrow(df)), split.key, drop = TRUE)
  units <- lapply(rows, function(i) {
    i <- i[order(date.all[i])]
    list(
      cf = cf.all[i],
      X = X.all[i, , drop = FALSE],
      fund.id = as.character(df$Fund.ID[i[1]]),
      key = paste(df$type[i[1]], df$Vintage[i[1]]),
      vintage = as.numeric(df$Vintage[i[1]]),
      size = if (has.size) df$Size[i[1]] else 1
    )
  })
  structure(list(units = units), class = "prepared.units")
}

# Horizon error of one prepared unit; same arithmetic as f1() in
# estim_model_optimized.R on the data-frame rows of that unit
unit.horizon.error <- function(u, max.month, par0, sdf.model = "linear") {
  X <- u$X[, names(par0), drop = FALSE]
  if (sdf.model == "exp.aff") {
    if ("Alpha" %in% colnames(X)) X[, "Alpha"] <- exp(1) - 1
    ret <- exp(cumsum(log(1 + X) %*% par0))
  } else {
    ret <- exp(cumsum(log(1 + (X %*% par0))))
  }
  switch(sdf.model,
    linear = , exp.aff = horizon.error(u$cf, ret, max.month),
    linear.duration = getNPVsDuration(u$cf, ret, max.month),
    linear.single.date = getNPVsSingleDate(u$cf, ret, max.month),
    stop("unknown sdf.model: ", sdf.model)
  )
}

# Funds of a prepared vintage (U_pair SE units): one prepared object per vintage
prepare.vintage.units <- function(df, factor.cols) {
  lapply(split(df, paste(df$type, df$Vintage), drop = TRUE), function(d) prepare.units(d, factor.cols))
}

# Horizon pricing error of one unit (B2): alias max.month 0 -> 1
horizon.error <- function(cf, ret, max.month) {
  getNPVs(cf, ret, max(as.integer(max.month), 1L))
}

# Cross-fund pair-product criterion U_{v,H} of one vintage (B5)
#   e: horizon errors of the vintage's funds, in the fund-level file's units
#      (EW: per unit of commitment; FW: per unit of commitment x fund size)
#   s: fund weights (EW: 1; FW: fund size), so that e = s * e.unit for FW.
#   U = [ (sum e)^2 - sum e^2 ] / [ (sum s)^2 - sum s^2 ]
#     = sum_{i != k} s_i s_k e.unit_i e.unit_k / sum_{i != k} s_i s_k,
#   unbiased for mu^2 + c under conditional independence (EW: ebar^2 - s^2/n).
u.pair <- function(e, s = rep(1, length(e))) {
  if (length(e) < 2) return(NA_real_)
  (sum(e)^2 - sum(e^2)) / (sum(s)^2 - sum(s^2))
}

# Normalized vintage weights w_v = V q_v / sum(q) for the vintages used by the
# U criterion (n_v >= 2), with q_v = 1 ("equal") or q_v = n_v ("n.funds").
# This is the single definition of the unit contributions shared by
# estimation and inference: Q(theta) = V^{-1} sum_v w_v U_v(theta), so the
# sandwich in helper/inference.R uses l_v = w_v U_v (review fix, 2026-11-16).
u.unit.weights <- function(n.funds, vintage.weights = "equal") {
  stopifnot(vintage.weights %in% c("equal", "n.funds"))
  ok <- n.funds >= 2
  q <- if (vintage.weights == "n.funds") n.funds else rep(1, length(n.funds))
  w <- rep(NA_real_, length(n.funds))
  w[ok] <- sum(ok) * q[ok] / sum(q[ok])
  names(w) <- names(n.funds)
  w
}

# Aggregate unit errors into the U criterion: weighted mean of U_v
#   vintage.weights = "equal" (q_v = 1) or "n.funds" (q_v = n_v);
#   vintages with fewer than two funds are dropped.
u.criterion <- function(e, vintage.key, s = rep(1, length(e)), vintage.weights = "equal") {
  u <- tapply(seq_along(e), vintage.key, function(i) u.pair(e[i], s[i]))
  n <- tapply(e, vintage.key, length)
  w <- u.unit.weights(n, vintage.weights)
  ok <- !is.na(w)
  if (!any(ok)) return(NA_real_)
  mean(w[ok] * u[ok])
}

# Typical parameter scales for the optimizer (optim/optimx parscale): alpha is
# a monthly return (~0.001-0.01), factor loadings are ~1. With equal scales
# L-BFGS-B effectively moved only alpha and stopped near the starting beta in
# the (alpha, beta) model (G1 finding, 2026-11-16; test in test_code_repairs.R).
optim.parscale <- function(par) ifelse(names(par) == "Alpha", 0.001, 1)

# ---------------------------------------------------------------------------
# Fast criterion (Step 4 Monte Carlo; decision of 18 November 2026: used in
# the Monte Carlo only, the empirical pipeline keeps unit.errors()).
#
# Same arithmetic as unit.horizon.error() for the linear SDF: all units that
# share a grid (the funds of a vintage start at t0 and end at T.max) share the
# benchmark path ret_t = prod G_u, so ret is computed once per grid instead of
# once per unit, and each unit only needs its nonzero cash flows:
#   e_i = (sum_t cf_it / ret_t) * mean(ret_1, ..., ret_min(m, n)).
# Equivalence with unit.errors() / u.criterion(): helper/test_mc.R (T1).
# ---------------------------------------------------------------------------

# df: padded unit rows (Date, CF, Fund.ID, type, Vintage, optional Size, and
# the factor columns), as passed to prepare.units(). Units in the same order
# as prepare.units() (split by Fund.ID).
prepare.fast <- function(df, factor.cols, split.key = df$Fund.ID) {
  factor.cols <- intersect(factor.cols, colnames(df))
  date.all <- as.numeric(df$Date)
  # common calendar: one row of factor values per date
  first <- !duplicated(date.all)
  cal <- data.frame(date = date.all[first], df[first, factor.cols, drop = FALSE])
  cal <- cal[order(cal$date), , drop = FALSE]
  X <- as.matrix(cal[, factor.cols, drop = FALSE])
  row <- match(date.all, cal$date)
  rows <- split(seq_len(nrow(df)), split.key, drop = TRUE)
  k0 <- vapply(rows, function(i) min(row[i]), 0L)
  n <- vapply(rows, length, 0L)
  # every unit must occupy the contiguous calendar rows k0, ..., k0 + n - 1
  ok <- vapply(seq_along(rows), function(j) identical(sort(row[rows[[j]]]), as.integer(k0[j] + seq_len(n[j]) - 1L)), TRUE)
  if (!all(ok)) stop("prepare.fast: units are not on contiguous calendar grids (pad them first)")
  # grids shared by units: (start row, length)
  grid.key <- paste(k0, n)
  grids <- unique(data.frame(k0 = k0, n = n, key = grid.key, stringsAsFactors = FALSE))
  unit.grid <- match(grid.key, grids$key)
  # nonzero cash flows: unit, position on the unit's grid
  cf.unit <- integer(); cf.pos <- integer(); cf.val <- numeric()
  for (j in seq_along(rows)) {
    i <- rows[[j]]; nz <- i[df$CF[i] != 0]
    cf.unit <- c(cf.unit, rep(j, length(nz))); cf.pos <- c(cf.pos, row[nz] - k0[j] + 1L); cf.val <- c(cf.val, df$CF[nz])
  }
  key <- vapply(rows, function(i) paste(df$type[i[1]], df$Vintage[i[1]]), "")
  groups <- sort(unique(key)) # vintage keys in u.criterion() order (tapply)
  group <- match(key, groups)
  structure(list(X = X, grids = grids[, c("k0", "n")], unit.grid = unit.grid,
    cf.unit = cf.unit, cf.pos = cf.pos, cf.row = k0[cf.unit] + cf.pos - 1L, cf.val = cf.val, n.units = length(rows),
    fund.id = names(rows), key = key, group = group, n.groups = length(groups), group.n = setNames(tabulate(group, length(groups)), groups),
    vintage = vapply(rows, function(i) as.numeric(df$Vintage[i[1]]), 0),
    size = if ("Size" %in% colnames(df)) vapply(rows, function(i) df$Size[i[1]], 0) else rep(1, length(rows))),
    class = "fast.units")
}

# Horizon errors of all units at par (named; RF enters with loading 1, as in
# unit.errors()). Non-finite where a gross return is not positive (log NaN),
# like the pipeline. The cumulative log return is taken once over the
# calendar; a grid starting at row k0 uses its differences to row k0 - 1
# (equal to the per-unit cumulative sum up to rounding, T1).
fast.unit.errors <- function(par, fp, max.month) {
  par0 <- c(RF = 1, par)
  g <- fp$grids
  # compound only from the earliest grid start: a non-positive gross return in
  # an earlier calendar month (outside every unit's grid) must not turn the
  # errors into NaN, as it does not in the per-unit pipeline computation
  r0 <- min(g$k0)
  L <- rep(NA_real_, nrow(fp$X) + 1L)
  L[r0:(nrow(fp$X) + 1L)] <- c(0, cumsum(log(1 + as.numeric(fp$X[r0:nrow(fp$X), names(par0), drop = FALSE] %*% par0))))
  m <- max(as.integer(max.month), 1L)
  base <- L[g$k0] # cumulative log return before the grid's first row
  A <- vapply(seq_len(nrow(g)), function(k) {
    j <- min(m, g$n[k])
    sum(exp(L[g$k0[k] + seq_len(j)] - base[k])) / j
  }, 0)
  disc <- exp(L[fp$cf.row + 1L] - base[fp$unit.grid[fp$cf.unit]])
  dcf <- numeric(fp$n.units)
  s <- rowsum(fp$cf.val / disc, fp$cf.unit, reorder = TRUE)
  dcf[as.integer(rownames(s))] <- s[, 1]
  setNames(dcf * A[fp$unit.grid], fp$fund.id)
}

# Horizon scale of every unit at par: the factor z with e_H = z * e_0 (the
# rescaling identity T1, relative to the inception NPV e_0 = the error at
# max.month 0), i.e. A_{v,H} / A_{v,0} in fast.unit.errors(); 1 at H = 0.
# The instrument of the IV comparator (simulation/mc_driver.R). Same
# compounding as fast.unit.errors() (test: e_H = z * e_0, helper/test_mc.R).
fast.horizon.scale <- function(par, fp, max.month) {
  par0 <- c(RF = 1, par)
  g <- fp$grids
  r0 <- min(g$k0)
  L <- rep(NA_real_, nrow(fp$X) + 1L)
  L[r0:(nrow(fp$X) + 1L)] <- c(0, cumsum(log(1 + as.numeric(fp$X[r0:nrow(fp$X), names(par0), drop = FALSE] %*% par0))))
  m <- max(as.integer(max.month), 1L)
  base <- L[g$k0]
  A <- vapply(seq_len(nrow(g)), function(k) sum(exp(L[g$k0[k] + seq_len(min(m, g$n[k]))] - base[k])) / min(m, g$n[k]), 0)
  A1 <- exp(L[g$k0 + 1L] - base) # the max.month 0 (= 1) scale
  setNames((A / A1)[fp$unit.grid], fp$fund.id)
}

# U criterion of fast.unit.errors() output: same value as u.criterion(e,
# fp$key, s, vintage.weights) (T1), with the vintage grouping precomputed.
# Fund weights s: equal (EW, the default, as the pipeline for EW data even
# when the file has a Size column); pass s = fp$size for FW (fund sizes).
fast.u.criterion <- function(e, fp, s = rep(1, length(e)), vintage.weights = "equal") {
  sum.by <- function(x) { out <- numeric(fp$n.groups); r <- rowsum(x, fp$group, reorder = TRUE); out[as.integer(rownames(r))] <- r[, 1]; out }
  u <- (sum.by(e)^2 - sum.by(e^2)) / (sum.by(s)^2 - sum.by(s^2))
  w <- u.unit.weights(fp$group.n, vintage.weights)
  ok <- !is.na(w)
  if (!any(ok)) return(NA_real_)
  mean(w[ok] * u[ok])
}

# Units of selected vintage keys (CV folds, jackknife, per-vintage SE units)
fast.subset <- function(fp, keep) {
  j <- which(keep)
  map <- match(seq_len(fp$n.units), j)
  sel <- fp$cf.unit %in% j
  used <- sort(unique(fp$unit.grid[j]))
  out <- fp
  out$grids <- fp$grids[used, , drop = FALSE]
  out$unit.grid <- match(fp$unit.grid[j], used)
  out$cf.unit <- map[fp$cf.unit[sel]]; out$cf.pos <- fp$cf.pos[sel]; out$cf.row <- fp$cf.row[sel]; out$cf.val <- fp$cf.val[sel]
  out$n.units <- length(j)
  for (f in c("fund.id", "key", "vintage", "size")) out[[f]] <- fp[[f]][j]
  groups <- sort(unique(out$key)); out$group <- match(out$key, groups); out$n.groups <- length(groups)
  out$group.n <- setNames(tabulate(out$group, length(groups)), groups)
  out
}

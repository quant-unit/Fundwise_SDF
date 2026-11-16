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

pad.all.units <- function(df, origin.convention = "first.cf", t.max = max(df$Date)) {
  t0 <- vintage.origins(df, origin.convention)
  units <- split(df, df$Fund.ID, drop = TRUE)
  out <- lapply(units, function(d) to.monthly.grid(d, t0[[paste(d$type[1], d$Vintage[1])]], t.max))
  as.data.frame(data.table::rbindlist(out, use.names = TRUE, fill = TRUE))
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

# Aggregate unit errors into the U criterion: vintage-weighted mean of U_v
#   vintage.weights = "equal" (q_v = 1) or "n.funds" (q_v = n_v);
#   vintages with fewer than two funds are dropped.
u.criterion <- function(e, vintage.key, s = rep(1, length(e)), vintage.weights = "equal") {
  u <- tapply(seq_along(e), vintage.key, function(i) u.pair(e[i], s[i]))
  n <- tapply(e, vintage.key, length)
  ok <- !is.na(u)
  if (!any(ok)) return(NA_real_)
  q <- if (vintage.weights == "n.funds") n[ok] else rep(1, sum(ok))
  sum(q * u[ok]) / sum(q)
}

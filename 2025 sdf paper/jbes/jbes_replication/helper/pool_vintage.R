######################
### Constant-composition vintage-year portfolio (VYP) pooling
### Revision plan B4 (constant composition) and B6 (FW normalization)
######################
#
# Input: fund-level cash flows of ONE vintage (and one strategy) with columns
#   Fund.ID, Date, CF, NAV, CON, DIS, Vintage, type, Size
# where CF/CON/DIS are flows and NAV is a stock, all per unit of commitment
# (EW) or per unit of commitment times fund Size in USD mn (FW).
#
# Pooled flow on date t:
#   EW:  CF_{v,t} = sum_i CF_{v,i,t} / N_v
#   FW:  CF_{v,t} = sum_i Size_i * cf_{v,i,t} / sum_i Size_i
# where the denominators are time-invariant (all funds of the vintage),
# and a fund without a transaction on date t contributes CF_{v,i,t} = 0.
# Because the denominator does not depend on t, the pooled NPV is the exact
# weighted average of fund NPVs for every discount curve:
#   PV_v(theta) = sum_i w_i PV_{v,i}(theta),  w_i = 1/N_v (EW) or Size_i/sum Size (FW).
#
# This replaces aggregate(CF ~ Date, df, mean), which divided by the number
# of funds reporting on date t (reporter-mean pooling).
#
# NAV is a stock: it is carried forward between a fund's first and last
# reporting date and is zero after the final date (the final NAV is booked
# as a liquidating distribution in CF on that date).

pool.vintage <- function(df, fund.size.weighting) {
  stopifnot(
    length(unique(df$Vintage)) == 1,
    length(unique(df$type)) == 1,
    all(c("Fund.ID", "Date", "CF", "NAV", "CON", "DIS", "Size") %in% colnames(df))
  )
  df$Date <- as.Date(df$Date)
  fund.ids <- unique(as.character(df$Fund.ID))
  sizes <- df$Size[!duplicated(as.character(df$Fund.ID))]
  n.funds <- length(fund.ids)
  denom <- if (fund.size.weighting) sum(sizes) else n.funds

  dates <- sort(unique(df$Date))

  # flows: sum over all funds on each date
  sum.by.date <- function(x) {
    s <- tapply(x, factor(as.character(df$Date), levels = as.character(dates)), sum)
    s[is.na(s)] <- 0
    as.numeric(s)
  }
  cf <- sum.by.date(df$CF)
  con <- sum.by.date(df$CON)
  dis <- sum.by.date(df$DIS)

  # stock: carry each fund's NAV forward within its reporting life
  nav <- numeric(length(dates))
  for (id in fund.ids) {
    d.i <- df[as.character(df$Fund.ID) == id, ]
    d.i <- d.i[order(d.i$Date), ]
    nav.i.last <- tapply(d.i$NAV, as.character(d.i$Date), function(x) x[length(x)])
    idx <- findInterval(as.numeric(dates), as.numeric(as.Date(names(nav.i.last))))
    alive <- idx > 0 & dates <= max(d.i$Date)
    nav[alive] <- nav[alive] + as.numeric(nav.i.last)[idx[alive]]
  }

  vintage <- df$Vintage[1]
  type <- df$type[1]
  data.frame(
    Date = dates,
    CF = cf / denom,
    NAV = nav / denom,
    CON = con / denom,
    DIS = dis / denom,
    type = type,
    Vintage = vintage,
    Fund.ID = paste0(vintage, "_", type),
    N.funds = n.funds,
    Commitment = sum(sizes),
    stringsAsFactors = FALSE
  )
}

# Pool every vintage of a fund-level data frame (single strategy `type`).
pool.all.vintages <- function(df, fund.size.weighting) {
  if (nrow(df) == 0) return(df)
  dfx <- split(df, df$Vintage, drop = TRUE)
  out <- do.call(rbind, lapply(dfx, pool.vintage, fund.size.weighting = fund.size.weighting))
  rownames(out) <- NULL
  out
}

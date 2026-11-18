######################
### Step 4 Monte Carlo: calibration of the deal model and sigma (decision (ii))
### Revision plan, Section I, "Step 4: scope and compute plan"
######################
#
# Moments per type (vintages up to 2010, EW funds, per unit of commitment):
#   duration   : median years from a fund's first to its last cash flow
#                (truncated at the data end, as in the data and the generator)
#   dispersion : pooled within-vintage sd of the fund-level H = 0 pricing
#                errors e_i, at the empirical single-factor estimate
#                (L2 on vintage portfolios, EW, H = 0) in the data and at
#                theta_0 in the simulation; the shared vintage component is
#                removed, so it measures fund-level noise, which sigma sets.
# Simulated funds commit no.deals units (one per deal): their errors are
# divided by no.deals.
source("simulation/mc_driver.R")

pooled.within.sd <- function(e, vintage) {
  r <- e - ave(e, vintage)
  sqrt(sum(r^2) / (length(e) - length(unique(vintage))))
}
fund.duration <- function(date.first, date.last) as.numeric(date.last - date.first) / 365.25

# Empirical moments
mc.empirical.moments <- function(type, public, beta.hat, max.vintage = 2010,
                                 file = "empirical/data_prepared_2026_c/preqin_cashflows_EW_North America_NAV.csv") {
  d <- read.csv(file); d$Date <- as.Date(d$Date)
  d <- d[d$type == type & d$Vintage <= max.vintage, ]
  nz <- d[d$CF != 0, ]
  dur <- tapply(nz$Date, nz$Fund.ID, function(x) fund.duration(min(x), max(x)))
  d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  all <- read.csv(file); all$Date <- as.Date(all$Date)
  p <- pad.all.units(d, "first.cf", max(all$Date)) # common data end of the file
  p <- merge(p, public, by = "Date", all.x = TRUE)
  fp <- prepare.fast(p, unique(c("RF", setdiff(colnames(public), "Date"))))
  e <- fast.unit.errors(c(MKT = beta.hat), fp, 1)
  c(duration = median(dur), duration.q1 = unname(quantile(dur, 0.25)), duration.q3 = unname(quantile(dur, 0.75)),
    dispersion = pooled.within.sd(e, fp$vintage), funds = length(e))
}

# Simulated moments of a design (average over replications seeds)
mc.sim.moments <- function(design, public, seeds = 1:5) {
  m <- vapply(seeds, function(s) {
    dr <- mc.draw(design, s, public)
    d <- dr$deals
    f <- paste(d$vintage, d$fund)
    first <- tapply(d$start, f, min); last <- tapply(d$end, f, max)
    dur <- fund.duration(dr$dates[first], dr$dates[last])
    fp <- mc.as.fast(dr, public, "fund")
    e <- fast.unit.errors(design$theta0, fp, 1) / design$no.deals
    c(duration = median(dur), dispersion = pooled.within.sd(e, fp$vintage))
  }, c(duration = 0, dispersion = 0))
  rowMeans(m)
}

# Acceptance tests of the Step 4 Monte Carlo (revision plan, Section I,
# "Step 4: scope and compute plan", implementation spec).
# Run from r_project/: LANG=en_US.UTF-8 Rscript helper/test_mc.R
#   T1  fast criterion (helper/objective.R) = unit.errors() / u.criterion()
#       and the vintage-portfolio criterion, within 1e-10 (relative)
#   T2-T4 follow with the generator and the driver.
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R"))))) # horizon.error (Rcpp)
source(here("helper/objective.R"))
n.fail <- 0L
check <- function(label, ok, detail = "") {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", label, if (nzchar(detail)) paste0("(", detail, ")"), "\n")
  if (!isTRUE(ok)) n.fail <<- n.fail + 1L
}

# Public factors and padded private data exactly as estim_model_optimized.R
# builds them (date normalization, Alpha column, grid to the common data end).
load.public <- function(file = here("empirical/data_prepared_2026/q_factors.csv")) {
  p <- read.csv(file); p$Date <- as.Date(p$Date)
  p$Date <- as.Date(format(p$Date, "%Y-%m-01"))
  p$Date <- seq.Date(p$Date[1], by = "month", length.out = nrow(p))
  p$Date <- as.Date(format(p$Date + 32, "%Y-%m-01")) - 1
  p$Alpha <- 1
  p
}
load.units <- function(private.file, public, types, max.vintage = 2010) {
  d <- read.csv(private.file)
  d$Date <- as.Date(d$Date)
  d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  d <- pad.all.units(d, "first.cf", max(d$Date))
  d$Fund.ID <- as.factor(as.character(d$Fund.ID))
  d <- merge(d, public, by = "Date", all.x = TRUE)
  d <- d[d$type %in% types & d$Vintage <= max.vintage, ]
  d$Fund.ID <- droplevels(as.factor(d$Fund.ID))
  d
}
factor.cols <- function(public) unique(c("RF", setdiff(colnames(public), "Date")))

# Reference: the pipeline's unit errors (unit.errors() in estim_model_optimized.R)
ref.errors <- function(par, prep, m) vapply(prep$units, unit.horizon.error, 0, max.month = m, par0 = c(RF = 1, par), sdf.model = "linear")
rel.diff <- function(a, b) max(abs(a - b)) / max(abs(b), .Machine$double.xmin)

# Parameter grid: single factor, (alpha, beta), MKT + EG, including points
# near the domain edge; horizons in months (0 is the alias of 1)
pars <- list(c(MKT = 1), c(MKT = 0.3), c(MKT = 2.6), c(MKT = 1.2, Alpha = 0.005), c(MKT = 0.8, Alpha = -0.01),
  c(MKT = 1, EG = 0.5), c(MKT = 2.6214, EG = 12.7806))
months <- c(0, 1, 30, 120, 360)

t1.compare <- function(label, df, public, criteria = c("L2", "U")) {
  fc <- factor.cols(public)
  prep <- prepare.units(df, fc)
  t <- system.time(fp <- prepare.fast(df, fc))[["elapsed"]]
  check(paste0("T1 ", label, ": same units in the same order"), identical(names(prep$units), fp$fund.id) &&
    identical(unname(vapply(prep$units, `[[`, "", "key")), unname(fp$key)), sprintf("%d units, prepared in %.1f s", fp$n.units, t))
  worst <- c(errors = 0, L2 = 0, U = 0, U.n = 0, U.size = 0, U.fast = 0)
  s <- vapply(prep$units, `[[`, 0, "size")
  key <- vapply(prep$units, `[[`, "", "key")
  for (par in pars) for (m in months) {
    a <- fast.unit.errors(par, fp, m); b <- ref.errors(par, prep, m)
    if (!all(is.finite(b))) { worst["errors"] <- max(worst["errors"], if (identical(is.finite(a), is.finite(b))) 0 else Inf); next }
    worst["errors"] <- max(worst["errors"], rel.diff(a, b))
    if ("L2" %in% criteria) worst["L2"] <- max(worst["L2"], rel.diff(sqrt(sum(a^2)) / length(a), sqrt(sum(b^2)) / length(b)))
    if ("U" %in% criteria) {
      worst["U"] <- max(worst["U"], rel.diff(u.criterion(a, fp$key, rep(1, length(a))), u.criterion(b, key, rep(1, length(b)))))
      worst["U.n"] <- max(worst["U.n"], rel.diff(u.criterion(a, fp$key, rep(1, length(a)), "n.funds"), u.criterion(b, key, rep(1, length(b)), "n.funds")))
      worst["U.size"] <- max(worst["U.size"], rel.diff(u.criterion(a, fp$key, fp$size), u.criterion(b, key, s)))
      for (vw in c("equal", "n.funds")) for (sz in list(rep(1, length(a)), fp$size)) {
        worst["U.fast"] <- max(worst["U.fast"], rel.diff(fast.u.criterion(a, fp, sz, vw), u.criterion(b, key, if (length(unique(sz)) == 1) rep(1, length(b)) else s, vw)))
      }
    }
  }
  for (k in c("errors", criteria, if ("U" %in% criteria) c("U.n", "U.size", "U.fast"))) {
    check(sprintf("T1 %s: %s within 1e-10", label, k), worst[[k]] <= 1e-10, sprintf("max relative difference %.1e", worst[[k]]))
  }
  # subsets of vintages (CV folds, SE units) give the same errors as preparing the subset
  keep <- fp$vintage %in% sort(unique(fp$vintage))[c(TRUE, FALSE)]
  sub.ref <- ref.errors(c(MKT = 1.1), prepare.units(df[df$Vintage %in% unique(fp$vintage[keep]), ], fc), 120)
  check(sprintf("T1 %s: vintage subset equals the subset prepared on its own", label),
    rel.diff(fast.unit.errors(c(MKT = 1.1), fast.subset(fp, keep), 120), sub.ref) <= 1e-10)
  invisible(fp)
}

public <- load.public()
ew <- here("empirical/data_prepared_2026_c/preqin_cashflows_EW_North America_NAV.csv")
fw <- here("empirical/data_prepared_2026_c/preqin_cashflows_FW_North America_NAV.csv")
vyp <- here("empirical/data_prepared_2026_c/preqin_cashflows_EW_VYP_North America_NAV.csv")
for (t in c("BO", "VC")) {
  t1.compare(paste(t, "EW funds"), load.units(ew, public, t), public)
  t1.compare(paste(t, "EW vintage portfolios"), load.units(vyp, public, t), public, criteria = "L2")
}
t1.compare("BO FW funds (fund-size weights)", load.units(fw, public, "BO"), public)

# Speed on real data (Buyout EW funds, U criterion)
d.bo <- load.units(ew, public, "BO"); prep <- prepare.units(d.bo, factor.cols(public)); fp <- prepare.fast(d.bo, factor.cols(public))
t.ref <- system.time(for (k in 1:20) u.criterion(ref.errors(c(MKT = 1 + k / 1e4), prep, 120), fp$key))[["elapsed"]] / 20
t.fast <- system.time(for (k in 1:200) fast.u.criterion(fast.unit.errors(c(MKT = 1 + k / 1e4), fp, 120), fp))[["elapsed"]] / 200
cat(sprintf("Speed, Buyout EW funds (%d units), U criterion at 10 years: pipeline %.2f ms, fast %.2f ms per evaluation (x%.0f)\n",
  fp$n.units, 1000 * t.ref, 1000 * t.fast, t.ref / t.fast))

cat("\n", n.fail, "test(s) failed\n")
if (n.fail > 0) quit(status = 1)

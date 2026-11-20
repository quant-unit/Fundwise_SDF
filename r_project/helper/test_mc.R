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

# ---------------------------------------------------------------------------
# T2: in-memory generator (simulation/mc_generate.R); T1 on a simulated
# replication (fast structure built directly = prepare.fast() on its rows =
# the pipeline's unit errors)
# ---------------------------------------------------------------------------
source(here("simulation/mc_generate.R"))
pub <- mc.public.factors(here("empirical/data_prepared_2026/q_factors.csv"))
shape <- mc.data.shape("BO", file = ew)
# padded rows of a replication, merged with the factors as the pipeline does
sim.rows <- function(dr, level) {
  d <- mc.as.rows(dr, level); d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  d <- pad.all.units(d, "first.cf", max(d$Date)); d$Fund.ID <- as.factor(as.character(d$Fund.ID))
  merge(d, pub, by = "Date", all.x = TRUE)
}
fc <- factor.cols(pub)
same.fast <- function(a, b) identical(a$fund.id, sub("_SIM$", "", b$fund.id)) || identical(sub("_SIM_SIM$", "_SIM", a$fund.id), b$fund.id)

d0 <- mc.design(shape, sigma = 0)
dr0 <- mc.draw(d0, 1, pub)
for (lev in c("fund", "vyp")) {
  f.direct <- mc.as.fast(dr0, pub, lev)
  e0 <- vapply(c(1, 120, 360), function(m) max(abs(fast.unit.errors(c(MKT = 1), f.direct, m))), 0)
  check(sprintf("T2 sigma = 0: every %s unit priced exactly at theta_0 (H = 0, 10, 30 years)", lev), max(e0) <= 1e-12, sprintf("max |e| %.1e", max(e0)))
}
dr <- mc.draw(mc.design(shape, sigma = 0.2), 1, pub)
check("T2 fund counts per vintage equal the design", identical(as.integer(table(unique(dr$deals[, c("vintage", "fund")])$vintage)), unname(shape)) &&
  identical(names(table(dr$deals$vintage)), names(shape)))
check("T2 every deal starts in or after January of its vintage and ends by T.max",
  all(dr$dates[dr$deals$start] >= as.Date(paste0(dr$deals$vintage, "-01-31"))) && all(dr$deals$end <= length(dr$dates)) && max(dr$dates) == as.Date("2021-12-31"))
check("T2 cash flows booked at quarter ends", all(as.integer(format(dr$dates[c(dr$deals$start, dr$deals$end)], "%m")) %% 3 == 0))
check("T2 same seed, same replication; other seed, other replication", identical(mc.draw(mc.design(shape), 7, pub), mc.draw(mc.design(shape), 7, pub)) &&
  !identical(mc.draw(mc.design(shape), 7, pub)$deals$mult, mc.draw(mc.design(shape), 8, pub)$deals$mult))
# mean-one shock: Z = m / R_{a,b}(theta_0), R recomputed here by direct products
G <- 1 + pub$RF[seq_along(dr$dates)] + pub$MKT[seq_along(dr$dates)]
# log Z = sum of k i.i.d. N(-s^2/2, s^2): standardized, it is N(0, 1) (the
# mean of Z itself is too heavy-tailed to test sharply)
for (sg in c(0.2, 0.3)) {
  x <- mc.draw(mc.design(shape, sigma = sg), 3, pub)$deals
  x <- x[x$end > x$start, ]
  k <- x$end - x$start; s2 <- log(1 + sg^2)
  z <- (log(x$mult / vapply(seq_len(nrow(x)), function(i) prod(G[(x$start[i] + 1):x$end[i]]), 0)) + k * s2 / 2) / sqrt(k * s2)
  check(sprintf("T2 sigma = %.1f: mean-one lognormal deal shock (standardized log shock: mean 0, sd 1, within 4 SEs)", sg),
    abs(mean(z)) <= 4 / sqrt(length(z)) && abs(sd(z) - 1) <= 4 / sqrt(2 * length(z)),
    sprintf("mean %.4f, sd %.4f, %d deals", mean(z), sd(z), length(z)))
}
w <- with(mc.draw(mc.design(shape, sigma = 0, shared.sd = 0.3), 5, pub)$deals, tapply(log(mult / vapply(seq_along(start), function(i) if (end[i] > start[i]) prod(G[(start[i] + 1):end[i]]) else 1, 0)), vintage, function(x) c(mean(x), sd(x))))
w <- do.call(rbind, w)
check("T2 shared vintage shock: constant within a vintage, sd across vintages near 0.3",
  max(w[, 2], na.rm = TRUE) < 1e-12 && abs(sd(w[, 1]) - 0.3) < 0.12, sprintf("sd across %d vintages %.3f", nrow(w), sd(w[, 1])))
# T1 on the simulated replication, fund level and vintage portfolios
for (lev in c("fund", "vyp")) {
  rows <- sim.rows(dr, lev)
  fp.rows <- prepare.fast(rows, fc); prep.rows <- prepare.units(rows, fc); fp.dir <- mc.as.fast(dr, pub, lev)
  worst <- c(direct = 0, pipeline = 0)
  for (par in pars) for (m in months) {
    a <- fast.unit.errors(par, fp.dir, m); b <- fast.unit.errors(par, fp.rows, m); r <- ref.errors(par, prep.rows, m)
    worst["direct"] <- max(worst["direct"], rel.diff(unname(a), unname(b)))
    worst["pipeline"] <- max(worst["pipeline"], rel.diff(unname(a), unname(r)))
  }
  check(sprintf("T1 simulated %s units: same units as the padded rows", lev), fp.dir$n.units == fp.rows$n.units && same.fast(fp.dir, fp.rows))
  check(sprintf("T1 simulated %s units: direct structure = prepare.fast() on the rows within 1e-10", lev), worst[["direct"]] <= 1e-10, sprintf("%.1e", worst[["direct"]]))
  check(sprintf("T1 simulated %s units: = pipeline unit errors within 1e-10", lev), worst[["pipeline"]] <= 1e-10, sprintf("%.1e", worst[["pipeline"]]))
  if (lev == "fund") check("T1 simulated funds: U criterion (fast) = u.criterion() on the pipeline errors within 1e-10",
    rel.diff(fast.u.criterion(fast.unit.errors(c(MKT = 1.2), fp.dir, 120), fp.dir), u.criterion(ref.errors(c(MKT = 1.2), prep.rows, 120), vapply(prep.rows$units, `[[`, "", "key"))) <= 1e-10)
}
# ---------------------------------------------------------------------------
# T4: checkpoints (simulation/mc_driver.R): a restart skips finished batches,
# a lost batch is recomputed identically, another specification is refused
# ---------------------------------------------------------------------------
source(here("simulation/mc_driver.R"))
spec <- list(design = mc.design(c("1990" = 3, "1991" = 4, "1993" = 2, "1995" = 5), sigma = 0.2),
             estimators = c("L2.vyp", "U.eq"), models = "MKT", horizons = c(0, 120))
out4 <- file.path(tempdir(), "mc_t4"); unlink(out4, recursive = TRUE)
first <- mc.run("t4", spec, 1:4, pub, out.dir = out4, batch.size = 2, workers = 2, verbose = FALSE)
check("T4 a run writes one checkpoint per batch", length(list.files(file.path(out4, "t4"))) == 2 && nrow(first) == 4 * 2 * 2)
mt <- file.info(file.path(out4, "t4", "batch_0001.rds"))$mtime
unlink(file.path(out4, "t4", "batch_0002.rds"))
again <- mc.run("t4", spec, 1:4, pub, out.dir = out4, batch.size = 2, workers = 2, verbose = FALSE)
drop.t <- function(d) { d$seconds <- NULL; attr(d, "provenance") <- NULL; d[order(d$rep, d$estimator, d$max.month), ] }
check("T4 a restart skips finished batches and recomputes the lost one identically",
  identical(file.info(file.path(out4, "t4", "batch_0001.rds"))$mtime, mt) && identical(drop.t(again), drop.t(first)))
spec2 <- spec; spec2$horizons <- c(0, 60)
check("T4 checkpoints of another specification are refused", inherits(tryCatch(mc.run("t4", spec2, 1:4, pub, out.dir = out4, batch.size = 2, workers = 1, verbose = FALSE), error = function(e) e), "error"))
unlink(out4, recursive = TRUE)

# extension: more replications of the same specification keep the finished batches
spec.r <- spec; spec.r$reps <- 2
unlink(out4, recursive = TRUE)
mc.run("t4", spec.r, 1:2, pub, out.dir = out4, batch.size = 2, workers = 1, verbose = FALSE)
mt <- file.info(file.path(out4, "t4", "batch_0001.rds"))$mtime
spec.r$reps <- 4
ext <- mc.run("t4", spec.r, 1:4, pub, out.dir = out4, batch.size = 2, workers = 2, verbose = FALSE)
check("T4 extending a design to more replications keeps its finished batches",
  identical(file.info(file.path(out4, "t4", "batch_0001.rds"))$mtime, mt) && identical(drop.t(ext), drop.t(first)))
# failures: a failing replication leaves no finished batch, only a _failed file
spec.f <- spec; spec.f$models <- "NOT_A_FACTOR"
res.f <- suppressWarnings(mc.run("t4f", spec.f, 1:2, pub, out.dir = out4, batch.size = 2, workers = 2, verbose = FALSE))
check("T4 a batch with failed replications is not marked finished (kept as _failed for inspection)",
  !file.exists(file.path(out4, "t4f", "batch_0001.rds")) && file.exists(file.path(out4, "t4f", "batch_0001_failed.rds")) && nrow(res.f) == 0)
unlink(out4, recursive = TRUE)

# fast criterion: a gross return <= 0 only before every unit's grid must not
# turn the errors into NaN (the pipeline compounds over each unit's grid only)
pre <- seq_len(min(mc.as.fast(dr, pub, "fund")$grids$k0) - 1L)
win <- setdiff(seq_along(dr$dates), pre)
# (no actual parameter point does this with the 1967-2021 factors, so the test
# makes one calendar month before every grid explosive: RF = -5 gives G < 0
# there; that month is unused, so the errors must equal those without it)
pub.bad <- pub; pub.bad$RF[pre[length(pre)]] <- -5
for (par.t in list(c(MKT = 1), c(MKT = 1.2, Alpha = 0.002))) for (m in c(1, 120)) {
  a <- fast.unit.errors(par.t, mc.as.fast(dr, pub.bad, "fund"), m); b <- fast.unit.errors(par.t, mc.as.fast(dr, pub, "fund"), m)
  if (!(all(is.finite(a)) && identical(unname(a), unname(b)))) { check("T1 errors unaffected by a non-positive gross return before every grid", FALSE); break }
}
check("T1 errors unaffected by a non-positive gross return before every grid (finite, identical)", all(is.finite(a)) && identical(unname(a), unname(b)))

# ---------------------------------------------------------------------------
# T5-T7: the Step 4 additions of 19 November (plan, Section J, "Gap found"):
# T5 persistent shared vintage shocks, T6 resampled factor paths, T7 the IV
# comparator. Designs without the new options are unchanged (spec hashes and
# draws of all 28 production designs were compared before and after the
# change; scratch check, plan log).
# ---------------------------------------------------------------------------
cons <- function(funds, count, last = 2005) setNames(rep(as.integer(funds), count), as.character(seq(last - count + 1, last)))
# the shared log shock of each vintage, recovered from sigma = 0 deals: log(mult / R_{a,b}(theta_0))
shared.log <- function(dr, G) with(dr$deals, tapply(log(mult / vapply(seq_along(start), function(i) if (end[i] > start[i]) prod(G[(start[i] + 1):end[i]]) else 1, 0)), vintage, mean))
check("T5 shared.rho = 0 leaves the design (and so its hash and draws) as without the option",
  identical(mc.design(shape, shared.sd = 0.3), mc.design(shape, shared.sd = 0.3, shared.rho = 0)))
dp <- mc.design(cons(2, 35), sigma = 0, shared.sd = 0.3, shared.rho = 0.6)
W <- sapply(1:150, function(s) shared.log(mc.draw(dp, s, pub), G))
# (tolerance of the mean: 4 Monte Carlo SEs; a replication's mean over 35 AR(1) vintages has sd about 0.3 * sqrt((1 + 0.6) / (1 - 0.6) / 35) = 0.10, so 0.0083 over 150 replications)
within <- with(mc.draw(dp, 1, pub)$deals, tapply(log(mult / vapply(seq_along(start), function(i) if (end[i] > start[i]) prod(G[(start[i] + 1):end[i]]) else 1, 0)), vintage, sd))
lag1 <- cor(as.vector(W[-1, ]), as.vector(W[-nrow(W), ]))
check("T5 persistent shock: constant within a vintage; marginal mean -c^2/2 and sd c; lag-1 correlation across vintages near rho",
  max(within, na.rm = TRUE) < 1e-12 && abs(mean(W) + 0.045) < 0.035 && abs(sd(as.vector(W)) - 0.3) < 0.02 && abs(lag1 - 0.6) < 0.04,
  sprintf("mean %.4f (-0.045), sd %.4f, lag-1 correlation %.3f", mean(W), sd(as.vector(W)), lag1))
# gaps: vintages two calendar years apart are correlated rho^2
dg <- mc.design(setNames(rep(2L, 20), as.character(seq(1968, 2006, by = 2))), sigma = 0, shared.sd = 0.3, shared.rho = 0.6)
Wg <- sapply(1:150, function(s) shared.log(mc.draw(dg, s, pub), G))
lag.g <- cor(as.vector(Wg[-1, ]), as.vector(Wg[-nrow(Wg), ]))
check("T5 persistent shock: vintages two years apart correlated near rho^2 (AR(1) in calendar years)", abs(lag.g - 0.36) < 0.05, sprintf("%.3f", lag.g))

path <- list(resample = "iid", from = "1967-01-31", to = "2021-12-31")
dpath <- mc.design(cons(2, 100), sigma = 0, path = path)
drp <- mc.draw(dpath, 4, pub)
src <- pub[pub$Date >= as.Date(path$from) & pub$Date <= as.Date(path$to), ]
fcols <- setdiff(colnames(pub), "Date")
check("T6 resampled path: monthly month ends from January of the first vintage to T.max, as the deal calendar",
  identical(drp$public$Date, drp$dates) && drp$dates[1] == as.Date("1906-01-31") && max(drp$dates) == as.Date("2021-12-31") &&
  all(as.integer(diff(as.integer(format(drp$dates, "%Y")) * 12 + as.integer(format(drp$dates, "%m")))) == 1L))
key <- function(d) do.call(paste, c(lapply(d[, fcols], function(x) sprintf("%.12g", x)), sep = "|"))
check("T6 resampled path: every month is a whole historical month (all factors and r_f together)", all(key(drp$public) %in% key(src)))
check("T6 resampled path: same seed, same path; other seed, other path",
  identical(mc.draw(dpath, 4, pub), drp) && !identical(mc.draw(dpath, 5, pub)$public$MKT, drp$public$MKT))
for (lev in c("fund", "vyp")) {
  f.p <- mc.as.fast(drp, pub, lev)
  e0 <- vapply(c(1, 120, 360), function(m) max(abs(fast.unit.errors(c(MKT = 1), f.p, m))), 0)
  check(sprintf("T6 sigma = 0 on a resampled path: every %s unit priced exactly at theta_0 (H = 0, 10, 30 years)", lev), max(e0) <= 1e-12, sprintf("max |e| %.1e", max(e0)))
}
# T1 on a resampled-path replication: the direct structure equals the pipeline's unit errors on its rows
drp2 <- mc.draw(mc.design(cons(3, 30), sigma = 0.2, path = path), 6, pub)
rows.p <- local({ d <- mc.as.rows(drp2, "fund"); d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  d <- pad.all.units(d, "first.cf", max(d$Date)); d$Fund.ID <- as.factor(as.character(d$Fund.ID)); merge(d, drp2$public, by = "Date", all.x = TRUE) })
prep.p <- prepare.units(rows.p, fc); fp.p <- mc.as.fast(drp2, pub, "fund")
worst.p <- max(unlist(lapply(pars[1:3], function(par) lapply(months, function(m) rel.diff(unname(fast.unit.errors(par, fp.p, m)), unname(ref.errors(par, prep.p, m)))))))
check("T1 resampled path (funds): direct structure = pipeline unit errors within 1e-10", worst.p <= 1e-10, sprintf("%.1e", worst.p))
rep.p <- mc.replication(mc.design(cons(5, 60), sigma = 0.13, path = path), 2, pub, c("L2.vyp", "U.eq", "IV"), "MKT", c(0, 180))
check("T6 the driver runs a resampled-path design: finite estimates and SEs for L2, U and IV",
  nrow(rep.p) == 6 && all(is.finite(rep.p$MKT)) && all(is.finite(rep.p$SE.MKT[!rep.p$at.bound])))

# T7: IV comparator. (a) the instrument is the T1 scale: e_H = z * e_0
fpv <- mc.as.fast(dr, pub, "vyp")
worst.z <- max(unlist(lapply(pars[1:3], function(par) lapply(c(60, 120, 360), function(m)
  rel.diff(unname(fast.unit.errors(par, fpv, m)), unname(fast.horizon.scale(par, fpv, m) * fast.unit.errors(par, fpv, 0)))))))
check("T7 instrument: e_H = z * e_0 within 1e-12 (rescaling identity), z = 1 at H = 0",
  worst.z <= 1e-12 && max(abs(fast.horizon.scale(c(MKT = 1.3), fpv, 0) - 1)) < 1e-14, sprintf("%.1e", worst.z))
# (b) exact pricing: the root is theta_0
r0 <- do.call(rbind, lapply(c(0, 60, 180), function(m) mc.fit(MC.ESTIMATORS$IV, mc.as.fast(dr0, pub, "vyp"), "MKT", m)))
check("T7 sigma = 0: the IV root is theta_0 within 1e-8 at H = 0, 5, 15 years", max(abs(r0$MKT - 1)) < 1e-8 && all(r0$n.roots >= 1),
  sprintf("max |beta - 1| %.1e; roots %s", max(abs(r0$MKT - 1)), paste(r0$n.roots, collapse = ",")))
# (c) estimate and SE against an independent computation on the pipeline's unit errors
rows.v <- sim.rows(dr, "vyp"); prep.v <- prepare.units(rows.v, fc)
for (m in c(0, 120)) {
  fit <- mc.fit(MC.ESTIMATORS$IV, fpv, "MKT", m)
  zz <- ref.errors(c(MKT = 1), prep.v, m) / ref.errors(c(MKT = 1), prep.v, 0)
  gb <- function(b) mean(zz * ref.errors(c(MKT = b), prep.v, 0))
  yrs <- fpv$vintage; V <- length(yrs); S <- ceiling(1.3 * sqrt(V)); s <- zz * ref.errors(c(MKT = fit$MKT), prep.v, 0)
  om <- 0; for (i in seq_len(V)) for (j in seq_len(V)) om <- om + max(0, 1 - abs(yrs[i] - yrs[j]) / S) * s[i] * s[j]
  Gd <- (gb(fit$MKT + 1e-4) - gb(fit$MKT - 1e-4)) / 2e-4
  se.ref <- sqrt(om / V / (V * Gd^2))
  check(sprintf("T7 H = %d years: IV moment zero at the estimate and SE = independent computation (pipeline errors, explicit kernel sum)", m / 12),
    abs(gb(fit$MKT)) < 1e-10 && abs(fit$SE.MKT / se.ref - 1) < 1e-6 && fit$hac.b == S / V,
    sprintf("gbar %.1e, SE %.5f vs %.5f, %d root(s)", gb(fit$MKT), fit$SE.MKT, se.ref, fit$n.roots))
}

# (d) no root in Theta (25 vintages of 5 funds, replication 4: the mean moment stays positive): the
# minimizer of gbar^2 is reported without SE, so it counts as not tested
nr <- mc.replication(mc.design(cons(5, 25), sigma = 0.13, max.holding.period = 11.5, path = path), 4, pub, "IV", "MKT", c(0, 180))
check("T7 no root of the moment in Theta: estimate = minimizer of gbar^2, no SE (not tested)",
  all(nr$n.roots == 0) && all(is.na(nr$SE.MKT)) && all(nr$fit.method == "min.square") && all(abs(nr$value) > 0.1),
  sprintf("|gbar| at the minimizer %s", paste(sprintf("%.3f", abs(nr$value)), collapse = ", ")))

# ---------------------------------------------------------------------------
# T9: Online Appendix D on the stored G1 samples (simulation/mc_appendix_d.R).
# (a) The reader's fast structure equals the pipeline's preparation of the same
# samples, built here the way estim_model_optimized.R builds it (unique Fund.ID
# per sample, all rows padded together to the file's common data end, public
# factors merged by date, simulated MKT and RF merged by (Date, sample_id)):
# same units, unit errors and L2 criterion within 1e-10. (b) The domain (the
# months of the sample's window with that sample's factor values) equals the
# pipeline's admissible months. (c) The alpha bounds equal the config
# loader's. (d) The inputs are those of the G1 table (md5 in the G1 results).
# Samples 1 and 2 of S2, S5 (simulated market) and Panel B S3 (truncated
# payoffs). With the same prepared data the estimator is the pipeline's (T3).
source(here("simulation/mc_appendix_d.R"))
source(here("helper/config_loader.R"))
ad.pub <- ad.public(); ad.fc <- unique(c("RF", setdiff(colnames(ad.pub), "Date")))
ad.cfg <- load_simulation_config(here(AD.CONFIG))
ad.pars <- list(c(MKT = 1, Alpha = 0), c(MKT = 0.7, Alpha = 0.004), c(MKT = 1.3, Alpha = -0.003), c(MKT = 2.5, Alpha = 0.01), c(MKT = -0.5, Alpha = 0.015))
for (id in c("big_n_v_50funds_alpha", "big_n_v_50funds_alpha_stdv30_shifted_mkt", "big_n_v_50funds_alpha_stdv30_additive")) {
  L <- ad.load(id)
  g1f <- list.files(sprintf("simulation/data_out_2026_g1/cache_q_factors_%s_simulated_cashflows_EW_VYP", id), "_m0_l0_cached_res\\.csv$", full.names = TRUE)
  g1 <- read.csv(g1f[length(g1f)]) # the latest run at H = 0 (the table's; S1 also has an earlier one)
  check(sprintf("T9 %s: input file is the G1 table's (md5)", id), identical(unique(g1$input.md5), L$md5), L$md5)
  ref <- as.data.frame(data.table::rbindlist(L$samples[c("1", "2")])); ref$sample <- NULL
  ref$Fund.ID <- as.factor(paste(ref$Fund.ID, ref$type, sep = "_"))
  ref <- pad.all.units(ref, "first.cf", L$t.max)
  df0 <- merge(ref, ad.pub, by = "Date", all.x = TRUE)
  if (!is.null(L$factors)) {
    sf <- read.csv(L$factors.file); sf$Date <- as.Date(sf$Date)
    df0$sample_id <- as.integer(sub(".*\\.(\\d+)$", "\\1", as.character(df0$type)))
    df0 <- merge(df0, sf[, c("Date", "sample_id", "MKT", "RF")], by = c("Date", "sample_id"), all.x = TRUE, suffixes = c(".hist", ".sim"))
    df0$MKT <- ifelse(is.na(df0$MKT.sim), df0$MKT.hist, df0$MKT.sim); df0$RF <- ifelse(is.na(df0$RF.sim), df0$RF.hist, df0$RF.sim)
    df0$MKT.sim <- df0$MKT.hist <- df0$RF.sim <- df0$RF.hist <- df0$sample_id <- NULL
  }
  df0 <- df0[df0$Vintage <= 2021, ]
  bounds <- unlist(scenario_to_estimation_params(ad.cfg$scenarios[[id]])[c("alpha_lower", "alpha_upper")])
  check(sprintf("T9 %s: alpha bounds = the config loader's", id), identical(unname(bounds), ad.alpha.bounds(id)), paste(ad.alpha.bounds(id), collapse = ", "))
  MC.SETTINGS$alpha.bounds <- ad.alpha.bounds(id)
  for (k in 1:2) {
    d.k <- df0[df0$type == L$samples[[as.character(k)]]$type[1], ]
    prep <- prepare.units(d.k, ad.fc)
    fp <- ad.prepare(L$samples[[as.character(k)]], L$t.max, ad.pub, if (is.null(L$factors)) NULL else L$factors[[as.character(k)]])
    worst <- c(errors = 0, L2 = 0)
    for (par in ad.pars) {
      a <- fast.unit.errors(par, fp, 0); b <- ref.errors(par, prep, 0)
      worst["errors"] <- max(worst["errors"], rel.diff(unname(a), unname(b)))
      worst["L2"] <- max(worst["L2"], rel.diff(sqrt(sum(a^2)) / length(a), sqrt(sum(b^2)) / length(b)))
    }
    check(sprintf("T9 %s sample %d: same units; unit errors and L2 criterion = pipeline preparation within 1e-10", id, k),
      identical(fp$fund.id, names(prep$units)) && max(worst) <= 1e-10, sprintf("%d units, max relative difference %.1e", fp$n.units, max(worst)))
    mo <- d.k[!duplicated(d.k$Date), c("Date", ad.fc)]; mo <- mo[order(mo$Date), ]
    cp <- admissible.constraints(mo, c("MKT", "Alpha"), MC.SETTINGS$gross.floor, MC.SETTINGS$mkt.bounds, MC.SETTINGS$alpha.bounds)
    cf <- mc.domain(fp, c("MKT", "Alpha"))
    check(sprintf("T9 %s sample %d: domain = the pipeline's admissible months (sample's own factor path)", id, k),
      identical(dim(cp$ui), dim(cf$ui)) && max(abs(cp$ui - cf$ui)) == 0 && max(abs(cp$ci - cf$ci)) == 0 && identical(cp$kind, cf$kind),
      sprintf("%d constraints; months %s to %s", nrow(cf$ui), format(min(mo$Date)), format(max(mo$Date))))
  }
}
ad.own <- sapply(AD.ROWS$scenario, ad.alpha.bounds)
ad.ref <- sapply(AD.ROWS$scenario, function(id) unname(unlist(scenario_to_estimation_params(ad.cfg$scenarios[[id]])[c("alpha_lower", "alpha_upper")])))
check("T9 all 19 scenarios: alpha bounds = the config loader's (S6 and S10 with upper bound 1%)", identical(unname(ad.own), unname(ad.ref)) &&
  all(ad.own[2, grepl("_mkt(2|8)", AD.ROWS$scenario)] == 0.01) && all(ad.own[2, !grepl("_mkt(2|8)", AD.ROWS$scenario)] == 0.02),
  sprintf("upper bounds %s", paste(sort(unique(ad.own[2, ])), collapse = ", ")))
MC.SETTINGS$alpha.bounds <- c(-0.02, 0.02)

# ---------------------------------------------------------------------------
# T10: path-correlated vintage shock (simulation/mc_generate.R, path.shock;
# log entry "hardening weak links 1 and 3"). (a) With a = 0 the draw equals
# the control's (same design without the shock); (b) at sigma = 0 the shock
# is constant within a vintage and equals a log R_v - K_v log mean(G^a),
# computed here independently from the path; (c) its mean over draws is one
# (unconditional pricing (P)); (d) given the path it is not one (X fails).
cons20 <- cons(20, 25)
ps <- function(a, sigma = 0.13) mc.design(cons20, sigma = sigma, max.holding.period = 11.5, path = path, path.shock = list(a = a, window.years = 10))
ctrl <- mc.design(cons20, sigma = 0.13, max.holding.period = 11.5, path = path)
check("T10 a = 0: the draw equals the control's (shock uses no random numbers)",
  identical(mc.draw(ps(0), 7, pub)$deals, mc.draw(ctrl, 7, pub)$deals) && identical(mc.draw(ps(0), 7, pub)$public, mc.draw(ctrl, 7, pub)$public))
d10 <- mc.draw(ps(0.615, sigma = 0), 3, pub)
G10 <- as.numeric(1 + d10$public$RF + d10$public$MKT)
src10 <- pub[pub$Date >= as.Date(path$from) & pub$Date <= as.Date(path$to), ]
lm10 <- log(mean((1 + src10$RF + src10$MKT)^0.615))
wv <- with(d10$deals, tapply(log(mult) - vapply(seq_along(start), function(i) sum(log(G10[seq_len(end[i] - start[i]) + start[i]])), 0), vintage, range))
w.ref <- vapply(as.integer(names(cons20)), function(v) { f <- match(as.Date(paste0(v, "-01-31")), d10$public$Date); K <- min(120L, nrow(d10$public) - f + 1L)
  0.615 * sum(log(G10[f:(f + K - 1L)])) - K * lm10 }, 0)
check("T10 sigma = 0: the shock is constant within a vintage and equals a log R_v - K_v log mean(G^a) (independent computation)",
  max(vapply(wv, function(r) diff(r), 0)) < 1e-10 && max(abs(vapply(wv, `[`, 0, 1) - w.ref)) < 1e-10,
  sprintf("max spread %.1e, max difference %.1e", max(vapply(wv, function(r) diff(r), 0)), max(abs(vapply(wv, `[`, 0, 1) - w.ref))))
ew <- vapply(1:400, function(sd) { d <- mc.draw(ps(0.615, sigma = 0), sd, pub); G <- as.numeric(1 + d$public$RF + d$public$MKT)
  i <- which(d$deals$vintage == 1990)[1]; exp(log(d$deals$mult[i]) - sum(log(G[seq_len(d$deals$end[i] - d$deals$start[i]) + d$deals$start[i]]))) }, 0)
check("T10 (P): the shock's mean over 400 draws is one within 4 Monte Carlo SEs (vintage 1990, a = 0.615)",
  abs(mean(ew) - 1) < 4 * sd(ew) / sqrt(length(ew)), sprintf("mean %.4f, SE %.4f, sd of log shock %.3f", mean(ew), sd(ew) / sqrt(length(ew)), sd(log(ew))))
check("T10 (X fails): given the path the shock is not one (spread across draws)", sd(log(ew)) > 0.25, sprintf("sd of the log shock %.3f (design target 0.306)", sd(log(ew))))

# T11: fund sizes and fund-size weighting (design B). (a) sizes are drawn last:
# the deals and the EW fast structures equal those of the design without
# sizes; (b) the FW fast structures equal prepare.fast() on the FW rows and the
# pipeline's unit errors, and the FW U criterion equals u.criterion() with
# sizes, within 1e-10; (c) the FW vintage pools equal pool.vintage() on the
# FW fund rows; (d) at sigma = 0 every FW unit is priced exactly at theta_0.
source(here("helper/pool_vintage.R"))
szd <- function(sigma = 0.13) mc.design(shape, sigma = sigma, max.holding.period = 11.5, sizes = list(sdlog = 1.19))
plain <- mc.design(shape, sigma = 0.13, max.holding.period = 11.5)
a11 <- mc.draw(szd(), 5, pub); b11 <- mc.draw(plain, 5, pub)
check("T11 sizes drawn last: deals and EW structures equal the design without sizes",
  identical(a11$deals, b11$deals) && identical(unclass(mc.as.fast(a11, pub, "fund")), unclass(mc.as.fast(b11, pub, "fund"))) &&
  identical(unclass(mc.as.fast(a11, pub, "vyp")), unclass(mc.as.fast(b11, pub, "vyp"))))
fw.rows <- function(dr, level) { d <- mc.as.rows(dr, level, fw = TRUE); d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  d <- pad.all.units(d, "first.cf", max(d$Date)); d$Fund.ID <- as.factor(as.character(d$Fund.ID)); merge(d, pub, by = "Date", all.x = TRUE) }
for (lev in c("fund", "vyp")) {
  rows <- fw.rows(a11, lev); fpr <- prepare.fast(rows, fc); prr <- prepare.units(rows, fc); fpd <- mc.as.fast(a11, pub, paste0(lev, ".fw"))
  worst <- c(direct = 0, pipeline = 0)
  for (par in pars[1:5]) for (m in months) { a <- fast.unit.errors(par, fpd, m); b <- fast.unit.errors(par, fpr, m); r <- ref.errors(par, prr, m)
    worst["direct"] <- max(worst["direct"], rel.diff(unname(a), unname(b))); worst["pipeline"] <- max(worst["pipeline"], rel.diff(unname(a), unname(r))) }
  check(sprintf("T11 FW %s units: direct structure = prepare.fast() on the FW rows and = pipeline unit errors within 1e-10", lev),
    fpd$n.units == fpr$n.units && max(worst) <= 1e-10, sprintf("%.1e / %.1e", worst[["direct"]], worst[["pipeline"]]))
  if (lev == "fund") {
    s.ref <- vapply(prr$units, `[[`, 0, "size"); key.ref <- vapply(prr$units, `[[`, "", "key")
    check("T11 FW U criterion (fast, sizes) = u.criterion() with the pipeline's sizes within 1e-10; sizes match the rows",
      rel.diff(fast.u.criterion(fast.unit.errors(c(MKT = 1.2), fpd, 120), fpd, fpd$size), u.criterion(ref.errors(c(MKT = 1.2), prr, 120), key.ref, s.ref)) <= 1e-10 &&
      max(abs(unname(fpd$size) - unname(s.ref))) < 1e-12)
  }
}
fr <- mc.as.rows(a11, "fund", fw = TRUE); fr$NAV <- 0; fr$CON <- 0; fr$DIS <- 0
pv <- do.call(rbind, lapply(split(fr, fr$Vintage), pool.vintage, fund.size.weighting = TRUE))
vr <- mc.as.rows(a11, "vyp", fw = TRUE)
mm <- merge(pv[, c("Vintage", "Date", "CF")], vr[, c("Vintage", "Date", "CF")], by = c("Vintage", "Date"), all = TRUE)
check("T11 FW vintage pools = pool.vintage(fund.size.weighting = TRUE) on the FW fund rows",
  !anyNA(mm) && max(abs(mm$CF.x - mm$CF.y)) < 1e-12, sprintf("%d dates, max difference %.1e", nrow(mm), max(abs(mm$CF.x - mm$CF.y))))
z11 <- mc.draw(szd(sigma = 0), 9, pub)
check("T11 sigma = 0: every FW fund and FW vintage portfolio is priced exactly at theta_0 (H = 0)",
  max(abs(fast.unit.errors(c(MKT = 1), mc.as.fast(z11, pub, "fund.fw"), 0))) < 1e-10 && max(abs(fast.unit.errors(c(MKT = 1), mc.as.fast(z11, pub, "vyp.fw"), 0))) < 1e-10)
r11 <- mc.replication(szd(), 3, pub, c("L2.vyp.fw", "U.eq.fw"), "MKT", c(0, 180))
check("T11 the driver fits FW L2 and FW U (interior, with SEs)", nrow(r11) == 4 && all(is.finite(r11$MKT)) && all(is.finite(r11$SE.MKT[!r11$at.bound])))

t.gen <- system.time(for (s in 1:20) { x <- mc.draw(mc.design(shape), s, pub); mc.as.fast(x, pub, "fund"); mc.as.fast(x, pub, "vyp") })[["elapsed"]] / 20
cat(sprintf("Generator, Buyout shape (%d funds, %d deals): %.2f s per replication incl. both fast structures; %.1f MB in memory\n",
  sum(shape), sum(shape) * 15, t.gen, as.numeric(object.size(list(dr, mc.as.fast(dr, pub, "fund")))) / 2^20))

cat("\n", n.fail, "test(s) failed\n")
if (n.fail > 0) quit(status = 1)

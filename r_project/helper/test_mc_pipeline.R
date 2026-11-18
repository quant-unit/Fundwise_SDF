# Acceptance test T3 of the Step 4 Monte Carlo (revision plan, Section I,
# implementation spec): the driver's estimates (simulation/mc_driver.R, fast
# criterion) equal the pipeline's (run_simulation_study() -> run_estimation(),
# estim_model_optimized.R) on the same replications, within 1e-6. Slow (the
# pipeline fits every model on fund-level data): several minutes.
# Run from r_project/: LANG=en_US.UTF-8 Rscript helper/test_mc_pipeline.R
suppressPackageStartupMessages(library(here))
source("simulation/mc_driver.R")
n.fail <- 0L
check <- function(label, ok, detail = "") {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", label, if (nzchar(detail)) paste0("(", detail, ")"), "\n")
  if (!isTRUE(ok)) n.fail <<- n.fail + 1L
}

seeds <- 1:2
horizons <- c(0, 120, 360)
models <- c("MKT", "Alpha", "EG")
pub <- mc.public.factors()
design <- mc.design(mc.data.shape("BO"), sigma = 0.2)

# 1) the replications as pipeline input files (one type per replication)
folder <- "t3_mc"
dir <- file.path("simulation/data_prepared_sim", folder)
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
draws <- lapply(seeds, function(s) mc.draw(design, s, pub))
for (lev in c("fund", "vyp")) {
  rows <- do.call(rbind, lapply(seq_along(seeds), function(i) mc.as.rows(draws[[i]], lev, type = paste0("SIM.", seeds[i]))))
  write.csv(rows, file.path(dir, paste0(folder, "_simulated_cashflows_", if (lev == "fund") "EW" else "EW_VYP", ".csv")), row.names = FALSE)
}

# 2) the pipeline on these files: one scenario per estimator
cfg <- file.path(tempdir(), "t3_mc.yaml")
scen <- function(id, weighting, error, vw = "equal") sprintf(
  "  %s:\n    estimation: {weighting: \"%s\", error_function: \"%s\", vintage_weights: \"%s\", factors_to_use: \"\"}\n    active: true\n    data_folder: %s\n",
  id, weighting, error, vw, folder)
writeLines(c("defaults:", "  dgp: {beta_MKT: 1, no_samples: 2, no_funds: 1, min_vin: 1985, max_vin: 2010} # required by validate_scenario(); data are not generated here",
  "  estimation:", "    sdf_model: \"linear\"", sprintf("    max_months: [%s]", paste(horizons, collapse = ", ")), "    lambdas: 0",
  "    kernel_bandwidth: \"llsw\"", "    alpha_lower: -0.02", "    alpha_upper: 0.02", "    do_cross_validation: false",
  "    part_to_keep: 1", "    no_partitions: 1", "    max_vintage: 2021", "scenarios:",
  scen("t3_L2.vyp", "EW_VYP", "L2_Lasso"), scen("t3_L2.fund", "EW", "L2_Lasso"),
  scen("t3_U.eq", "EW", "U_pair"), scen("t3_U.n", "EW", "U_pair", "n.funds")), cfg)
out <- file.path(tempdir(), "t3_out")
unlink(out, recursive = TRUE)
source("run_simulation_study.R")
t.pipe <- system.time(invisible(capture.output(res <- run_simulation_study(config_path = cfg, generate_data = FALSE, estimate = TRUE, data_out_folder = out))))[["elapsed"]]
check("T3 pipeline ran all four estimators", all(vapply(res, function(x) isTRUE(x$success), TRUE)) && length(res) == 4)
pipe <- do.call(rbind.all.columns.mc, lapply(list.dirs(out, recursive = FALSE), function(d) {
  x <- do.call(rbind.all.columns.mc, lapply(list.files(d, "_cached_res\\.csv$", full.names = TRUE), read.csv, stringsAsFactors = FALSE))
  x$estimator <- if (grepl("U_pair_qn.funds", d)) "U.n" else if (grepl("U_pair", d)) "U.eq" else if (grepl("EW_VYP", d)) "L2.vyp" else "L2.fund"
  x
}))
pipe <- pipe[pipe$Factor %in% models, ]
pipe$rep <- as.integer(sub("SIM\\.", "", pipe$Type))
pipe$Coef <- ifelse(pipe$Factor == "MKT", NA, vapply(seq_len(nrow(pipe)), function(i) if (pipe$Factor[i] == "MKT") NA_real_ else pipe[[pipe$Factor[i]]][i], 0))
pipe$SE.Coef <- ifelse(pipe$Factor == "MKT", NA, vapply(seq_len(nrow(pipe)), function(i) if (pipe$Factor[i] == "MKT") NA_real_ else pipe[[paste0("SE.", pipe$Factor[i])]][i], 0))

# 3) the driver on the same draws
t.mc <- system.time(mc <- do.call(rbind, lapply(seeds, function(s) mc.replication(design, s, pub, models = models, horizons = horizons))))[["elapsed"]]
m <- merge(mc, pipe, by.x = c("rep", "estimator", "model", "max.month"), by.y = c("rep", "estimator", "Factor", "max.month"), suffixes = c(".mc", ".pipe"))
check("T3 every driver row has a pipeline row", nrow(m) == nrow(mc) && nrow(mc) == length(seeds) * 4 * length(models) * length(horizons), sprintf("%d of %d", nrow(m), nrow(mc)))
d.est <- pmax(abs(m$MKT.mc - m$MKT.pipe), abs(m$Coef.mc - m$Coef.pipe), na.rm = TRUE)
check("T3 estimates equal the pipeline's within 1e-6", max(d.est) <= 1e-6, sprintf("max |difference| %.1e", max(d.est)))
check("T3 bound flags equal", identical(as.logical(m$at.bound.mc), as.logical(m$at.bound.pipe)), sprintf("%d of %d rows on a bound", sum(m$at.bound.mc), nrow(m)))
int <- !m$at.bound.mc
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-12)
# SEs at the driver's own estimates differ from the pipeline's only through
# the optimizer end points (about 1e-7 apart; a Wald statistic of beta = 1 has
# relative sensitivity 2 d(theta) / (theta - 1), large near the null), so the
# sandwich code is compared at the pipeline's estimates themselves
d.se.own <- pmax(rel(m$SE.MKT.mc, m$SE.MKT.pipe), rel(m$SE.Coef.mc, m$SE.Coef.pipe), na.rm = TRUE)[int]
fps <- lapply(draws, function(dr) list(fund = mc.as.fast(dr, pub, "fund"), vyp = mc.as.fast(dr, pub, "vyp")))
at.pipe <- do.call(rbind, lapply(which(int), function(i) {
  est <- MC.ESTIMATORS[[m$estimator[i]]]
  par <- if (m$model[i] == "MKT") c(MKT = m$MKT.pipe[i]) else setNames(c(m$MKT.pipe[i], m$Coef.pipe[i]), c("MKT", m$model[i]))
  cv <- mc.sandwich(est, fps[[match(m$rep[i], seeds)]][[est$level]], par, m$max.month[i])
  data.frame(se.mkt = cv$SE[["SE.MKT"]], se.coef = if (length(par) > 1) cv$SE[[2]] else NA_real_, wald = cv$Wald.stat, S = cv$hac.S, V = cv$hac.V)
}))
# Tolerance of the derived quantities: the SEs come from finite-difference
# Hessians (delta = 1e-4), which magnify criterion differences by about
# 1 / delta^2 = 1e8; with criteria equal to about 1e-13 (T1) the SEs can agree
# only to about 1e-5 (the first run with 1e-6 failed at 1.8e-6 for that reason)
tol.se <- 1e-5
d.se <- pmax(rel(at.pipe$se.mkt, m$SE.MKT.pipe[int]), rel(at.pipe$se.coef, m$SE.Coef.pipe[int]), na.rm = TRUE)
check("T3 sandwich SEs at the pipeline's estimates equal the pipeline's within 1e-5 (relative)", max(d.se) <= tol.se, sprintf("max %.1e; at the driver's own estimates %.1e", max(d.se), max(d.se.own)))
check("T3 bandwidth and vintage count equal", all(at.pipe$S == m$hac.S.pipe[int]) && all(at.pipe$V == m$hac.V.pipe[int]))
w <- rel(at.pipe$wald, m$Wald.stat.MKT_1[int])
check("T3 Wald statistics at the pipeline's estimates equal within 1e-5 (relative)", max(w, na.rm = TRUE) <= tol.se, sprintf("max %.1e", max(w, na.rm = TRUE)))
by.est <- aggregate(d.est ~ estimator, data.frame(estimator = m$estimator, d.est = d.est), max)
cat("largest estimate difference by estimator:", paste(sprintf("%s %.1e", by.est$estimator, by.est$d.est), collapse = "; "), "\n")
cat(sprintf("time: pipeline %.0f s (with jackknife and validation errors), driver %.0f s, %d fits\n", t.pipe, t.mc, nrow(mc)))

unlink(dir, recursive = TRUE); unlink(out, recursive = TRUE)
cat("\n", n.fail, "test(s) failed\n")
if (n.fail > 0) quit(status = 1)

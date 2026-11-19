# Step 4: Online Appendix D (the DLP12 comparison, scenarios S1-S10, both
# panels) under the r5 estimator (revision plan, log entries "S1--S10 of
# Online Appendix D: scope" and "S1--S10 under r5").
# Run from r_project/ with a clean, committed tree:
#   LANG=en_US.UTF-8 Rscript simulation/mc_appendix_d.R [scenario ...]
#
# Route (C): the G1 simulated samples (simulation/data_prepared_sim/, the
# pipeline's EW vintage-portfolio cash flows of 16-17 November, generated after
# the B1 timing repair) are re-estimated with the Step 4 estimator (mc.fit():
# L2 on vintage portfolios, (alpha, beta), H = 0, the r5 domain, the Bartlett
# sandwich, bound flags). Each sample is prepared as estim_model_optimized.R
# prepares a simulation input: Fund.ID made unique per sample, the monthly grid
# from the vintage origin (first cash flow of the sample's vintage) to the
# common data end of the file, the public factors merged by date, MKT and RF
# replaced by the sample's simulated path where the scenario simulates the
# market. The alpha bounds are each scenario's (config/simulation_scenarios.yaml,
# with the config default); the beta range and the gross floor are r5's.
# Acceptance test: T9 in helper/test_mc.R.
#
# Output: simulation/data_out_mc/appendix_d/<scenario>.rds (one row per sample,
# provenance; a finished scenario is skipped) and
# results/mc_step4/appendix_d_summary.csv (one row per scenario), read by
# simulation/mc_paper.R for the table.
source("simulation/mc_driver.R")
suppressPackageStartupMessages(library(data.table))

AD.DIR <- "simulation/data_out_mc/appendix_d"
AD.SIM <- "simulation/data_prepared_sim"
AD.CONFIG <- "config/simulation_scenarios.yaml"
AD.MAX.VINTAGE <- 2021 # estimation default max_vintage of the simulation config

# The 19 rows: Panel A (multiplicative shocks; S4 equals S3 there and is
# omitted) and Panel B (the "_additive" twins), labels as in the G1 table
AD.BASE <- c(S1 = "base_case_zero_alpha", S2 = "big_n_v_50funds_alpha", S3 = "big_n_v_50funds_alpha_stdv30",
  S4 = "big_n_v_50funds_alpha_stdv30_shifted", S5 = "big_n_v_50funds_alpha_stdv30_shifted_mkt",
  S6 = "big_n_v_50funds_alpha_stdv30_shifted_mkt2", S7 = "big_n_v_50funds_alpha_stdv30_shifted_mkt4",
  S8 = "big_n_v_50funds_alpha_stdv30_shifted_mkt6", S9 = "big_n_v_50funds_alpha_stdv30_shifted_mkt7",
  S10 = "big_n_v_50funds_alpha_stdv30_shifted_mkt8")
AD.ROWS <- rbind(data.frame(panel = "A", row = setdiff(names(AD.BASE), "S4"), scenario = AD.BASE[setdiff(names(AD.BASE), "S4")]),
                 data.frame(panel = "B", row = names(AD.BASE), scenario = paste0(AD.BASE, "_additive")))
rownames(AD.ROWS) <- NULL

# Each scenario's alpha bounds: the scenario's estimation setting, else the
# config default, else none (as helper/config_loader.R resolves them)
ad.alpha.bounds <- function(id, config = yaml::read_yaml(AD.CONFIG)) {
  est <- config$scenarios[[id]]$estimation; def <- config$defaults$estimation
  b <- function(k, dflt) { x <- est[[k]] %||% def[[k]]; if (is.null(x)) dflt else as.numeric(x) }
  c(b("alpha_lower", -Inf), b("alpha_upper", Inf))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Public factors as estim_model_optimized.R loads them (dates normalized to
# month ends, Alpha column)
ad.public <- function(file = "empirical/data_prepared_2026/q_factors.csv") {
  p <- read.csv(file)
  colnames(p) <- gsub("_World", "", colnames(p))
  p$Date <- as.Date(p$Date)
  p$Date <- as.Date(format(p$Date, "%Y-%m-01"))
  p$Date <- seq.Date(p$Date[1], by = "month", length.out = nrow(p))
  p$Date <- as.Date(format(p$Date + 32, "%Y-%m-01")) - 1
  p$Alpha <- 1
  p
}

# All samples of one scenario: the VYP cash flows split by sample, the common
# data end, and the simulated factor paths (NULL on the historical path)
ad.load <- function(id) {
  f <- file.path(AD.SIM, id, paste0(id, "_simulated_cashflows_EW_VYP.csv"))
  x <- fread(f, showProgress = FALSE)
  x[, Date := as.Date(Date)]
  x[, sample := as.integer(sub(".*\\.(\\d+)$", "\\1", type))]
  ff <- file.path(AD.SIM, id, paste0(id, "_simulated_factors.csv"))
  # read.csv as the pipeline (fread parses a few decimals one unit in the last place differently)
  sf <- if (file.exists(ff)) { y <- read.csv(ff); y$Date <- as.Date(y$Date); split(y, y$sample_id) } else NULL
  list(id = id, file = f, md5 = unname(tools::md5sum(f)), t.max = max(x$Date),
       samples = split(as.data.frame(x), x$sample), factors = sf,
       factors.file = if (is.null(sf)) NA_character_ else ff, factors.md5 = if (is.null(sf)) NA_character_ else unname(tools::md5sum(ff)))
}

# One sample as the pipeline prepares it, as a fast structure
ad.prepare <- function(d, t.max, public, sim.factors = NULL) {
  d$sample <- NULL
  d$Fund.ID <- paste(d$Fund.ID, d$type, sep = "_")
  d <- pad.all.units(d, "first.cf", t.max)
  if (!is.null(sim.factors)) { # simulated MKT (and RF) where the path has them, historical elsewhere
    i <- match(as.numeric(public$Date), as.numeric(sim.factors$Date))
    for (k in intersect(c("MKT", "RF"), colnames(sim.factors))) public[[k]] <- ifelse(is.na(i), public[[k]], sim.factors[[k]][i])
  }
  d <- merge(d, public, by = "Date", all.x = TRUE)
  d <- d[d$Vintage <= AD.MAX.VINTAGE, ]
  prepare.fast(d, unique(c("RF", setdiff(colnames(public), "Date"))))
}

# Fit every sample of one scenario (forked workers; the alpha bounds are set
# in MC.SETTINGS before forking, so every worker uses them)
ad.run <- function(id, public, workers = n.parallel.workers(1000)) {
  t0 <- Sys.time()
  L <- ad.load(id)
  MC.SETTINGS$alpha.bounds <<- ad.alpha.bounds(id)
  ks <- as.integer(names(L$samples))
  fit1 <- function(k) tryCatch({
    fp <- ad.prepare(L$samples[[as.character(k)]], L$t.max, public, if (is.null(L$factors)) NULL else L$factors[[as.character(k)]])
    data.frame(sample = k, mc.fit(MC.ESTIMATORS$L2.vyp, fp, "Alpha", 0), n.units = fp$n.units, error = NA_character_)
  }, error = function(e) data.frame(sample = k, error = conditionMessage(e)))
  rows <- do.call(rbind.all.columns.mc, parallel::mclapply(ks, fit1, mc.cores = workers, mc.preschedule = FALSE))
  rows <- rows[order(rows$sample), ]
  prov <- list(git.commit = system("git rev-parse --short HEAD", intern = TRUE),
               git.dirty = length(system("git status --porcelain --untracked-files=no", intern = TRUE)) > 0,
               input.file = L$file, input.md5 = L$md5, factors.file = L$factors.file, factors.md5 = L$factors.md5,
               t.max = L$t.max, settings = MC.SETTINGS, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  list(rows = rows, provenance = prov)
}

# One row per scenario: the columns of the G1 table (H = 0), the share of
# estimates with an active constraint of Theta (any, and alpha's upper bound)
ad.summary <- function(res, row) {
  r <- res$rows; ok <- is.na(r$error)
  r <- r[ok, ]; a <- 100 * r$Coef # monthly alpha in percent
  data.frame(panel = row$panel, row = row$row, scenario = row$scenario, samples = nrow(res$rows), failed = sum(!ok),
    true.beta = 1, mean.beta = mean(r$MKT), sd.beta = sd(r$MKT), p25.beta = unname(quantile(r$MKT, 0.25)), p75.beta = unname(quantile(r$MKT, 0.75)),
    true.alpha = 0, mean.alpha = mean(a), sd.alpha = sd(a), p25.alpha = unname(quantile(a, 0.25)), p75.alpha = unname(quantile(a, 0.75)),
    on.bound = mean(r$at.bound), on.alpha.upper = mean(grepl("Alpha.upper", r$active)), on.floor = mean(grepl("floor", r$active)),
    on.beta.bound = mean(grepl("MKT\\.", r$active)),
    alpha.lower = res$provenance$settings$alpha.bounds[1], alpha.upper = res$provenance$settings$alpha.bounds[2],
    seconds = res$provenance$seconds, git.commit = res$provenance$git.commit, git.dirty = res$provenance$git.dirty,
    input.md5 = res$provenance$input.md5)
}

if (sys.nframe() == 0L) {
  ids <- commandArgs(TRUE)
  if (!length(ids)) ids <- AD.ROWS$scenario
  stopifnot(all(ids %in% AD.ROWS$scenario))
  dir.create(AD.DIR, recursive = TRUE, showWarnings = FALSE)
  public <- ad.public()
  for (id in ids) {
    out <- file.path(AD.DIR, paste0(id, ".rds"))
    if (file.exists(out)) { cat(format(Sys.time(), "%H:%M:%S"), id, "exists, skipped\n"); next }
    res <- ad.run(id, public)
    saveRDS(res, out)
    cat(sprintf("%s %-52s %4d samples, %d failed, %4.0f s, commit %s%s\n", format(Sys.time(), "%H:%M:%S"), id, nrow(res$rows),
      sum(!is.na(res$rows$error)), res$provenance$seconds, res$provenance$git.commit, if (res$provenance$git.dirty) " (DIRTY)" else ""))
  }
  # the summary of every finished scenario
  done <- AD.ROWS[file.exists(file.path(AD.DIR, paste0(AD.ROWS$scenario, ".rds"))), ]
  s <- do.call(rbind, lapply(seq_len(nrow(done)), function(i) ad.summary(readRDS(file.path(AD.DIR, paste0(done$scenario[i], ".rds"))), done[i, ])))
  dir.create("results/mc_step4", showWarnings = FALSE, recursive = TRUE)
  write.csv(s, "results/mc_step4/appendix_d_summary.csv", row.names = FALSE)
  cat("summary:", nrow(s), "of", nrow(AD.ROWS), "rows\n")
}

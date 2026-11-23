# Integration regressions for the cache/reporting review fixes.
# Run from r_project/: LANG=en_US.UTF-8 Rscript helper/test_review_fixes.R
# All data, caches and plots are synthetic and live in a temporary project.
repo <- normalizePath(getwd())
stopifnot(file.exists(file.path(repo, "estim_model_optimized.R")))
stage <- tempfile("review-fixes-")
dir.create(stage)
for (dir in c("helper", "empirical/data_prepared_2026", "empirical/fixture", "config", "simulation")) {
  dir.create(file.path(stage, dir), recursive = TRUE, showWarnings = FALSE)
}
file.copy(list.files(file.path(repo, "helper"), pattern = "\\.R$", full.names = TRUE), file.path(stage, "helper"))
file.copy(file.path(repo, "simulation/analyze_simulation_bias.R"), file.path(stage, "simulation"))
file.create(file.path(stage, "estim_model_optimized.R"))
writeLines("Date,RF,MKT\n2000-01-31,0.001,0.01\n2000-02-29,0.001,0.02", file.path(stage, "empirical/data_prepared_2026/q_factors.csv"))
writeLines("Fund.ID,Date,CF,type,Vintage\nA,2000-01-31,-1,BO,2000\nA,2000-02-29,2,BO,2000", file.path(stage, "empirical/fixture/preqin_cashflows_EW_NAV.csv"))
config.file <- file.path(stage, "config/custom.yaml")
writeLines("max_vintage: 2010", config.file)
setwd(stage)
source("helper/estim_model_core.R")
n.fail <- 0L
check <- function(label, ok) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", label, "\n")
  if (!isTRUE(ok)) n.fail <<- n.fail + 1L
}

# Replace only the numerical engine; exercise the actual runner, identity,
# cache writer metadata and resume decision. Numerical wiring is tested in
# helper/test_code_repairs.R and in the real-data smoke run.
engine.calls <- 0L
fixture.beta <- 2
real.source <- source
fixture.engine <- function(file, ...) {
  engine.calls <<- engine.calls + 1L
  meta <- get("result.identity.used", envir = .GlobalEnv)
  df.res <<- cbind(data.frame(MKT = fixture.beta, Factor = "MKT", Type = "BO", value = 1,
    max.month = max.months[1], lambda = lambdas[1], CV.key = "ALL",
    datetime = "2026-11-16 12:00:00", SE.MKT = .01, at.bound = FALSE, git.commit = "0012345"), meta)
  path <- file.path(data.out.folder, paste0("cache_", public.filename, "_", cache.folder.tag))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  write.csv(df.res, file.path(path, paste0(engine.calls, "_cached_res.csv")), row.names = FALSE)
  invisible(NULL)
}
source <- function(file, ...) {
  if (identical(file, "estim_model_optimized.R")) fixture.engine(file, ...) else real.source(file, ...)
}
params <- list(data_prepared_folder = "fixture", weighting = "EW", factors_to_use = "MKT",
  max_vintage = 2010, max_months = 0, data_out_folder = file.path(stage, "results"),
  config_path = config.file, do_parallel = FALSE, verbose = FALSE)
r <- do.call(run_estimation, params)
cache <- file.path(params$data_out_folder, paste0("cache_q_factors_", r$cache_folder_tag))
check("fresh estimation writes readable current results", r$success && nrow(read.result.cache(cache)) == 1)
check("cache provenance preserves numeric-only commit identifiers", identical(read.result.cache(cache)$git.commit, "0012345"))
invisible(do.call(run_estimation, params))
check("unchanged estimation resumes", engine.calls == 1)
params$max_vintage <- 2005
invisible(do.call(run_estimation, params))
check("changed vintage cutoff is re-estimated", engine.calls == 2)
ambiguous <- tryCatch(read.result.cache(cache, verbose = FALSE), error = function(e) e)
check("analysis refuses to mix different current estimator settings", inherits(ambiguous, "error"))
params$kernel_bandwidth <- 6
invisible(do.call(run_estimation, params))
check("changed inference settings are re-estimated", engine.calls == 3)
public.file <- file.path(stage, "empirical/data_prepared_2026/q_factors.csv")
cat("2000-03-31,0.001,0.03\n", file = public.file, append = TRUE)
invisible(do.call(run_estimation, params))
check("changed public factors are re-estimated", engine.calls == 4)
writeLines("max_vintage: 2005", config.file)
invisible(do.call(run_estimation, params))
check("editing the configuration file without changing settings resumes (config is provenance only)", engine.calls == 4)
stored <- do.call(rbind, lapply(list.files(cache, "_cached_res\\.csv$", full.names = TRUE), read.csv, stringsAsFactors = FALSE))
check("result files store data paths relative to the project root",
  !any(is.absolute.path(na.omit(c(stored$input.path, stored$public.path, stored$config.path)))))
# a row as written before 18 November 2026: absolute paths, configuration in the hash
current <- read.result.cache(cache, verbose = FALSE)
legacy <- current[1, ]
legacy$input.path <- file.path(normalizePath(stage), legacy$input.path)
legacy$public.path <- file.path("/elsewhere/clone/r_project", legacy$public.path)
legacy$cache.identity <- "former-hash-including-config"
legacy <- current.result.rows(legacy, verbose = FALSE)
check("rows with absolute paths and the former identity stay current with the same identity",
  nrow(legacy) == 1 && identical(legacy$cache.identity, current$cache.identity[1]) && identical(legacy$input.path, current$input.path[1]))
summary.file <- file.path(stage, "saved-summary.csv")
write.csv(df.res, summary.file, row.names = FALSE)
check("current summaries remain readable", nrow(read.result.summary(summary.file)) == 1)
input.file <- file.path(stage, "empirical/fixture/preqin_cashflows_EW_NAV.csv")
cat("A,2000-03-31,0,BO,2000\n", file = input.file, append = TRUE)
check("analysis rejects rows for changed current input", nrow(read.result.cache(cache, verbose = FALSE)) == 0)
stale.summary <- tryCatch(read.result.summary(summary.file), error = function(e) e)
check("saved summaries cannot bypass input validation", inherits(stale.summary, "error"))
unlink(input.file)
check("analysis rejects rows with missing current input", nrow(read.result.cache(cache, verbose = FALSE)) == 0)
missing <- tryCatch(do.call(run_estimation, params), error = function(e) e)
check("missing input stops estimation before resume", inherits(missing, "error"))

# Use the real simulation runner and bias reader. Only its numerical engine
# is replaced above; the stored cache namespaces and metadata are real.
real.source(file.path(repo, "run_simulation_study.R"))
check("qualified cache tags stay canonical when partitioned",
  result.cache.tag("base_U_pair_qn.funds_origin_predetermined", "U_pair", "predetermined", "n.funds", 1) ==
    "base_part1_U_pair_qn.funds_origin_predetermined")
cfg <- yaml::read_yaml(file.path(repo, "config/simulation_scenarios.yaml"))
cfg$scenarios <- list(namespace = list(active = TRUE, data_folder = "auto",
  estimation = list(weighting = "EW", error_function = "U_pair", origin_convention = "predetermined",
    vintage_weights = "n.funds", factors_to_use = "MKT", max_months = 0)))
yaml::write_yaml(cfg, config.file)
dir.create("simulation/data_prepared_sim/namespace", recursive = TRUE)
writeLines("Fund.ID,Date,CF,type,Vintage\nA,2000-01-31,-1,S1,2000", "simulation/data_prepared_sim/namespace/namespace_simulated_cashflows_EW.csv")
sc <- load_simulation_config(config.file)$scenarios$namespace
p <- scenario_to_estimation_params(sc)
p$config_path <- config.file
p$data_out_folder <- file.path(stage, "simulation-results")
p$verbose <- FALSE
baseline <- p
baseline$error_function <- "L2_Lasso"; baseline$origin_convention <- "first.cf"; baseline$vintage_weights <- "equal"
fixture.beta <- 9
invisible(do.call(run_estimation, baseline))
fixture.beta <- 2
correct <- do.call(run_estimation, p)
invisible(capture.output(sim <- run_simulation_study(config_path = config.file, estimate = FALSE,
  analyze = TRUE, data_out_folder = p$data_out_folder, verbose = FALSE)))
check("analysis-only selects the requested estimator cache", identical(as.numeric(attr(sim, "bias_summary")$est_beta_MKT), 2))
cfg$scenarios$namespace$estimation$no_partitions <- 2
yaml::write_yaml(cfg, config.file)
invisible(capture.output(parts <- run_simulation_study(config_path = config.file, estimate = TRUE,
  analyze = TRUE, data_out_folder = p$data_out_folder, verbose = FALSE)))
check("partitioned analysis reads every option-specific partition", nrow(attr(parts, "bias_summary")) == 2)
part.cache <- file.path(p$data_out_folder, paste0("cache_q_factors_", parts[[1]]$cache_folder_tags[1]))
factor.file <- "simulation/data_prepared_sim/namespace/namespace_simulated_factors.csv"
writeLines("Date,sample_id,MKT,RF\n2000-01-31,1,0.02,0.001", factor.file)
check("adding a simulated factor override invalidates existing results", nrow(read.result.cache(part.cache, verbose = FALSE)) == 0)
unlink(factor.file)

# Report fixtures pass through real analysis and the real table exporters.
source <- real.source
suppressPackageStartupMessages(library(dplyr))
here <- function(...) {p <- file.path(...); if (startsWith(p, "/")) p else file.path(stage, p)}
file.copy(file.path(repo, "results"), stage, recursive = TRUE)
source.internally <- FALSE
prefix <- "q_factors_preqin_"; suffix <- "EW_VYP"
data.out.folder <- file.path(stage, "reporting")
report.cache <- file.path(data.out.folder, paste0("cache_", prefix, suffix))
dir.create(report.cache, recursive = TRUE)
meta <- make.result.identity("simulation/data_prepared_sim/namespace/namespace_simulated_cashflows_EW.csv",
  public.file, list(test = "report"), config.file)
rows <- do.call(rbind, lapply(c("MKT", "Alpha"), function(fac) do.call(rbind,
  lapply(c("ALL", "F1", "F2"), function(fold) cbind(data.frame(
    MKT = 1.2, Alpha = if (fac == "Alpha") .02 else NA_real_,
    SE.MKT = .01, SE.MKT.indep = .015, SE.Alpha = .001, SE.Alpha.indep = .0015,
    Factor = fac, Type = "BO", max.month = 1, lambda = 0, CV.key = fold,
    datetime = "2026-11-16 12:00:00", at.bound = fold == "ALL", validation.error = .1, validation.error.H0 = .1,
    git.commit = if (fold == "ALL") "ALL_COMMIT" else "CV_COMMIT", git.dirty = FALSE,
    scenario.id = if (fold == "ALL") "ALL_SCENARIO" else "CV_SCENARIO", data.prepared.folder = "fixture"), meta)))))
write.csv(rows, file.path(report.cache, "report_cached_res.csv"), row.names = FALSE)
pdf(file.path(stage, "analysis.pdf")); plot(0, 0, type = "n")
invisible(capture.output(source(file.path(repo, "results/analyze_result.R"))))
dev.off()
cv <- read.csv(file.path(report.cache, "0_cross_validation_summary.csv"))
check("CV summaries preserve contributing commits and scenarios",
  all(c("git.commit", "scenario.id") %in% names(cv)) && all(grepl("CV_COMMIT", cv$git.commit)) && all(grepl("CV_SCENARIO", cv$scenario.id)))
for (ex in parse(file.path(repo, "results/plot_empirical_estimates.R"))) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.call(ex[[3]]) && identical(ex[[3]][[1]], as.name("function"))) eval(ex, .GlobalEnv)
}
asym <- read.csv(file.path(report.cache, "0_asymptotic_inference_summary.csv"))
asym$weighting <- "EW"; asym$method <- "Asymptotic"; asym$vintage <- 2010
cv$weighting <- "EW"; cv$method <- "Cross-Validation"
all.rows <- bind_rows(asym, cv)
.export_latex_tables(all.rows, "BO", c("MKT", "Alpha"), file.path(stage, "horizon"))
horizon <- readLines(file.path(stage, "horizon_tables.tex"))
check("horizon headers include full-sample and CV provenance", any(grepl("CV_COMMIT", horizon, fixed = TRUE)) && any(grepl("ALL_COMMIT", horizon, fixed = TRUE)))
.export_vintage_latex_tables(asym, "BO", c("MKT", "Alpha"), 2010, file.path(stage, "cutoff"))
cutoff <- readLines(file.path(stage, "cutoff_tables.tex"))
check("cutoff tables mark the boundary coefficient cells", sum(grepl("^0.1 yr &", cutoff) & grepl("dagger", cutoff, fixed = TRUE)) == 3)
direct <- readLines(file.path(report.cache, "empirical_tables.tex"))
check("direct analysis tables mark boundaries and suppress invalid SEs",
  any(grepl("0.020$^{\\dagger}$", direct, fixed = TRUE) & grepl("--", direct, fixed = TRUE) & !grepl("0.001", direct, fixed = TRUE)))

# Regression checks for Phase 1.2 fixes: parallel worker fallback, do_parallel
# YAML propagation, and invalid simulation scenario handling without next-in-tryCatch.
source("helper/parallel_workers.R")
check("n.parallel.workers falls back to 1L when detectCores() is NA/non-finite",
  identical(n.parallel.workers(8, cores = NA_integer_, verbose = FALSE), 1L) &&
  identical(n.parallel.workers(8, cores = Inf, verbose = FALSE), 1L))
source("helper/config_loader.R")
sim.p.off <- scenario_to_estimation_params(list(id = "x", data_folder = "x",
  estimation = list(weighting = "EW", sdf_model = "linear", error_function = "L2_Lasso",
                    max_months = 0, lambdas = 0, kernel_bandwidth = "llsw", do_parallel = FALSE)))
emp.p.off <- empirical_scenario_to_params(list(id = "y", weighting = "EW", do_cross_validation = FALSE, do_parallel = FALSE))
emp.p.def <- empirical_scenario_to_params(list(id = "z", weighting = "EW", do_cross_validation = FALSE))
check("config_loader propagates do_parallel for simulation and empirical scenarios",
  identical(sim.p.off$do_parallel, FALSE) && identical(emp.p.off$do_parallel, FALSE) && identical(emp.p.def$do_parallel, TRUE))
bad.cfg <- list(defaults = list(dgp = list(), estimation = list()),
  scenarios = list(broken = list(active = TRUE, dgp = list(), estimation = list())))
bad.yaml <- file.path(stage, "config/bad_sim.yaml")
yaml::write_yaml(bad.cfg, bad.yaml)
bad.res <- tryCatch(invisible(capture.output(r.bad <- run_simulation_study(config_path = bad.yaml,
  generate_data = FALSE, estimate = TRUE, analyze = FALSE, verbose = FALSE))), error = function(e) e)
check("run_simulation_study skips invalid scenario cleanly without next-in-tryCatch error",
  !inherits(bad.res, "error") && isFALSE(r.bad$broken$success) && nzchar(r.bad$broken$error))

setwd(repo)
cat(n.fail, "review-fix regression(s) failed; artifacts:", stage, "\n")
if (n.fail > 0) quit(status = 1)

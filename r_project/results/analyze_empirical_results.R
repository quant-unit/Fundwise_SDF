# Analysis step between estimation and the empirical exhibits. Builds the summaries that
# the exhibit scripts read (0_asymptotic_inference_summary.csv, 0_cross_validation_summary.csv,
# 0_best_models_summary.csv, the *_sumabs.csv files and empirical_tables.tex) by running
# results/analyze_result.R on every cache folder of the empirical result folders, each in its
# own R process. Run from r_project/ after run_empirical_study.R:
#   LANG=en_US.UTF-8 Rscript results/analyze_empirical_results.R [result folder ...]
# Default result folders: results/data_out_2026-emp-c, -allfunds, -navtrim20.
# A folder without current cached rows (for example, results of a scenario not rerun after a
# change of the prepared data) is reported and skipped, and its old summaries stay in place;
# the exhibit scripts reject such summaries as stale. Any other error stops the run.
args <- commandArgs(TRUE)
folders <- if (length(args)) args else
  c("results/data_out_2026-emp-c", "results/data_out_2026-emp-c-allfunds", "results/data_out_2026-emp-c-navtrim20")
prefix <- "q_factors_preqin_"
rscript <- file.path(R.home("bin"), "Rscript")
skipped <- character(0); done <- 0L
for (out in folders) {
  caches <- list.dirs(out, recursive = FALSE)
  caches <- caches[startsWith(basename(caches), paste0("cache_", prefix))]
  for (d in caches) {
    if (!length(list.files(d, "_cached_res\\.csv$"))) next
    suffix <- sub(paste0("^cache_", prefix), "", basename(d))
    code <- sprintf("library(here); source.internally <- FALSE; data.out.folder <- %s; prefix <- %s; suffix <- %s; source(here('results/analyze_result.R'))",
      deparse(out), deparse(prefix), deparse(suffix))
    log <- suppressWarnings(system2(rscript, c("-e", shQuote(code)), stdout = TRUE, stderr = TRUE))
    status <- attr(log, "status")
    if (is.null(status) || status == 0) { done <- done + 1L; next }
    if (any(grepl("No current cached results", log, fixed = TRUE))) { skipped <- c(skipped, d); next }
    cat(log, sep = "\n")
    stop("analysis failed for ", d)
  }
}
message("analyzed ", done, " cache folder(s)")
if (length(skipped)) message("skipped (no current cached rows): ", paste(skipped, collapse = "; "))

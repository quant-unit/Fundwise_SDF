# Sample-robustness driver (19 November 2026; revision plan, log entries
# "robustness runs: scope" and "robustness runs"). Run from r_project/ with a
# clean, committed tree and LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/robustness_runs/run_robustness.R" <log file>
# Needs empirical/data_prepared_2026_c_navtrim20/ (empirical/prep_nav_threshold.R).
# Loads only the function definitions of run_empirical_study.R, as the G1
# driver does (archive/g1_runs/run_g1_empirical.R).
for (.e in parse("run_empirical_study.R")) {
  if (is.call(.e) && (identical(.e[[1]], as.name("source")) ||
      (identical(.e[[1]], as.name("<-")) && is.call(.e[[3]]) && identical(.e[[3]][[1]], as.name("function"))))) eval(.e, globalenv())
}
.ids <- as.vector(outer(c("ew_no_cv_north_america", "fw_no_cv_north_america", "u_ew_no_cv_north_america", "u_fw_no_cv_north_america"),
                        c("allfunds", "navtrim20"), paste, sep = "_"))
.log <- commandArgs(TRUE)[1]
cat(format(Sys.time()), "robustness start, commit", system("git rev-parse --short HEAD", intern = TRUE),
    if (length(system("git status --porcelain --untracked-files=no", intern = TRUE))) "WITH UNCOMMITTED CHANGES" else "(clean tree)", "\n", file = .log, append = TRUE)
for (.id in .ids) {
  .t0 <- Sys.time()
  .out <- capture.output(.r <- tryCatch(run_empirical_study(scenario_ids = .id, verbose = TRUE), error = function(e) e))
  .ok <- if (inherits(.r, "error")) paste("ERROR:", conditionMessage(.r)) else paste(vapply(.r, function(x) if (isTRUE(x$success)) "ok" else paste("FAILED:", x$error), ""), collapse = ";")
  cat(sprintf("%s %-40s %6.0f s  %s\n", format(Sys.time(), "%H:%M:%S"), .id, as.numeric(difftime(Sys.time(), .t0, units = "secs")), .ok), file = .log, append = TRUE)
}
cat(format(Sys.time()), "robustness done\n", file = .log, append = TRUE)

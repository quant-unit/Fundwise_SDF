# Archived Gate G1 driver (17 November 2026). Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_runs/run_g1_empirical.R" <log file>
# Gate G1 empirical driver: loads only the function definitions of
# run_empirical_study.R (its bottom block would run a hard-coded list)
for (.e in parse("run_empirical_study.R")) {
  if (is.call(.e) && (identical(.e[[1]], as.name("source")) ||
      (identical(.e[[1]], as.name("<-")) && is.call(.e[[3]]) && identical(.e[[3]][[1]], as.name("function"))))) eval(.e, globalenv())
}
.g1.ids <- c(
  # core: single- and multi-factor, no CV and CV, EW and FW; U_{v,H}
  "ew_no_cv_north_america", "fw_no_cv_north_america", "ew_cv_north_america", "fw_cv_north_america",
  "u_ew_no_cv_north_america", "u_ew_no_cv_north_america_qn", "u_fw_no_cv_north_america", "u_fw_no_cv_north_america_qn",
  # NAV haircuts and vintage cutoffs
  "fw_no_cv_north_america_NC0", "ew_no_cv_north_america_NC0", "fw_no_cv_north_america_NC25", "ew_no_cv_north_america_NC25",
  "fw_no_cv_north_america_NC50", "ew_no_cv_north_america_NC50", "fw_no_cv_north_america_NC75", "ew_no_cv_north_america_NC75",
  "fw_vintage_sweep_north_america", "ew_vintage_sweep_north_america", "fw_vintage_sweep_north_america_NC50", "ew_vintage_sweep_north_america_NC50",
  # Fama-French 3 and L1 robustness
  "ff3_fw_no_cv", "ff3_ew_no_cv", "ff3_fw_cv", "ff3_ew_cv", "ff3_fw_no_cv_alpha", "ff3_ew_no_cv_alpha", "ff3_fw_cv_alpha", "ff3_ew_cv_alpha",
  "ff3_fw_no_cv_ALL", "ff3_ew_no_cv_ALL", "ff3_fw_cv_ALL", "ff3_ew_cv_ALL", "ff3_fw_no_cv_ALL_3f", "ff3_ew_no_cv_ALL_3f", "ff3_fw_cv_ALL_3f", "ff3_ew_cv_ALL_3f",
  "fw_no_cv_north_america_l1", "ew_no_cv_north_america_l1", "fw_cv_north_america_l1", "ew_cv_north_america_l1")
.g1.log <- commandArgs(TRUE)[1]
cat(format(Sys.time()), "G1 start, commit", system("git rev-parse --short HEAD", intern = TRUE), "\n", file = .g1.log, append = TRUE)
for (.id in .g1.ids) {
  .t0 <- Sys.time()
  .out <- capture.output(.r <- tryCatch(run_empirical_study(scenario_ids = .id, verbose = TRUE), error = function(e) e))
  .ok <- if (inherits(.r, "error")) paste("ERROR:", conditionMessage(.r)) else paste(vapply(.r, function(x) if (isTRUE(x$success)) "ok" else paste("FAILED:", x$error), ""), collapse = ";")
  cat(sprintf("%s %-38s %6.0f s  %s\n", format(Sys.time(), "%H:%M:%S"), .id, as.numeric(difftime(Sys.time(), .t0, units = "secs")), .ok), file = .g1.log, append = TRUE)
}
cat(format(Sys.time()), "G1 done\n", file = .g1.log, append = TRUE)

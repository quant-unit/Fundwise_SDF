# Archived Gate G1 driver (17 November 2026). Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_runs/run_g1_simulation_retry.R" <log file>
# Gate G1: rerun the four failed Appendix D scenarios (data already
# generated at 1eb3628; estimation resumes from cached horizons), retrying
# up to three times, then the bias analysis over all 20 scenarios.
source("run_simulation_study.R")
.base <- c("base_case_zero_alpha", "big_n_v_50funds_alpha", "big_n_v_50funds_alpha_stdv30", "big_n_v_50funds_alpha_stdv30_shifted",
  "big_n_v_50funds_alpha_stdv30_shifted_mkt", "big_n_v_50funds_alpha_stdv30_shifted_mkt2", "big_n_v_50funds_alpha_stdv30_shifted_mkt4",
  "big_n_v_50funds_alpha_stdv30_shifted_mkt6", "big_n_v_50funds_alpha_stdv30_shifted_mkt7", "big_n_v_50funds_alpha_stdv30_shifted_mkt8")
.all <- c(.base, paste0(.base, "_additive"))
.redo <- c("big_n_v_50funds_alpha_stdv30_shifted_mkt2", "big_n_v_50funds_alpha_stdv30_shifted_mkt2_additive",
  "base_case_zero_alpha", "base_case_zero_alpha_additive")
.log <- commandArgs(TRUE)[1]; .out <- "simulation/data_out_2026_g1"
cat(format(Sys.time()), "G1 simulation retry, commit", system("git rev-parse --short HEAD", intern = TRUE), "\n", file = .log, append = TRUE)
for (.id in .redo) for (.try in 1:3) {
  .t0 <- Sys.time()
  .r <- tryCatch(run_simulation_study(scenario_ids = .id, generate_data = FALSE, estimate = TRUE,
    analyze = FALSE, data_out_folder = .out), error = function(e) e)
  .ok <- if (inherits(.r, "error")) paste("ERROR:", conditionMessage(.r)) else paste(vapply(.r, function(x) if (isTRUE(x$success)) "ok" else paste("FAILED:", x$error), ""), collapse = ";")
  cat(sprintf("%s %-50s try %d %6.0f s  %s\n", format(Sys.time(), "%H:%M:%S"), .id, .try, as.numeric(difftime(Sys.time(), .t0, units = "secs")), .ok), file = .log, append = TRUE)
  if (.ok == "ok") break
}
.r <- tryCatch(run_simulation_study(scenario_ids = .all, generate_data = FALSE, estimate = FALSE, analyze = TRUE, data_out_folder = .out), error = function(e) e)
cat(format(Sys.time()), "analysis:", if (inherits(.r, "error")) conditionMessage(.r) else "ok", "\n", file = .log, append = TRUE)
cat(format(Sys.time()), "G1 simulation retry done\n", file = .log, append = TRUE)

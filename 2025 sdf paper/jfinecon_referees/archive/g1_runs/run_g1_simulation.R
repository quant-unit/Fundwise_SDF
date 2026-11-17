# Archived Gate G1 driver (17 November 2026). Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_runs/run_g1_simulation.R" <log file>
# Gate G1 simulation driver: Appendix D (DLP12 replication) S1-S10 with the
# B1 multiplicative shocks and their additive-shock twins.
source("run_simulation_study.R")
.base <- c("base_case_zero_alpha", "big_n_v_50funds_alpha", "big_n_v_50funds_alpha_stdv30", "big_n_v_50funds_alpha_stdv30_shifted",
  "big_n_v_50funds_alpha_stdv30_shifted_mkt", "big_n_v_50funds_alpha_stdv30_shifted_mkt2", "big_n_v_50funds_alpha_stdv30_shifted_mkt4",
  "big_n_v_50funds_alpha_stdv30_shifted_mkt6", "big_n_v_50funds_alpha_stdv30_shifted_mkt7", "big_n_v_50funds_alpha_stdv30_shifted_mkt8")
.ids <- c(.base, paste0(.base, "_additive"))
.log <- commandArgs(TRUE)[1]; .out <- "simulation/data_out_2026_g1"
cat(format(Sys.time()), "G1 simulation start, commit", system("git rev-parse --short HEAD", intern = TRUE), "\n", file = .log, append = TRUE)
for (.id in .ids) {
  .t0 <- Sys.time()
  .o <- capture.output(.r <- tryCatch(run_simulation_study(scenario_ids = .id, generate_data = TRUE, estimate = TRUE,
    analyze = FALSE, data_out_folder = .out), error = function(e) e))
  .ok <- if (inherits(.r, "error")) paste("ERROR:", conditionMessage(.r)) else paste(vapply(.r, function(x) if (isTRUE(x$success)) "ok" else paste("FAILED:", x$error), ""), collapse = ";")
  cat(sprintf("%s %-48s %6.0f s  %s\n", format(Sys.time(), "%H:%M:%S"), .id, as.numeric(difftime(Sys.time(), .t0, units = "secs")), .ok), file = .log, append = TRUE)
}
.o <- capture.output(.r <- tryCatch(run_simulation_study(scenario_ids = .ids, generate_data = FALSE, estimate = FALSE, analyze = TRUE, data_out_folder = .out), error = function(e) e))
cat(format(Sys.time()), "analysis:", if (inherits(.r, "error")) conditionMessage(.r) else "ok", "\n", file = .log, append = TRUE)
cat(format(Sys.time()), "G1 simulation done\n", file = .log, append = TRUE)

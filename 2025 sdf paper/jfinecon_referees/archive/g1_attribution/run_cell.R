# Attribution cell: run this checkout's estimator on one data folder
# args: data_prepared_folder, output folder, log file
.a <- commandArgs(TRUE); .data.folder <- .a[1]; .outdir <- .a[2]; .log <- .a[3]; .ids <- if (length(.a) > 3) .a[-(1:3)] else c("ew_no_cv_north_america", "fw_no_cv_north_america")
for (.e in parse("run_empirical_study.R")) {
  if (is.call(.e) && (identical(.e[[1]], as.name("source")) ||
      (identical(.e[[1]], as.name("<-")) && is.call(.e[[3]]) && identical(.e[[3]][[1]], as.name("function"))))) eval(.e, globalenv())
}
.cfg <- readLines("config/empirical_scenarios.yaml")
.cfg <- sub('^(  data_prepared_folder: ).*$', paste0('\\1"', .data.folder, '"'), .cfg)
.cfg.path <- file.path(.outdir, "config.yaml"); dir.create(.outdir, showWarnings = FALSE, recursive = TRUE); writeLines(.cfg, .cfg.path)
.n.fail <- 0
for (.id in .ids) {
  .t0 <- Sys.time()
  .o <- capture.output(.r <- tryCatch(run_empirical_study(scenario_ids = .id, config_path = .cfg.path, data_out_folder = .outdir, verbose = TRUE), error = function(e) e))
  writeLines(.o, file.path(.outdir, paste0(.id, "_stdout.log")))
  .ok <- if (inherits(.r, "error")) paste("ERROR:", conditionMessage(.r)) else paste(vapply(.r, function(x) if (isTRUE(x$success)) "ok" else paste("FAILED:", x$error), ""), collapse = ";")
  cat(sprintf("%s %s %s %s %.0f s %s\n", format(Sys.time(), "%H:%M:%S"), basename(getwd()), .data.folder, .id, as.numeric(difftime(Sys.time(), .t0, units = "secs")), .ok), file = .log, append = TRUE)
  if (.ok != "ok") .n.fail <- .n.fail + 1
}
if (.n.fail > 0) quit(status = 1) # stop reproduce.sh (set -e)

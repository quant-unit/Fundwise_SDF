# Archived Gate G1 driver (17 November 2026). Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_runs/run_g1_analysis.R" <log file>
# Gate G1 analysis driver (run from r_project/): 1) summarize every cache
# folder of the G1 empirical run with results/analyze_result.R; 2) rerun the
# North America figure/table block of results/plot_empirical_estimates.R on
# the new results, writing to results/figures_bovc_us_g1.
suppressMessages(library(here))
.log <- commandArgs(TRUE)[1]
.lg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", file = .log, append = TRUE)
.root <- "results/data_out_2026-emp-c"
.dirs <- sub("^cache_", "", list.files(here(.root), pattern = "^cache_"))
for (.d in .dirs) {
  .pre <- if (startsWith(.d, "ff3_factors_preqin_")) "ff3_factors_preqin_" else "q_factors_preqin_"
  .env <- new.env(parent = globalenv())
  assign("source.internally", FALSE, globalenv()); assign("prefix", .pre, .env)
  assign("suffix", substring(.d, nchar(.pre) + 1), .env); assign("data.out.folder", .root, .env)
  .msg <- tryCatch({ invisible(capture.output(suppressWarnings(sys.source(here("results/analyze_result.R"), envir = .env)))); "ok" },
                   error = function(e) paste("ERROR:", conditionMessage(e)))
  .lg("analyze", .d, .msg)
}
# 2) plots and tables
setwd(here("results"))
.ex <- parse("plot_empirical_estimates.R", keep.source = TRUE)
.src <- vapply(attr(.ex, "srcref"), function(r) paste(as.character(r), collapse = "\n"), "")
.start <- which(grepl('out.folder <- "figures_bovc_us"', .src, fixed = TRUE))[1]
for (.i in seq_along(.ex)) {
  .e <- .ex[[.i]]
  is.def <- is.call(.e) && identical(.e[[1]], as.name("<-")) && is.call(.e[[3]]) && identical(.e[[3]][[1]], as.name("function"))
  is.lib <- is.call(.e) && as.character(.e[[1]])[1] %in% c("library", "source")
  if (is.def || is.lib) eval(.e, globalenv())
}
file.folder <- "data_out_2026-emp-c"; out.folder <- "figures_bovc_us_g1"
dir.create(out.folder, showWarnings = FALSE)
for (.i in seq(.start + 1, length(.ex))) {
  .e <- .ex[[.i]]
  if (is.call(.e) && identical(.e[[1]], as.name("<-")) && as.character(.e[[2]]) %in% c("file.folder", "out.folder")) next
  if (!is.call(.e) || !grepl("^plot_", as.character(.e[[1]])[1])) next
  .of <- tryCatch(eval(.e$output_file, globalenv()), error = function(e) "?")
  .msg <- tryCatch({ invisible(capture.output(suppressWarnings(eval(.e, globalenv())))); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
  .lg("plot", basename(.of), .msg)
}
.lg("analysis done")

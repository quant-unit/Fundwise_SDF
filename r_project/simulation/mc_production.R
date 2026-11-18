# Step 4 Monte Carlo: production runs (revision plan, Section I, "Step 4:
# scope and compute plan", order of work). Run from r_project/ with a clean,
# committed tree (every batch records commit and uncommitted changes):
#   LANG=en_US.UTF-8 Rscript simulation/mc_production.R [design ...]
# Without arguments: grid G, then the two cores, then the margins M1-M7, each
# at the replications of config/mc_designs.yaml. Finished batches are skipped,
# so the script can be restarted at any time. Progress: one line per design
# in simulation/data_out_mc/production.log.
source("simulation/mc_driver.R")
pub <- mc.public.factors()
D <- mc.load.designs()
ids <- commandArgs(TRUE)
if (!length(ids)) ids <- c(grep("^g_", names(D), value = TRUE), "core_BO", "core_VC", grep("^m[0-9]", names(D), value = TRUE))
stopifnot(all(ids %in% names(D)))
log.file <- "simulation/data_out_mc/production.log"
dir.create(dirname(log.file), showWarnings = FALSE, recursive = TRUE)
lg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n", file = log.file, append = TRUE)
dirty <- length(system("git status --porcelain --untracked-files=no", intern = TRUE)) > 0
lg("start", length(ids), "designs at commit", system("git rev-parse --short HEAD", intern = TRUE), if (dirty) "WITH UNCOMMITTED CHANGES" else "(clean tree)")
for (id in ids) {
  t0 <- Sys.time()
  rows <- tryCatch(mc.run(id, D[[id]], seq_len(D[[id]]$reps), pub, verbose = FALSE), error = function(e) { lg(id, "ERROR:", conditionMessage(e)); NULL })
  if (!is.null(rows)) lg(sprintf("%-14s %5d replications, %6d rows, %5.0f s", id, length(unique(rows$rep)), nrow(rows), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
lg("done")

# Step 4 Monte Carlo: criterion targets of every production design
# (mc.targets() in simulation/mc_analyze.R), on replications 1-200; for the
# two cores also on 1-100 and 101-200, to measure the targets' Monte Carlo
# precision. Run from r_project/ (writes simulation/data_out_mc/targets.rds):
#   LANG=en_US.UTF-8 Rscript simulation/mc_targets.R
source("simulation/mc_driver.R"); source("simulation/mc_analyze.R")
pub <- mc.public.factors(); D <- mc.load.designs()
w <- n.parallel.workers(48)
out <- list()
for (id in names(D)) {
  t0 <- Sys.time()
  tg <- mc.targets(D[[id]], pub, reps = 1:200, workers = w)
  tg$design <- id; tg$half <- "all"
  if (id %in% c("core_BO", "core_VC")) for (h in list(a = 1:100, b = 101:200)) {
    x <- mc.targets(D[[id]], pub, reps = h, workers = w); x$design <- id; x$half <- if (h[1] == 1) "1-100" else "101-200"
    tg <- rbind(tg, x)
  }
  out[[id]] <- tg
  cat(format(Sys.time(), "%H:%M:%S"), id, nrow(tg), "cells,", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}
targets <- do.call(rbind, out); rownames(targets) <- NULL
saveRDS(list(targets = targets, commit = system("git rev-parse --short HEAD", intern = TRUE),
             dirty = length(system("git status --porcelain --untracked-files=no", intern = TRUE)) > 0), "simulation/data_out_mc/targets.rds")
cat("done\n")

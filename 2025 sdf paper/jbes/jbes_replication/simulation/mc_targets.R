# Step 4 Monte Carlo: criterion targets of every production design
# (mc.targets() in simulation/mc_analyze.R), on replications 1-200; for the
# two cores also on 1-100 and 101-200, to measure the targets' Monte Carlo
# precision. Run from r_project/:
#   LANG=en_US.UTF-8 Rscript simulation/mc_targets.R              -> targets.rds (every design)
#   LANG=en_US.UTF-8 Rscript simulation/mc_targets.R design ...   -> targets_additions.rds (those designs only;
#     the additions of 19 November, so the targets of the first production stay as committed)
#   LANG=en_US.UTF-8 Rscript simulation/mc_targets.R --out=targets_NAME.rds design ...
#     -> that file (later additions, so earlier target files are never overwritten)
# All in simulation/data_out_mc/; simulation/mc_report.R reads every targets*.rds there.
source("simulation/mc_driver.R"); source("simulation/mc_analyze.R")
pub <- mc.public.factors(); D <- mc.load.designs()
args <- commandArgs(TRUE)
out.arg <- grep("^--out=", args, value = TRUE)
ids <- setdiff(args, out.arg)
stopifnot(all(ids %in% names(D)), length(out.arg) <= 1, !length(out.arg) || length(ids))
out.name <- if (length(out.arg)) sub("^--out=", "", out.arg) else if (length(ids)) "targets_additions.rds" else "targets.rds"
stopifnot(grepl("^targets.*\\.rds$", out.name))
out.file <- file.path("simulation/data_out_mc", out.name)
if (length(out.arg) && file.exists(out.file)) stop(out.file, " exists; choose another --out name")
if (!length(ids)) ids <- names(D)
w <- n.parallel.workers(48)
out <- list()
for (id in ids) {
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
spec.md5 <- setNames(vapply(ids, function(id) mc.spec.md5(D[[id]]), ""), ids)
saveRDS(list(targets = targets, spec.md5 = spec.md5, public.md5 = value.md5(pub), settings = MC.SETTINGS,
             commit = system("git rev-parse --short HEAD", intern = TRUE),
             dirty = length(system("git status --porcelain --untracked-files=no", intern = TRUE)) > 0), out.file)
cat("done\n")

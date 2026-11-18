# One-time migration (18 November 2026). Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/housekeeping_2026-11-18/migrate_cv_markers.R" [result folder]
# cache.identity no longer hashes the configuration file and uses relative
# data paths (helper/result_cache.R). CV completion markers are named by the
# identity, so for every current result row's former identity this copies
# its markers to the name under the new identity (old markers are kept).
# Without it, resuming a CV scenario would refit its folds.
suppressMessages(source("helper/result_cache.R"))
root <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "results/data_out_2026-emp-c"
n.new <- 0L; n.have <- 0L
for (d in list.dirs(root, recursive = FALSE)) {
  files <- list.files(d, "_cached_res\\.csv$", full.names = TRUE)
  if (!length(files)) next
  raw <- do.call(rbind, lapply(files, function(f) {
    x <- read.result.csv(f)
    if (!all(RESULT.IDENTITY.COLUMNS %in% names(x))) return(NULL)
    x[, RESULT.IDENTITY.COLUMNS, drop = FALSE]
  }))
  if (is.null(raw)) next
  raw <- raw[result.identities.current(raw), , drop = FALSE] # data files unchanged
  if (!nrow(raw)) next
  pairs <- unique(data.frame(old = raw$cache.identity, new = result.identity.hash(raw)))
  for (i in seq_len(nrow(pairs))) {
    old <- list.files(d, paste0("^cv_complete_", pairs$old[i], "_"), full.names = TRUE)
    target <- file.path(d, sub(pairs$old[i], pairs$new[i], basename(old), fixed = TRUE))
    n.have <- n.have + sum(file.exists(target))
    n.new <- n.new + sum(file.copy(old[!file.exists(target)], target[!file.exists(target)]))
  }
}
cat("markers copied:", n.new, "| already present:", n.have, "\n")

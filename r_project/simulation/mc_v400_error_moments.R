# Moments of the Monte Carlo estimation error beta-hat - 1 for equal-vintage U_{v,H}
# (q_v = 1, single-factor MKT) in the growing-V designs on resampled factor paths.
# Source for the Online Supplement statement that beta-hat - 1 is already close to
# Gaussian at V = 400 (Appendix D, studentization diagnosis).
# Reads the stored production batches (simulation/data_out_mc/<design>/batch_*.rds);
# no estimation. Run from r_project/:  LANG=en_US.UTF-8 Rscript simulation/mc_v400_error_moments.R
# Writes results/mc_step4/v400_error_moments.csv.

designs <- c("v_V25_n5", "v_V100_n5", "v_V400_n5", "v_V25_n20", "v_V100_n20", "v_V400_n20")
rows <- list()
for (d in designs) {
  f <- list.files(file.path("simulation/data_out_mc", d), pattern = "^batch_.*\\.rds$", full.names = TRUE)
  if (!length(f)) stop("no stored batches for design ", d)
  b <- lapply(f, readRDS)
  r <- do.call(rbind, lapply(b, `[[`, "rows"))
  commits <- unique(vapply(b, function(x) x$provenance$git.commit, ""))
  dirty <- any(vapply(b, function(x) isTRUE(x$provenance$git.dirty), TRUE))
  for (m in c(0, 180)) {
    x <- r[r$estimator == "U.eq" & r$model == "MKT" & r$max.month == m & is.finite(r$MKT), ]
    stopifnot(!anyDuplicated(x$rep))
    e <- x$MKT - 1
    ec <- e - mean(e)
    rows[[length(rows) + 1L]] <- data.frame(
      design = d, H = m / 12, reps = nrow(x), mean_err = mean(e), sd_err = sd(e),
      sd_times_sqrtV = sd(e) * sqrt(as.integer(sub("^v_V([0-9]+)_n[0-9]+$", "\\1", d))),
      kurtosis_err = mean(ec^4) / mean(ec^2)^2, skewness_err = mean(ec^3) / mean(ec^2)^1.5,
      batch_commits = paste(commits, collapse = ";"), batch_dirty = dirty, stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
write.csv(out, "results/mc_step4/v400_error_moments.csv", row.names = FALSE)
print(out, row.names = FALSE, digits = 4)
cat("wrote results/mc_step4/v400_error_moments.csv\n")

# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/compare_r3_r4.R" <out.rds>
# G1 rerun with the new inference (r4 rows: kernel.bandwidth "llsw") vs r3
# (bandwidth 12): point estimates and independent-units SEs must be identical.
suppressMessages(source("helper/result_cache.R")); invisible(capture.output(suppressMessages(source("helper/getNPVs.R"))))
root <- "results/data_out_2026-emp-c"; out <- list(); cmp <- list()
par.cols <- c("MKT", "ME", "IA", "ROE", "EG", "Alpha", "SMB", "HML")
for (d in list.dirs(root, recursive = FALSE)) {
  raw <- Reduce(rbind.all.columns, lapply(Sys.glob(file.path(d, "*_cached_res.csv")), read.csv))
  r3 <- dedupe.result.rows(raw[raw$pipeline.version %in% "2026-11-16-r3" & as.character(raw$kernel.bandwidth) == "12", ])
  r4 <- read.result.cache(d, verbose = FALSE)
  k <- c("Type", "Factor", "max.month", "CV.key", "lambda")
  m <- merge(r4, r3, by = k, suffixes = c(".n", ".o"))
  pc <- intersect(par.cols, sub("\\.n$", "", grep("\\.n$", names(m), value = TRUE)))
  dpar <- if (nrow(m)) max(sapply(pc, function(v) max(abs(m[[paste0(v, ".n")]] - m[[paste0(v, ".o")]]), na.rm = TRUE)), na.rm = TRUE) else NA
  dind <- if (nrow(m) && "SE.MKT.indep.n" %in% names(m)) max(abs(m$SE.MKT.indep.n - m$SE.MKT.indep.o), na.rm = TRUE) else NA
  out[[basename(d)]] <- data.frame(folder = basename(d), r3 = nrow(r3), r4 = nrow(r4), matched = nrow(m), max.dpar = dpar, max.dindep = dind,
    S = paste(sort(unique(r4$hac.S)), collapse = ","), jk.max = suppressWarnings(max(r4$jk.selfcheck, na.rm = TRUE)))
  if (nrow(m)) cmp[[basename(d)]] <- data.frame(folder = basename(d), m[, c(k, "MKT.n", "SE.MKT.o", "SE.MKT.n", "Wald.p.value.MKT_1.o", "Wald.p.value.MKT_1.n")])
}
o <- do.call(rbind, out); rownames(o) <- NULL
cat(sprintf("folders %d | all r4 rows matched to r3: %s | max |point estimate change| %.1e | max |indep SE change| %.1e\n", nrow(o),
  all(o$matched == o$r4 & o$r4 > 0), max(o$max.dpar, na.rm = TRUE), max(o$max.dindep, na.rm = TRUE)))
cat("rows r3 / r4 total:", sum(o$r3), "/", sum(o$r4), "| folders with unmatched rows:", sum(o$matched != o$r4), "\n")
cat("bandwidths S seen:", paste(sort(unique(unlist(strsplit(o$S, ",")))), collapse = " "), "| jackknife self-check max:", signif(max(o$jk.max[is.finite(o$jk.max)]), 3), "\n")
print(o[o$matched != o$r4 | o$max.dpar > 0, ], row.names = FALSE)
c <- do.call(rbind, cmp); c <- c[c$CV.key == "ALL", ]
c$ratio <- c$SE.MKT.n / c$SE.MKT.o
cat(sprintf("\nfull-sample rows %d: SE ratio new/old median %.2f, range %.2f-%.2f\n", nrow(c), median(c$ratio, na.rm = TRUE), min(c$ratio, na.rm = TRUE), max(c$ratio, na.rm = TRUE)))
rej.o <- c$Wald.p.value.MKT_1.o < 0.05; rej.n <- c$Wald.p.value.MKT_1.n < 0.05
cat(sprintf("Wald test (beta = 1, others 0) at 5%%: rejected r3 %d, r4 %d of %d; no longer rejected %d, newly rejected %d\n",
  sum(rej.o, na.rm = TRUE), sum(rej.n, na.rm = TRUE), sum(!is.na(rej.o) & !is.na(rej.n)), sum(rej.o & !rej.n, na.rm = TRUE), sum(!rej.o & rej.n, na.rm = TRUE)))
saveRDS(list(o = o, c = c), commandArgs(TRUE)[1])

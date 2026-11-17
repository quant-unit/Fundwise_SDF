# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/table_diff.R" 
# r3 tables (commit 90e5421) vs r4 tables (working tree): SEs and 5% significance of full-sample
# estimates. Each data CSV holds all factor rows of its type and criterion; H = 1 month
# rows duplicate H = 0 and are dropped.
f <- c(BO.L2 = "empirical_BO_US_MKT_data.csv", VC.L2 = "empirical_VC_US_MKT_data.csv", BO.L1 = "empirical_BO_US_MKT_L1_data.csv", VC.L1 = "empirical_VC_US_MKT_L1_data.csv")
k <- c("Type", "Factor", "max.month", "weighting")
rows <- list()
for (nm in names(f)) {
  p <- file.path("results/figures_bovc_us_g1", f[[nm]])
  o <- read.csv(text = system(sprintf("git show 90e5421:r_project/%s", p), intern = TRUE)); n <- read.csv(p)
  stopifnot(!"p.MKT" %in% names(o), "p.MKT" %in% names(n))
  o <- o[o$method == "Asymptotic" & o$max.month != 1, ]; n <- n[n$method == "Asymptotic" & n$max.month != 1, ]
  m <- merge(o, n, by = k, suffixes = c(".o", ".n")); stopifnot(nrow(m) == nrow(n), max(abs(m$MKT.o - m$MKT.n)) == 0)
  for (par in c("MKT", "Coef")) {
    est <- m[[paste0(par, ".n")]]; so <- m[[paste0("SE.", par, ".o")]]; sn <- m[[paste0("SE.", par, ".n")]]; pn <- m[[paste0("p.", par)]]
    ok <- is.finite(est) & is.finite(so) & is.finite(sn)
    r3 <- abs(est / so) > qnorm(0.975); r4 <- pn < 0.05
    rows[[length(rows) + 1]] <- data.frame(set = nm, par = par, n = sum(ok), se.ratio.med = median((sn / so)[ok]), se.ratio.min = min((sn / so)[ok]), se.ratio.max = max((sn / so)[ok]),
      sig.r3 = sum(r3[ok]), sig.r4 = sum(r4[ok], na.rm = TRUE), lost = sum((r3 & !r4)[ok], na.rm = TRUE), gained = sum((!r3 & r4)[ok], na.rm = TRUE), p.na = sum(is.na(pn[ok])))
  }
  # MKT row of the single-factor model, by weighting and horizon (for the log)
  s <- m[m$Factor == "MKT", ]; s <- s[order(s$weighting, s$max.month), ]
  cat("\n", nm, "single-factor beta: weighting H SE.r3 SE.r4 t.r4 p.fixedb\n")
  print(data.frame(w = s$weighting, H = s$max.month / 12, SE.r3 = round(s$SE.MKT.o, 3), SE.r4 = round(s$SE.MKT.n, 3), t = round(s$MKT.n / s$SE.MKT.n, 2), p = signif(s$p.MKT, 2)), row.names = FALSE)
}
out <- do.call(rbind, rows); num <- sapply(out, is.numeric); out[num] <- lapply(out[num], round, 2); print(out, row.names = FALSE)

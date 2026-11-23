# Step 4 Monte Carlo: two inference checks quoted in the paper that are not in
# summary.csv (audit of 23 November 2026; the numbers were typed into the
# exhibit notes without a script).
# Run from r_project/:  LANG=en_US.UTF-8 Rscript simulation/mc_inference_checks.R
#
# 1. ar_replay: the centered Bartlett fixed-b moment test of the paper's
#    identification-robust IV sets (results/empirical_iv.R, iv.ar.set()),
#    evaluated at the true beta = 1, H = 0, on the regenerated draws of the
#    core, persistent-shock (m4p) and deal-level path-shock (pd_, design A')
#    designs. Draws are regenerated from the seeds (mc.draw(design, rep)); the
#    script first checks that the regenerated draw of replication 1 reproduces
#    the stored IV estimate of iv_core_BO exactly. Also reports the uncentered
#    variant (scores not demeaned) for comparison.
# 2. v400_tails: for the resampled-path V = 400 designs, U_{v,H} (q_v = 1)
#    single-factor Wald test of beta = 1 from the stored batches: rejection
#    rate, left/right rejections, median N_eff, and left/right rejections when
#    the sandwich SE is replaced by the Monte Carlo SD of the estimates.
# Writes results/mc_step4/inference_checks.csv (simulated results only).
suppressMessages({ source("simulation/mc_driver.R"); if (!exists("fixedb.crit.t")) source("helper/fixedb.R") })
D <- mc.load.designs(); pub <- mc.public.factors()
GIT <- sprintf("%s%s", system("git rev-parse --short HEAD", intern = TRUE),
  if (length(system("git status --porcelain --untracked-files=no", intern = TRUE))) " (uncommitted changes)" else "")

# ---- 1. AR moment test at beta = 1, H = 0 (as iv.ar.set() in results/empirical_iv.R) ----
ar.t <- function(fp, center = TRUE) {
  e <- fast.unit.errors(c(MKT = 1), fp, 0); V <- fp$n.units; S <- hac.bandwidth.llsw(V)
  K <- outer(fp$vintage, fp$vintage, bartlett.weight, bw = S); g <- mean(e); s <- if (center) e - g else e
  c(t = abs(g) / sqrt(sum(K * outer(s, s)) / V / V), crit = fixedb.crit.t(0.05, S / V))
}
r1 <- mc.read("iv_core_BO"); r1 <- r1$MKT[r1$rep == 1 & r1$estimator == "IV" & r1$max.month == 0]
fp1 <- mc.as.fast(mc.draw(D[["iv_core_BO"]]$design, 1, pub), pub, "vyp")
chk <- abs(r1 - mc.fit.iv(MC.ESTIMATORS$IV, fp1, "MKT", 0)$MKT)
cat("check: regenerated draw reproduces the stored IV estimate of iv_core_BO rep 1:", signif(chk, 3), "\n")
stopifnot(chk < 1e-10)
workers <- if (exists("n.parallel.workers")) n.parallel.workers(1000, verbose = FALSE) else 1L
rows <- list()
for (id in c("core_BO", "core_VC", "m4p_shared", "m4p_shared2", "pd_V25_n20_c", "pd_V25_n20_2c", "pd_V100_n20_c", "pd_V100_n20_2c")) {
  n <- D[[id]]$reps
  m <- do.call(rbind, parallel::mclapply(seq_len(n), function(k) {
    fp <- mc.as.fast(mc.draw(D[[id]]$design, k, pub), pub, "vyp")
    c(ar.t(fp, TRUE), t.unc = ar.t(fp, FALSE)[["t"]])
  }, mc.cores = workers))
  rows[[length(rows) + 1]] <- data.frame(check = "ar_replay", design = id, H = 0, reps = n, crit = m[1, "crit"],
    reject = mean(m[, "t"] > m[, "crit"]), reject.uncentered = mean(m[, "t.unc"] > m[, "crit"]),
    left = NA, right = NA, median.neff = NA, left.mcsd = NA, right.mcsd = NA, stringsAsFactors = FALSE)
  cat(sprintf("ar_replay %-16s reps %4d reject %5.1f%% (uncentered %4.1f%%)\n", id, n, 100 * tail(rows, 1)[[1]]$reject, 100 * tail(rows, 1)[[1]]$reject.uncentered))
}

# ---- 2. V = 400 tails of the U_{v,H} Wald test (stored batches) ----
for (id in c("v_V400_n5", "v_V400_n20")) {
  r <- mc.read(id)
  for (m in c(0, 180)) {
    x <- r[r$estimator == "U.eq" & r$model == "MKT" & r$max.month == m & is.finite(r$SE.MKT), ]
    cr <- fixedb.crit.t(0.05, unique(x$hac.b)); t <- (x$MKT - 1) / x$SE.MKT; t2 <- (x$MKT - 1) / sd(x$MKT)
    rows[[length(rows) + 1]] <- data.frame(check = "v400_tails", design = id, H = m / 12, reps = nrow(x), crit = cr,
      reject = mean(abs(t) > cr), reject.uncentered = NA, left = sum(t < -cr), right = sum(t > cr),
      median.neff = median(x$neff.MKT), left.mcsd = sum(t2 < -cr), right.mcsd = sum(t2 > cr), stringsAsFactors = FALSE)
    cat(sprintf("v400_tails %-10s H%2d reject %5.1f%% left/right %d/%d median neff %.1f | MC-SD left/right %d/%d\n",
      id, m / 12, 100 * mean(abs(t) > cr), sum(t < -cr), sum(t > cr), median(x$neff.MKT), sum(t2 < -cr), sum(t2 > cr)))
  }
}
out <- do.call(rbind, rows); out$git <- GIT
write.csv(out, "results/mc_step4/inference_checks.csv", row.names = FALSE)
cat("wrote results/mc_step4/inference_checks.csv\n")

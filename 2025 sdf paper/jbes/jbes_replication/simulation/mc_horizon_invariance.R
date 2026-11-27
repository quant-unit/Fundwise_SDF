# Size and power of the specification tests of helper/horizon_invariance.R on Step 4
# Monte Carlo designs (common random numbers: replication r is mc.draw(design, r), as in
# simulation/mc_production.R, so the U estimates equal the production ones).
#   null (T7 holds, c_v = 0):           core_BO, core_VC (EW, U.eq); s_core_BO, s_core_VC (FW, U.eq.fw)
#   shared vintage shocks (c_v > 0):    m4_shared (c), m4_shared2 (2c), m4p_shared2 (2c, AR(1) 0.6),
#                                       u_vc_shared (VC shape, sigma_W = 0.50)
#   factor-path dependence ((X) fails): pc_V25_n20_2c (Design A')
# Horizon-invariance test at H = 15 yr (primary) and 12.5 yr; two-moment test at H = 0.
# Run from r_project/:  LANG=en_US.UTF-8 Rscript simulation/mc_horizon_invariance.R
# Writes simulation/data_out_mc/horizon_invariance/<design>.rds (replication rows, gitignored)
# and results/mc_step4/horizon_invariance_mc.csv (summary).
suppressMessages({
  source("simulation/mc_driver.R")
  if (!exists("fixedb.crit.t")) source("helper/fixedb.R")
  source("helper/horizon_invariance.R")
  source("helper/parallel_workers.R")
})
pub <- mc.public.factors()
D <- mc.load.designs()
RUNS <- list(
  list(d = "core_BO", est = "U.eq", reps = 500), list(d = "core_VC", est = "U.eq", reps = 500),
  list(d = "s_core_BO", est = "U.eq.fw", reps = 500), list(d = "s_core_VC", est = "U.eq.fw", reps = 500),
  list(d = "m4_shared", est = "U.eq", reps = 250), list(d = "m4_shared2", est = "U.eq", reps = 250),
  list(d = "m4p_shared2", est = "U.eq", reps = 250), list(d = "u_vc_shared", est = "U.eq", reps = 250),
  list(d = "pc_V25_n20_2c", est = "U.eq", reps = 250))
H.TEST <- c(150, 180)
out.dir <- "simulation/data_out_mc/horizon_invariance"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

one.rep <- function(r, design, est.name) {
  est <- MC.ESTIMATORS[[est.name]]
  dr <- mc.draw(design, r, pub)
  fp <- mc.as.fast(dr, pub, est$level)
  f0 <- mc.fit(est, fp, "MKT", 0)
  hi <- do.call(rbind, lapply(H.TEST, function(m) u.horizon.contrast(est, fp, 0, m, f0 = f0)))
  tm <- u.two.moment.test(est, fp)
  data.frame(rep = r, hi, tm.p.max = tm$p.max, tm.set.share = tm$set.share, stringsAsFactors = FALSE)
}

summ <- list()
for (run in RUNS) {
  f <- file.path(out.dir, paste0(run$d, ".rds"))
  if (file.exists(f) && identical(readRDS(f)$reps, seq_len(run$reps))) {
    rows <- readRDS(f)$rows
  } else {
    t0 <- proc.time()[["elapsed"]]
    res <- parallel::mclapply(seq_len(run$reps), one.rep, design = D[[run$d]]$design, est.name = run$est,
      mc.cores = n.parallel.workers(run$reps), mc.preschedule = FALSE)
    bad <- !vapply(res, is.data.frame, TRUE)
    if (any(bad)) stop(run$d, ": ", sum(bad), " replications failed, e.g. ", as.character(res[[which(bad)[1]]])[1])
    rows <- do.call(rbind, res)
    saveRDS(list(design = run$d, estimator = run$est, reps = seq_len(run$reps), rows = rows), f)
    cat(sprintf("%s: %d reps in %.0f s\n", run$d, run$reps, proc.time()[["elapsed"]] - t0))
  }
  for (m in H.TEST) {
    x <- rows[rows$m1 == m, ]
    ok <- is.finite(x$p)
    tm <- rows[rows$m1 == H.TEST[1], ] # one two-moment result per replication
    summ[[length(summ) + 1]] <- data.frame(design = run$d, estimator = run$est, H = m / 12, reps = run$reps,
      bound.share = mean(x$bound), tested = sum(ok),
      reject.5 = mean(x$p[ok] < 0.05), reject.10 = mean(x$p[ok] < 0.10),
      bias0 = mean(x$beta0[!x$bound]) - 1, bias1 = mean(x$beta1[!x$bound]) - 1, # true beta_MKT = 1 in every design
      mean.diff = mean(x$diff[!x$bound]), sd.diff = sd(x$diff[!x$bound]),
      se.ratio = mean(x$se.diff[ok]) / sd(x$diff[ok]), median.corr = median(x$corr[ok]),
      tm.reject.5 = mean(tm$tm.p.max < 0.05, na.rm = TRUE), tm.tested = sum(is.finite(tm$tm.p.max)),
      stringsAsFactors = FALSE)
  }
}
summ <- do.call(rbind, summ)
summ$git.commit <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1], error = function(e) NA)
write.csv(summ, "results/mc_step4/horizon_invariance_mc.csv", row.names = FALSE)
cat("wrote results/mc_step4/horizon_invariance_mc.csv\n")
print(format(summ[, setdiff(names(summ), "git.commit")], digits = 3), row.names = FALSE)

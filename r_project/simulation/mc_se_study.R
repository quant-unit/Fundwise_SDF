# Step 4 Monte Carlo: which standard error for U_{v,H}? (plan log entry
# "Step 4, first analysis": with q_v = n_v the sandwich SEs are too small).
# On the Buyout core replications (same draws, the stored estimates; no
# refit of the full-sample estimate), single-factor beta, horizons 0 and 15
# years, for U (q_v = n_v and q_v = 1) and L2 on vintage portfolios:
#   - Bartlett sandwich with S in {7 (= ceil(1.3 sqrt V), the estimator),
#     10, 14, V}, fixed-b critical value at b = S / V;
#   - independent vintages (no kernel), normal critical value;
#   - leave-one-vintage-out jackknife (refit by Brent on the feasible
#     interval of each subsample), critical values N(0,1) and t(V - 1).
# Reports the mean SE / Monte Carlo sd and the rejection rate of the true
# beta at 5%. Run from r_project/:
#   LANG=en_US.UTF-8 Rscript simulation/mc_se_study.R [replications] [workers]
source("simulation/mc_driver.R")
args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
workers <- if (length(args) >= 2) as.integer(args[2]) else 4
pub <- mc.public.factors(); D <- mc.load.designs(); spec <- D$core_BO
rows <- mc.read("core_BO")
rows <- rows[rows$model == "MKT" & rows$max.month %in% c(0, 180) & rows$estimator %in% c("U.n", "U.eq", "L2.vyp") & rows$rep <= R, ]
S.grid <- c(7, 10, 14, 25)

one.rep <- function(s) {
  dr <- mc.draw(spec$design, s, pub)
  fps <- list(fund = mc.as.fast(dr, pub, "fund"), vyp = mc.as.fast(dr, pub, "vyp"))
  d <- rows[rows$rep == s, ]
  do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    est <- MC.ESTIMATORS[[d$estimator[i]]]; fp <- fps[[est$level]]; m <- d$max.month[i]
    par <- c(MKT = d$MKT[i])
    out <- data.frame(rep = s, estimator = d$estimator[i], H = m / 12, beta = par[["MKT"]])
    for (S in S.grid) {
      MC.SETTINGS$bandwidth <<- S
      cv <- mc.sandwich(est, fp, par, m)
      out[[paste0("se.S", S)]] <- unname(cv$SE[1]); out[[paste0("b.S", S)]] <- cv$hac.b
    }
    MC.SETTINGS$bandwidth <<- "llsw"
    # independent vintages: S -> 0 (own-unit products only)
    units.cv <- {
      MC.SETTINGS$bandwidth <<- 0.5 # Bartlett weight 0 for distinct vintages, 1 for the same vintage
      x <- mc.sandwich(est, fp, par, m); MC.SETTINGS$bandwidth <<- "llsw"; x
    }
    out$se.indep <- unname(units.cv$SE[1])
    # leave-one-vintage-out jackknife (vintage groups; units of a vintage dropped together)
    vint <- sort(unique(fp$vintage))
    th <- vapply(vint, function(v) {
      sub <- fast.subset(fp, fp$vintage != v)
      con <- mc.domain(sub, "MKT"); br <- feasible.interval(con)
      f <- mc.criterion(est, sub, m)
      optimize(function(b) f(c(MKT = b)), br, tol = 1e-8)$minimum
    }, 0)
    V <- length(vint)
    out$se.jk <- sqrt((V - 1) / V * sum((th - mean(th))^2)); out$V <- V
    out
  }))
}
t0 <- Sys.time()
res <- do.call(rbind, parallel::mclapply(sort(unique(rows$rep)), function(s) tryCatch(one.rep(s), error = function(e) NULL), mc.cores = workers, mc.preschedule = FALSE))
cat("replications:", length(unique(res$rep)), "of", length(unique(rows$rep)), "|", round(as.numeric(difftime(Sys.time(), t0, units = "mins"))), "min\n")

rej <- function(t, crit) mean(abs(t) > crit, na.rm = TRUE)
summ <- do.call(rbind, lapply(split(res, list(res$estimator, res$H), drop = TRUE), function(d) {
  sdb <- sd(d$beta); t <- function(se) (d$beta - 1) / se
  x <- data.frame(estimator = d$estimator[1], H = d$H[1], reps = nrow(d), sd = sdb)
  for (S in S.grid) {
    b <- d[[paste0("b.S", S)]][1]
    x[[paste0("ratio.S", S)]] <- mean(d[[paste0("se.S", S)]]) / sdb
    x[[paste0("rej.S", S)]] <- rej(t(d[[paste0("se.S", S)]]), fixedb.crit.t(0.05, b))
  }
  x$ratio.indep <- mean(d$se.indep) / sdb; x$rej.indep <- rej(t(d$se.indep), qnorm(0.975))
  x$ratio.jk <- mean(d$se.jk) / sdb; x$rej.jk.norm <- rej(t(d$se.jk), qnorm(0.975)); x$rej.jk.t <- rej(t(d$se.jk), qt(0.975, d$V[1] - 1))
  x
}))
dir.create("simulation/data_out_mc/se_study", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(rows = res, summary = summ, commit = system("git rev-parse --short HEAD", intern = TRUE)), "simulation/data_out_mc/se_study/se_study.rds")
print(summ, digits = 3, row.names = FALSE)

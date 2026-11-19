# Step 4 Monte Carlo: extension rule (revision plan, Section I, implementation
# spec): "a margin is extended to 1,000 replications only if a statement in
# the paper rests on a difference smaller than two Monte Carlo SEs".
# Run from r_project/:  LANG=en_US.UTF-8 Rscript simulation/mc_extension.R
# Every statement drawn from the margins (250 replications) and the grid (500)
# in the plan log entry "Step 4, first analysis" is tested against its
# threshold: z = (estimate - threshold) / Monte Carlo SE, signed so that the
# statement holds for z > 0; the statement needs more replications if
# |z| < 2. Differences between estimators use the paired SE (same draws).
# (g) The persistent-shock statements of the additions (log entry "Step 4,
# additions: production and first look"): m4p_shared(2) against m4p_ind(2),
# the same draws with independent shocks, paired by replication.
# Writes results/mc_step4/extension_check.csv.
source("simulation/mc_driver.R")
D <- mc.load.designs()
margins <- grep("^m[0-9]+_", names(D), value = TRUE) # the first production's margins (not m4p)
grid <- grep("^g_", names(D), value = TRUE)
rows <- lapply(setNames(c(margins, grid), c(margins, grid)), mc.read)
out <- list()
add <- function(statement, design, H, estimate, se, threshold, direction) {
  z <- direction * (estimate - threshold) / se
  out[[length(out) + 1]] <<- data.frame(statement = statement, design = design, H = H, estimate = estimate, mcse = se,
                                         threshold = threshold, z = z, holds = z > 0, needs.more = abs(z) < 2)
}
sel <- function(r, e, m, h) r[r$estimator == e & r$model == m & r$max.month == h, ]
for (id in margins) {
  r <- rows[[id]]; th <- D[[id]]$design$theta0
  ests <- unique(r$estimator); hs <- sort(unique(r$max.month)); models <- unique(r$model)
  for (h in hs[hs <= 240]) {
    # (a) U's (alpha, beta) minimizer on a bound in under 5% of replications (0-20 years)
    for (e in intersect(c("U.eq", "U.n"), ests)) if ("Alpha" %in% models) {
      b <- sel(r, e, "Alpha", h)$at.bound; p <- mean(b)
      add("U (alpha, beta) on a bound < 5%", id, h / 12, p, sqrt(max(p, 0.01) * (1 - max(p, 0.01)) / length(b)), 0.05, -1)
    }
    # (b) U.n over-rejects the true beta (single factor): rejection > 5%
    if ("U.n" %in% ests && "MKT" %in% models) {
      d <- sel(r, "U.n", "MKT", h); t <- (d$MKT - th[["MKT"]]) / d$SE.MKT; ok <- is.finite(t) & is.finite(d$hac.b)
      p <- mean(fixedb.pvalue(t[ok]^2, d$hac.b[ok], 1) < 0.05)
      add("U (q_v = n_v) rejects the true beta in more than 5%", id, h / 12, p, sqrt(0.05 * 0.95 / sum(ok)), 0.05, +1)
    }
  }
  # (c) at H = 0 the L2 bias exceeds the U bias (single factor; paired), except M6
  if (all(c("L2.vyp", "U.eq") %in% ests) && "MKT" %in% models) {
    a <- sel(r, "L2.vyp", "MKT", 0); b <- sel(r, "U.eq", "MKT", 0); m <- merge(a[, c("rep", "MKT")], b[, c("rep", "MKT")], by = "rep")
    d <- abs(m$MKT.x - th[["MKT"]]) - abs(m$MKT.y - th[["MKT"]])
    add("H = 0: |L2 error| > |U error| on average (paired)", id, 0, mean(d), sd(d) / sqrt(nrow(m)), 0, +1)
  }
  # (d) M4: the U target moves away from theta_0 (bias of U.eq at H = 0 and 30 years)
  if (grepl("^m4", id)) for (h in c(0, 360)) {
    x <- sel(r, "U.eq", "MKT", h)$MKT
    add("M4: U (q_v = 1) biased at H = 0 and 30 years", id, h / 12, abs(mean(x) - th[["MKT"]]), sd(x) / sqrt(length(x)), 0, +1)
  }
  # (e) M6: L2 overestimates the second loading at H = 0; U (q_v = 1) within 0.05 of it at 15 years
  if (grepl("^m6", id)) {
    f <- setdiff(names(th), "MKT"); x <- sel(r, "L2.vyp", f, 0)$Coef
    add("M6: L2 overestimates the second loading at H = 0", id, 0, mean(x) - th[[f]], sd(x) / sqrt(length(x)), 0, +1)
    x <- sel(r, "U.eq", f, 180)$Coef
    add("M6: U (q_v = 1) second loading within 0.05 at 15 years", id, 15, abs(mean(x) - th[[f]]), sd(x) / sqrt(length(x)), 0.05, -1)
  }
}
# (f) grid: the bias of beta at H = 0 falls from 5 to 200 funds per vintage, for each V
for (V in c(10, 20, 39)) {
  a <- rows[[sprintf("g_n5_V%d", V)]]; b <- rows[[sprintf("g_n200_V%d", V)]]
  x <- sel(a, "L2.vyp", "Alpha", 0)$MKT - 1; y <- sel(b, "L2.vyp", "Alpha", 0)$MKT - 1
  add("G: |bias| at H = 0 smaller with 200 than with 5 funds per vintage", sprintf("V = %d", V), 0, abs(mean(x)) - abs(mean(y)),
      sqrt(var(x) / length(x) + var(y) / length(y)), 0, +1)
}
# (g) persistent shared shocks (rho_W = 0.6) against independent ones, paired:
# every estimator rejects the true beta more often at H = 0 and 15 years; the
# sd of the estimate is larger (SE of the sd difference from a paired bootstrap)
rej <- function(d) { t <- (d$MKT - 1) / d$SE.MKT; p <- rep(NA, nrow(d)); ok <- is.finite(t) & is.finite(d$hac.b)
  p[ok] <- fixedb.pvalue(t[ok]^2, d$hac.b[ok], 1) < 0.05; p }
for (k in c("", "2")) {
  ind <- mc.read(paste0("m4p_ind", k)); per <- mc.read(paste0("m4p_shared", k)); id <- paste0("m4p_shared", k)
  for (e in c("L2.vyp", "U.eq", "U.n", "IV")) for (h in c(0, 180)) {
    f <- function(x) { d <- sel(x, e, "MKT", h); data.frame(rep = d$rep, rj = rej(d), b = d$MKT) }
    m <- merge(f(ind), f(per), by = "rep"); mt <- m[!is.na(m$rj.x) & !is.na(m$rj.y), ]
    d <- mt$rj.y - mt$rj.x
    add(paste0("persistence raises rejection of the true beta (", e, ", paired)"), id, h / 12, mean(d), sd(d) / sqrt(nrow(mt)), 0, +1)
    set.seed(1); bs <- replicate(500, { i <- sample(nrow(m), replace = TRUE); sd(m$b.y[i]) - sd(m$b.x[i]) })
    add(paste0("persistence widens the sd of the estimate (", e, ", paired)"), id, h / 12, sd(m$b.y) - sd(m$b.x), sd(bs), 0, +1)
  }
}
res <- do.call(rbind, out)
dir.create("results/mc_step4", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "results/mc_step4/extension_check.csv", row.names = FALSE)
cat(nrow(res), "checks;", sum(!res$holds), "statements do not hold as written;", sum(res$needs.more), "within 2 Monte Carlo SEs\n")
print(res[res$needs.more | !res$holds, ], digits = 3, row.names = FALSE)

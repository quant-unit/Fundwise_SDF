# Step 4 SE diagnostic, analysis of diag_se.rds (output of diag_se.R). Run from r_project/:
#   LANG=en_US.UTF-8 Rscript <this file> <diag_se.rds>
source("helper/fixedb.R"); source("helper/inference.R")
x <- readRDS(commandArgs(TRUE)[1]); sc <- x$scores; rows <- x$rows
options(width = 200)
bart.dm <- function(s, v, S) { d <- s - mean(s); W <- outer(v, v, bartlett.weight, bw = S); as.numeric(t(d) %*% W %*% d) / length(s) }
gauss.bench <- function(v, S, sd.v, reps = 20000, seed = 1) {
  set.seed(seed); b <- S / length(v)
  r <- t(replicate(reps, { s <- rnorm(length(v), 0, sd.v); n <- length(s)
    c(t = sqrt(n) * mean(s) / sqrt(bart.dm(s, v, S)), psi = bart.dm(s, v, S), g = sqrt(n) * mean(s)) }))
  c(size = mean(fixedb.q1.pvalue(r[, "t"]^2, b) < 0.05), psi.ratio = mean(r[, "psi"]) / var(r[, "g"]),
    psi.cv = sd(r[, "psi"]) / mean(r[, "psi"]))
}
out <- list(); lagc <- list()
for (id in unique(sc$design)) for (m in c(0, 180)) {
  a <- sc[sc$design == id & sc$max.month == m, ]
  r <- rows[rows$design == id & rows$max.month == m, ]; r <- r[order(r$rep), ]
  S <- unique(r$hac.S); b <- unique(r$hac.b); stopifnot(length(S) == 1, length(b) == 1)
  per <- split(a, a$rep); per <- per[as.character(r$rep)]
  n <- vapply(per, nrow, 0)
  g <- vapply(per, function(p) sqrt(nrow(p)) * mean(p$s), 0)
  hbar <- vapply(per, function(p) mean(p$h), 0)
  psi <- vapply(per, function(p) bart.dm(p$s, p$vint, S), 0)
  tLM <- g / sqrt(psi)
  t.act <- (r$MKT - 1) / r$SE.MKT
  lin <- -vapply(per, function(p) mean(p$s), 0) / hbar # one Newton step from beta = 1
  # per-vintage score variance across replications (heteroskedasticity) and Gaussian benchmarks
  vs <- tapply(a$s, a$vint, var); vv <- as.numeric(names(vs))
  gb.iid <- gauss.bench(vv, S, rep(1, length(vv)))
  gb.het <- gauss.bench(vv, S, sqrt(as.numeric(vs)))
  z <- (a$s - ave(a$s, a$vint)) / sqrt(ave(a$s, a$vint, FUN = var))
  # cross-vintage correlation of scores across replications, by year gap
  wide <- reshape(a[, c("rep", "vint", "s")], idvar = "rep", timevar = "vint", direction = "wide")
  C <- cor(wide[, -1], use = "pairwise.complete.obs"); yrs <- as.numeric(sub("s.", "", colnames(C)))
  gap <- abs(outer(yrs, yrs, "-"))
  lagc[[paste(id, m)]] <- sapply(1:5, function(k) mean(C[gap == k]))
  out[[paste(id, m)]] <- data.frame(design = id, H = m / 12, reps = nrow(r), V = median(n), S = S, b = round(b, 3),
    size.actual = mean(fixedb.q1.pvalue(t.act^2, b) < 0.05, na.rm = TRUE),
    size.score = mean(fixedb.q1.pvalue(tLM^2, b) < 0.05),
    gauss.iid.size = gb.iid[["size"]], gauss.het.size = gb.het[["size"]],
    psi.ratio = mean(psi) / var(g), gauss.iid.psi.ratio = gb.iid[["psi.ratio"]], gauss.het.psi.ratio = gb.het[["psi.ratio"]],
    psi.cv = sd(psi) / mean(psi), gauss.iid.psi.cv = gb.iid[["psi.cv"]], gauss.het.psi.cv = gb.het[["psi.cv"]],
    neff.var = sum(vs)^2 / sum(vs^2), top.share = max(vs) / sum(vs),
    kurt.score = mean(z^4), cor.t = cor(t.act, -sign(hbar) * tLM),
    sd.ratio.lin = sd(lin) / sd(r$MKT - 1), h.cv = sd(hbar) / abs(mean(hbar)),
    se.cv = sd(r$SE.MKT) / mean(r$SE.MKT), cor.se.err = cor(r$SE.MKT, abs(r$MKT - 1)),
    t.q975 = unname(quantile(abs(t.act), 0.95)), tLM.q975 = unname(quantile(abs(tLM), 0.95)), crit = fixedb.crit.t(0.05, b))
}
res <- do.call(rbind, out); rownames(res) <- NULL
num <- vapply(res, is.numeric, TRUE); res[num] <- lapply(res[num], round, 3)
print(t(res))
cat("\ncross-vintage correlation of scores, year gaps 1-5:\n"); print(round(do.call(rbind, lagc), 3))
cat("\nper-vintage score variance shares, core_BO H=0 (vintage: share, n funds):\n")
a <- sc[sc$design == "core_BO" & sc$max.month == 0, ]; vs <- tapply(a$s, a$vint, var); print(round(sort(vs / sum(vs), decreasing = TRUE)[1:6], 3))
a <- sc[sc$design == "w_V25_n15" & sc$max.month == 0, ]; vs <- tapply(a$s, a$vint, var); cat("w_V25_n15 H=0 top shares:\n"); print(round(sort(vs / sum(vs), decreasing = TRUE)[1:6], 3))

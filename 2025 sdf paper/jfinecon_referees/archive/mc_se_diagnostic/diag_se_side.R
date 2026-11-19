# Step 4 SE diagnostic: rejections by side; SE at the estimate against at the truth. Run from r_project/:
#   LANG=en_US.UTF-8 Rscript <this file> <diag_se.rds>
source("helper/fixedb.R"); source("helper/inference.R")
x <- readRDS(commandArgs(TRUE)[1]); sc <- x$scores; rows <- x$rows; options(width = 200)
bart.dm <- function(s, v, S) { d <- s - mean(s); W <- outer(v, v, bartlett.weight, bw = S); as.numeric(t(d) %*% W %*% d) / length(s) }
out <- list()
for (id in unique(sc$design)) for (m in c(0, 180)) {
  a <- sc[sc$design == id & sc$max.month == m, ]; r <- rows[rows$design == id & rows$max.month == m, ]; r <- r[order(r$rep), ]
  S <- unique(r$hac.S); b <- unique(r$hac.b); per <- split(a, a$rep)[as.character(r$rep)]
  n <- vapply(per, nrow, 0); gb <- vapply(per, function(p) mean(p$s), 0); hb <- vapply(per, function(p) mean(p$h), 0)
  psi <- vapply(per, function(p) bart.dm(p$s, p$vint, S), 0)
  se.lm <- sqrt(psi / n) / abs(hb)            # SE of the linearized estimator, all at beta = 1
  tlm <- -gb / hb / se.lm                      # signed like beta-hat - 1
  ta <- (r$MKT - 1) / r$SE.MKT
  rej.a <- fixedb.q1.pvalue(ta^2, b) < 0.05; rej.l <- fixedb.q1.pvalue(tlm^2, b) < 0.05
  d <- r$MKT - 1; lr <- log(r$SE.MKT / se.lm)
  # Wald with the SE at beta = 1 but the actual estimate, and with the actual SE but the linearized estimate
  t.mix1 <- d / se.lm; t.mix2 <- (-gb / hb) / r$SE.MKT
  out[[paste(id, m)]] <- data.frame(design = id, H = m / 12,
    rej.actual.up = mean(rej.a & d > 0), rej.actual.down = mean(rej.a & d < 0),
    rej.lm.up = mean(rej.l & tlm > 0), rej.lm.down = mean(rej.l & tlm < 0),
    rej.est.with.se.at.truth = mean(fixedb.q1.pvalue(t.mix1^2, b) < 0.05),
    rej.lin.est.with.actual.se = mean(fixedb.q1.pvalue(t.mix2^2, b) < 0.05),
    med.se.ratio = median(exp(lr)), slope.logSEratio.on.err = unname(coef(lm(lr ~ d))[2]), sd.err = sd(d),
    se.ratio.when.rejected = median(exp(lr[rej.a])), se.ratio.not.rejected = median(exp(lr[!rej.a])))
}
res <- do.call(rbind, out); rownames(res) <- NULL; num <- vapply(res, is.numeric, TRUE); res[num] <- lapply(res[num], round, 3); print(t(res))

# Check of the external reviews (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/domains_eg_h0.R" <stall.rds from stall_check.R>
# Second review's D2 proposal: floor 0.25 on the monthly gross benchmark return plus box [-4,4].
# Constrained multi-start for the three MKT+EG H=0 cells with a degenerate lower minimum,
# and the VC FW refit without vintage 2003 (the contaminated jackknife entry).
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R"))))); source(here("helper/objective.R")); source(here("helper/parallel_workers.R"))
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
res <- readRDS(commandArgs(TRUE)[1])
cells <- list(c("FW", "BO"), c("EW", "VC"), c("FW", "VC"))
dat <- lapply(cells, function(cl) {
  a <- read.csv(here(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_VYP_North America_NAV.csv", cl[1]))); a$Date <- as.Date(a$Date)
  v <- a[a$type == cl[2] & a$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
  U <- prepare.units(merge(pad.all.units(v, "first.cf", max(a$Date)), q, by = "Date", all.x = TRUE), c("RF", "MKT", "EG"))$units
  s <- q[q$Date >= min(v$Date) & q$Date <= max(a$Date), ]
  x <- res[res$w == cl[1] & res$Type == cl[2] & res$Factor == "EG" & res$H == 0, ]
  list(U = U, s = s, rep = c(x$MKT.row, x$Coef.row), edge = c(x$MKT.best, x$Coef.best), vint = sapply(U, `[[`, "vintage"))
})
domains <- list(none = c(g = -Inf, box = Inf), box4 = c(-Inf, 4), floor.25 = c(0.25, Inf), floor.25_box4 = c(0.25, 4), floor.40 = c(0.40, Inf), floor.50 = c(0.50, Inf))
fit <- function(d, dom, U) {
  minG <- function(p) min(1 + d$s$RF + p[1] * d$s$MKT + p[2] * d$s$EG)
  crit <- function(p) { if (any(abs(p) > dom[2]) || minG(p) < dom[1]) return(1e10)
    e <- vapply(U, function(u) unit.horizon.error(u, 0, c(RF = 1, MKT = p[1], EG = p[2])), 0); if (any(!is.finite(e))) 1e10 else sqrt(sum(e^2)) / length(e) }
  set.seed(7); st <- rbind(d$rep, c(1, 0), cbind(runif(40, -2, 6), runif(40, -10, 15)))
  st <- st[apply(st, 1, crit) < 1e9, , drop = FALSE]
  f <- t(apply(st, 1, function(s0) { o <- optim(s0, crit, method = "Nelder-Mead", control = list(reltol = 1e-12, maxit = 5000)); c(o$par, o$value) }))
  b <- f[which.min(f[, 3]), ]; c(MKT = b[1], EG = b[2], Q = b[3], minG = minG(b[1:2]), starts = nrow(st))
}
jobs <- expand.grid(cell = 1:3, dom = names(domains), loo = c(FALSE, TRUE), stringsAsFactors = FALSE); jobs <- jobs[!jobs$loo | jobs$cell == 3, ]
out <- parallel::mclapply(seq_len(nrow(jobs)), function(k) { j <- jobs[k, ]; d <- dat[[j$cell]]; U <- if (j$loo) d$U[d$vint != 2003] else d$U
  r <- fit(d, domains[[j$dom]], U); Qrep <- if (j$loo) NA else { e <- vapply(U, function(u) unit.horizon.error(u, 0, c(RF = 1, MKT = d$rep[1], EG = d$rep[2])), 0); sqrt(sum(e^2)) / length(e) }
  data.frame(cell = paste(cells[[j$cell]], collapse = " "), sample = if (j$loo) "without 2003" else "full", domain = j$dom, MKT = r[1], EG = r[2], Q = r[3], Q.reported = Qrep, gain.vs.reported = if (j$loo) NA else (Qrep - r[3]) / Qrep, minG = r[4]) },
  mc.cores = n.parallel.workers(nrow(jobs)))
o <- do.call(rbind, out); num <- sapply(o, is.numeric); o[num] <- lapply(o[num], signif, 4); print(o, row.names = FALSE)
cat("reported points: "); for (d in dat) cat(sprintf("(%.3f, %.3f) ", d$rep[1], d$rep[2])); cat("\n")

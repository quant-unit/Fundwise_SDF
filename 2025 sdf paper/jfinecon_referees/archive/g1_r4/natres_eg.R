# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/natres_eg.R" 
# Is the reported NATRES MKT+EG (FW, H0) estimate a minimum? Compare objective at
# the reported point with a Nelder-Mead restart and a multi-start search.
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R")))))
suppressMessages(source(here("helper/result_cache.R"))); source(here("helper/objective.R"))
rr <- read.result.cache("results/data_out_2026-emp-c/cache_q_factors_preqin_FW_VYP_North America", verbose = FALSE)
row <- function(f, H) rr[rr$Type == "NATRES" & rr$CV.key == "ALL" & rr$max.month == H & rr$lambda == 0 & rr$Factor == f, ]
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
all <- read.csv(here("empirical/data_prepared_2026_c/preqin_cashflows_FW_VYP_North America_NAV.csv")); all$Date <- as.Date(all$Date)
T.max <- max(all$Date)
v <- all[all$type == "NATRES" & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
U <- prepare.units(merge(pad.all.units(v, "first.cf", T.max), q, by = "Date", all.x = TRUE), c("RF", "MKT", "EG"))$units
cat("units:", length(U), "\n")
Q <- function(b, m = 0) { e <- vapply(U, function(u) horizon.error(u$cf, exp(cumsum(log(1 + u$X %*% c(1, b)))), m), 0); if (any(!is.finite(e))) 1e10 else sum(e^2) }
r1 <- row("MKT", 0); b1 <- optimize(function(b) Q(c(b, 0)), c(-1, 4), tol = 1e-10)$minimum
cat(sprintf("reproduce 1-factor: beta %.6f (row %.6f), Q %.7f (row %.7f)\n", b1, r1$MKT, Q(c(b1, 0)), r1$value))
for (H in c(0, 1)) {
  r2 <- row("EG", H); p <- c(r2$MKT, r2$EG)
  cat(sprintf("\nH=%d reported: MKT %.4f EG %.4f, Q(row value) %.7f, Q recomputed %.7f\n", H, p[1], p[2], r2$value, Q(p, H)))
  nm <- optim(p, Q, m = H, method = "Nelder-Mead")
  cat(sprintf("NM restart:  MKT %.4f EG %.4f, Q %.7f, move %.3f\n", nm$par[1], nm$par[2], nm$value, max(abs(nm$par - p))))
  set.seed(1); st <- cbind(runif(200, -2, 6), runif(200, -6, 6))
  ms <- t(apply(st, 1, function(s) { o <- optim(s, Q, m = H, method = "Nelder-Mead", control = list(reltol = 1e-12, maxit = 5000)); c(o$par, o$value) }))
  best <- ms[which.min(ms[, 3]), ]
  cat(sprintf("multi-start best of 200: MKT %.4f EG %.4f, Q %.7f\n", best[1], best[2], best[3]))
  near <- ms[ms[, 3] < best[3] + 1e-4, , drop = FALSE]
  cat(sprintf("starts within 1e-4 of best: %d; their MKT range %.3f to %.3f, EG range %.3f to %.3f\n", nrow(near), min(near[,1]), max(near[,1]), min(near[,2]), max(near[,2])))
}

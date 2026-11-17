# Run from r_project/:  Rscript "../2025 sdf paper/jfinecon_referees/archive/check_horizon_reweighting.R" VC 360
# Horizon error = PV x A (exact), and what drives the 30-year VC estimate.
# VC, EW, North America, repaired pools, vintages <= 2010, single-factor MKT.
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R")))))
source(here("helper/objective.R"))
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv"))
q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
all <- read.csv(here("empirical/data_prepared_2026_c/preqin_cashflows_EW_VYP_North America_NAV.csv")); all$Date <- as.Date(all$Date)
T.max <- max(all$Date)
ty <- commandArgs(TRUE)[1]; H <- as.integer(commandArgs(TRUE)[2])
vyp <- all[all$type == ty & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]
vyp$Fund.ID <- paste(vyp$Fund.ID, vyp$type, sep = "_")
units.for <- function(end) {
  p <- if (end == "data") pad.all.units(vyp, "first.cf", T.max)
       else do.call(rbind, lapply(split(vyp, vyp$Fund.ID), function(d) pad.all.units(d, "first.cf", max(d$Date))))
  prepare.units(merge(p, q, by = "Date", all.x = TRUE), c("RF", "MKT"))$units
}
U <- list(a = units.for("data"), b = units.for("own"))
pv.A <- function(u, beta, m) {
  ret <- exp(cumsum(log(1 + u$X %*% c(1, beta))))
  c(PV = sum(u$cf / ret), A = mean(ret[seq_len(min(max(m, 1), length(ret)))]),
    e = horizon.error(u$cf, ret, m), n = length(ret))
}
tab <- function(units, beta, m) do.call(rbind, lapply(units, pv.A, beta = beta, m = m))
Q <- function(units, beta, m) sum(tab(units, beta, m)[, "e"]^2)

# 1) identity e = PV x A, every vintage, both endpoints, several beta and H
gap <- 0
for (k in names(U)) for (b in c(-0.5, 0.5, 1, 1.5, 2.5)) for (m in c(0, 60, 180, 360)) {
  x <- tab(U[[k]], b, m); gap <- max(gap, abs(x[, "e"] - x[, "PV"] * x[, "A"]) / pmax(1e-12, abs(x[, "e"])))
}
cat(sprintf("%s: identity e = PV*A, max relative gap %.1e (2 endpoints x 5 beta x 4 horizons x %d vintages)\n", ty, gap, length(U$a)))
cat(sprintf("PV identical under both endpoints: max |diff| %.1e\n",
  max(abs(tab(U$a, 1.2, 0)[, "PV"] - tab(U$b, 1.2, 0)[, "PV"]))))

# 2) estimates: H = 0 (A = 1: pure pricing), and H under each endpoint
est <- function(units, m) optimize(function(b) Q(units, b, m), c(-1, 4), tol = 1e-10)$minimum
b0 <- est(U$a, 0); ba <- est(U$a, H); bb <- est(U$b, H)
cat(sprintf("beta: H=0 %.3f | H=%d m: data end %.3f, own end %.3f\n", b0, H, ba, bb))

# 3) slope of Q under the data end at the own-end estimate, split into
#    pricing term  sum 2 A^2 PV dPV  and scale term  sum 2 PV^2 A dA
d <- 1e-6; x <- tab(U$a, bb, H); xp <- tab(U$a, bb + d, H); xm <- tab(U$a, bb - d, H)
dPV <- (xp[, "PV"] - xm[, "PV"]) / (2 * d); dA <- (xp[, "A"] - xm[, "A"]) / (2 * d)
t1 <- sum(2 * x[, "A"]^2 * x[, "PV"] * dPV); t2 <- sum(2 * x[, "PV"]^2 * x[, "A"] * dA)
cat(sprintf("data end, slope at beta = %.3f: pricing term %.4g, scale term %.4g (positive slope: Q falls as beta falls)\n", bb, t1, t2))

# 4) vintages: window months, A at beta = 1 under both endpoints, share of Q
xa <- tab(U$a, 1, H); xb <- tab(U$b, 1, H)
v <- data.frame(vintage = sapply(U$a, `[[`, "vintage"), months.own = xb[, "n"], months.data = xa[, "n"],
  PV = round(xa[, "PV"], 2), A.own = round(xb[, "A"], 2), A.data = round(xa[, "A"], 2),
  share.Q.own = round(100 * xb[, "e"]^2 / sum(xb[, "e"]^2), 1), share.Q.data = round(100 * xa[, "e"]^2 / sum(xa[, "e"]^2), 1))
print(v[order(v$vintage), ], row.names = FALSE)

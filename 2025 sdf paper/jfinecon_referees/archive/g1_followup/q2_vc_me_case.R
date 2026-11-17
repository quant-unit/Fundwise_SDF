# G1 follow-up (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_followup/q2_vc_me_case.R"
# The case cited in the G1 open item: VC EW, MKT + ME, H = 20 years.
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R")))))
source(here("helper/objective.R"))
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
all <- read.csv(here("empirical/data_prepared_2026_c/preqin_cashflows_EW_VYP_North America_NAV.csv")); all$Date <- as.Date(all$Date)
v <- all[all$type == "VC" & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
U <- prepare.units(merge(pad.all.units(v, "first.cf", max(all$Date)), q, by = "Date", all.x = TRUE), c("RF", "MKT", "ME"))$units
m <- 240
Q <- function(p, units) { e <- vapply(units, function(u) horizon.error(u$cf, exp(cumsum(log(1 + u$X %*% c(1, p)))), m), 0); if (any(!is.finite(e))) 1e10 else sum(e^2) }
fit <- function(units, start) { o <- optim(start, Q, units = units, control = list(reltol = 1e-14, maxit = 5000)); optim(o$par, Q, units = units, method = "BFGS", control = list(reltol = 1e-14))$par }
g1 <- read.csv(here("results/figures_bovc_us_g1/empirical_VC_US_Multi_data.csv")); g1 <- g1[g1$method == "Asymptotic" & g1$Factor == "ME" & g1$weighting == "EW" & g1$max.month == m, ]
p <- fit(U, c(1, 0)); cat(sprintf("full sample: MKT %.4f ME %.4f | G1 table: MKT %.4f ME %.4f (SE.MKT %.3f, SE.ME %.3f, indep %.3f)\n", p[1], p[2], g1$MKT, g1$Coef, g1$SE.MKT, g1$SE.Coef, g1$SE.Coef.indep))
pj <- t(sapply(seq_along(U), function(i) fit(U[-i], p))); V <- length(U)
se <- sqrt((V - 1) / V * colSums(sweep(pj, 2, colMeans(pj))^2))
cat(sprintf("jackknife SE: MKT %.3f, ME %.3f | ME leave-one-out range %.3f to %.3f (vintage removed at max: %d)\n", se[1], se[2], min(pj[, 2]), max(pj[, 2]), sapply(U, `[[`, "vintage")[which.max(pj[, 2])]))

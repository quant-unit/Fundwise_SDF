# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/q3_corr_llsw.R" 
# Q3 check: correlation of (beta, alpha) from the sandwich used for the G1 SEs
# (helper/inference.R, unit loss e^2), now with the bandwidth rule S = ceil(1.3 sqrt V) = 7; SEs checked against the r4 G1 table.
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R")))))
source(here("helper/objective.R")); source(here("helper/inference.R"))
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1; q$Alpha <- 1
for (w in c("EW", "FW")) {
  all <- read.csv(here(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_VYP_North America_NAV.csv", w))); all$Date <- as.Date(all$Date)
  v <- all[all$type == "BO" & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
  U <- prepare.units(merge(pad.all.units(v, "first.cf", max(all$Date)), q, by = "Date", all.x = TRUE), c("RF", "MKT", "Alpha"))$units
  g1 <- read.csv(here("results/figures_bovc_us_g1/empirical_BO_US_Alpha_data.csv")); g1 <- g1[g1$method == "Asymptotic" & g1$Factor == "Alpha" & g1$weighting == w, ]
  for (m in c(0, 60, 120, 180, 240)) {
    gg <- g1[g1$max.month == m, ]; par <- c(MKT = gg$MKT, Alpha = gg$Coef)
    loss <- function(p, u) horizon.error(u$cf, exp(cumsum(log(1 + u$X %*% c(1, p)))), m)^2
    s <- sandwich.covariance(U, par, loss, sapply(U, `[[`, "vintage"), bw = "llsw")
    cat(sprintf("BO %s H=%3d: sandwich corr(beta, alpha) %.3f | SE beta %.4f (G1 %.4f), SE alpha %.5f (G1 %.5f)\n", w, m,
      s$CoVM[1, 2] / sqrt(s$CoVM[1, 1] * s$CoVM[2, 2]), s$SE[1], gg$SE.MKT, s$SE[2], gg$SE.Coef))
  }
}

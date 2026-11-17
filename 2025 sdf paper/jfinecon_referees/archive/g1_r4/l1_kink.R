# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/l1_kink.R" 
# L1 single-factor: is the estimate at a kink (one vintage with error ~ 0), the same vintage across horizons?
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R"))))); source(here("helper/objective.R"))
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
for (ty in c("BO", "VC")) {
  d <- read.csv(here(sprintf("results/figures_bovc_us_g1/empirical_%s_US_MKT_L1_data.csv", ty))); d <- d[d$method == "Asymptotic" & d$Factor == "MKT" & d$max.month != 1, ]
  for (w in c("EW", "FW")) {
    all <- read.csv(here(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_VYP_North America_NAV.csv", w))); all$Date <- as.Date(all$Date)
    v <- all[all$type == ty & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
    U <- prepare.units(merge(pad.all.units(v, "first.cf", max(all$Date)), q, by = "Date", all.x = TRUE), c("RF", "MKT"))$units
    vint <- sapply(U, `[[`, "vintage")
    for (k in which(d$weighting == w)) {
      b <- d$MKT[k]; H <- d$max.month[k]
      e <- vapply(U, function(u) unit.horizon.error(u, H, c(RF = 1, MKT = b)), 0)
      L <- function(bb) mean(abs(vapply(U, function(u) unit.horizon.error(u, H, c(RF = 1, MKT = bb)), 0)))
      i <- order(abs(e))
      cat(sprintf("%s %s H=%3d beta %.5f SE %.1e | smallest |e| %.1e (vintage %d), next %.1e | median |e| %.3f | L(b)-L(b+-0.01) %.1e %.1e\n",
        ty, w, H, b, d$SE.MKT[k], abs(e[i[1]]), vint[i[1]], abs(e[i[2]]), median(abs(e)), L(b) - L(b - 0.01), L(b) - L(b + 0.01)))
    }
  }
}

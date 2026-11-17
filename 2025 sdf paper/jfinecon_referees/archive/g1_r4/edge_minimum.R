# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/edge_minimum.R" <stall.rds from stall_check.R>
# The lower MKT+EG minimum at H = 0: does it price funds, or discount cash flows away?
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R"))))); source(here("helper/objective.R"))
res <- readRDS(commandArgs(TRUE)[1])
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
for (cell in list(c("FW", "BO"), c("EW", "VC"), c("FW", "VC"))) {
  x <- res[res$w == cell[1] & res$Type == cell[2] & res$Factor == "EG" & res$H == 0, ]
  all <- read.csv(here(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_VYP_North America_NAV.csv", cell[1]))); all$Date <- as.Date(all$Date)
  v <- all[all$type == cell[2] & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
  U <- prepare.units(merge(pad.all.units(v, "first.cf", max(all$Date)), q, by = "Date", all.x = TRUE), c("RF", "MKT", "EG"))$units
  pts <- list(reported = c(x$MKT.row, x$Coef.row), lowest = c(x$MKT.best, x$Coef.best))
  for (lab in names(pts)) { pt <- pts[[lab]]
    p0 <- c(RF = 1, MKT = pt[1], EG = pt[2])
    g <- unlist(lapply(U, function(u) as.numeric(u$X[, names(p0)] %*% p0)))   # monthly SDF return
    e <- vapply(U, function(u) unit.horizon.error(u, 0, p0), 0)
    first <- vapply(U, function(u) u$cf[which(u$cf != 0)[1]], 0)                 # first nonzero cash flow
    yr1 <- vapply(U, function(u) { ret <- exp(cumsum(log(1 + u$X[, names(p0)] %*% p0))); r <- ret / ret[1]; sum(u$cf[seq_len(min(12, length(u$cf)))] / r[seq_len(min(12, length(u$cf)))]) }, 0)
    cat(sprintf("%s %s %-8s MKT %.3f EG %.3f | SDF return: mean %.3f/month, annual growth x%.1f | RMS e %.3f | RMS first CF %.3f | corr(e, first-year discounted CF) %.3f\n",
      cell[1], cell[2], lab, pt[1], pt[2], mean(g), exp(12 * mean(log(1 + g))), sqrt(mean(e^2)), sqrt(mean(first^2)), cor(e, yr1)))
  }
}

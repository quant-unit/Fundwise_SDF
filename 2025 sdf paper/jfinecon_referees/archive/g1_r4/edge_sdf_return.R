# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/edge_sdf_return.R" <stall.rds from stall_check.R>
suppressPackageStartupMessages(library(here))
res <- readRDS(commandArgs(TRUE)[1])
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(q$Date)
for (cell in list(c("FW", "BO"), c("EW", "VC"), c("FW", "VC"))) {
  x <- res[res$w == cell[1] & res$Type == cell[2] & res$Factor == "EG" & res$H == 0, ]
  for (lab in c("reported", "lowest")) {
    pt <- if (lab == "reported") c(x$MKT.row, x$Coef.row) else c(x$MKT.best, x$Coef.best)
    s <- q[q$Date >= as.Date("1969-01-01") & q$Date <= as.Date("2025-12-31"), ]   # wide window; the binding month is printed
    g <- 1 + s$RF + pt[1] * s$MKT + pt[2] * s$EG
    i <- which.min(g)
    cat(sprintf("%s %s %-8s MKT %6.3f EG %6.3f | min gross SDF return %.4f in %s | SD monthly %.3f\n", cell[1], cell[2], lab, pt[1], pt[2], g[i], format(s$Date[i], "%Y-%m"), sd(g)))
  }
}

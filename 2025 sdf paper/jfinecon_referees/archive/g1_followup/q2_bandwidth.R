# G1 follow-up (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_followup/q2_bandwidth.R" <out.rds of q2_jackknife.R>
# Does the Bartlett sandwich shrink because the scores sum to zero at the
# estimate? Psi(bw) = g' W(bw) g / n with g_v = d e_v^2 / d beta at beta-hat.
# Check 1: sqrt(Psi(12) / Psi(indep)) must equal the G1 ratio SE.MKT / SE.MKT.indep.
# Check 2: the ratio as a function of the bandwidth.
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R")))))
source(here("helper/objective.R")); source(here("helper/inference.R"))
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
r <- readRDS(commandArgs(TRUE)[1]); out <- list()
for (w in c("EW", "FW")) {
  all <- read.csv(here(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_VYP_North America_NAV.csv", w))); all$Date <- as.Date(all$Date)
  for (ty in c("BO", "VC")) {
    v <- all[all$type == ty & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
    U <- prepare.units(merge(pad.all.units(v, "first.cf", max(all$Date)), q, by = "Date", all.x = TRUE), c("RF", "MKT"))$units
    vin <- sapply(U, `[[`, "vintage")
    for (m in c(0, 120, 240, 360)) {
      x <- r[r$Type == ty & r$w == w & r$H == m, ]; b <- x$beta
      e2 <- function(bb) vapply(U, function(u) horizon.error(u$cf, exp(cumsum(log(1 + u$X %*% c(1, bb)))), m)^2, 0)
      g <- (e2(b + 1e-6) - e2(b - 1e-6)) / 2e-6
      psi <- function(bw) { W <- if (is.infinite(bw)) matrix(1, length(g), length(g)) else outer(vin, vin, bartlett.weight, bw = bw); as.numeric(t(g) %*% W %*% g) }
      p0 <- sum(g^2)
      out[[length(out) + 1]] <- data.frame(Type = ty, w = w, H = m, sum.g.rel = abs(sum(g)) / sqrt(p0), g1.ratio = x$se.bartlett / x$se.indep,
        bw3 = sqrt(psi(3) / p0), bw6 = sqrt(psi(6) / p0), bw12 = sqrt(psi(12) / p0), bw24 = sqrt(psi(24) / p0), bwInf = sqrt(max(psi(Inf), 0) / p0))
    }
  }
}
o <- do.call(rbind, out); cat(sprintf("check 1: max |sqrt(Psi12/Psi0) - G1 ratio| = %.4f\n", max(abs(o$bw12 - o$g1.ratio))))
o[, 4:10] <- round(o[, 4:10], 3); print(o, row.names = FALSE)

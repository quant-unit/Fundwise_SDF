# G1 follow-up (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_followup/q2_jackknife.R" <out.rds>   (then q2_bandwidth.R <out.rds>)
# Q2: effective number of vintages and leave-one-vintage-out jackknife SE vs
# the G1 sandwich SEs. Single-factor MKT, VYP, North America, vintages <= 2010,
# endpoint = common data end (as in G1), L2 criterion (argmin of sum e^2).
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R")))))
source(here("helper/objective.R"))
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
res <- list()
for (w in c("EW", "FW")) {
  all <- read.csv(here(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_VYP_North America_NAV.csv", w))); all$Date <- as.Date(all$Date)
  T.max <- max(all$Date)
  for (ty in c("BO", "VC")) {
    v <- all[all$type == ty & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
    U <- prepare.units(merge(pad.all.units(v, "first.cf", T.max), q, by = "Date", all.x = TRUE), c("RF", "MKT"))$units
    err <- function(b, m, units) vapply(units, function(u) horizon.error(u$cf, exp(cumsum(log(1 + u$X %*% c(1, b)))), m), 0)
    Q <- function(b, m, units) { e <- err(b, m, units); if (any(!is.finite(e))) 1e10 else sum(e^2) }
    fit <- function(m, units) optimize(Q, c(-1, 4), m = m, units = units, tol = 1e-10)$minimum
    for (m in c(0, 60, 120, 180, 240, 300, 360)) {
      b <- fit(m, U); e <- err(b, m, U); sh <- e^2 / sum(e^2)
      d <- 1e-6; g <- (err(b + d, m, U)^2 - err(b - d, m, U)^2) / (2 * d) # score contributions at the estimate
      bj <- vapply(seq_along(U), function(i) fit(m, U[-i]), 0); V <- length(U)
      res[[length(res) + 1]] <- data.frame(Type = ty, w = w, H = m, V = V, beta = b,
        se.jk = sqrt((V - 1) / V * sum((bj - mean(bj))^2)), neff.obj = 1 / sum(sh^2), neff.score = sum(g^2)^2 / sum(g^4),
        top.vintage = sapply(U, `[[`, "vintage")[which.max(sh)], top.share = max(sh), jk.range = diff(range(bj)))
    }
  }
}
r <- do.call(rbind, res)
g1 <- do.call(rbind, lapply(c("BO", "VC"), function(ty) { d <- read.csv(here(sprintf("results/figures_bovc_us_g1/empirical_%s_US_MKT_data.csv", ty))); d <- d[d$method == "Asymptotic" & d$Factor == "MKT", ]
  data.frame(Type = d$Type, w = d$weighting, H = d$max.month, beta.g1 = d$MKT, se.bartlett = d$SE.MKT, se.indep = d$SE.MKT.indep) }))
r <- merge(r, g1, by = c("Type", "w", "H"))
cat(sprintf("check: max |beta - G1 table| = %.1e\n", max(abs(r$beta - r$beta.g1))))
r$ratio.jk.bartlett <- r$se.jk / r$se.bartlett
x <- r[order(r$Type, r$w, r$H), c("Type", "w", "H", "V", "beta", "se.bartlett", "se.indep", "se.jk", "ratio.jk.bartlett", "neff.obj", "neff.score", "top.vintage", "top.share")]
num <- sapply(x, is.numeric) & !names(x) %in% c("H", "V", "top.vintage"); x[num] <- lapply(x[num], round, 3); print(x, row.names = FALSE)
saveRDS(r, commandArgs(TRUE)[1])

# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/jackknife_vc_fw_eg.R" 
# Leave-one-vintage-out refits of VC FW MKT+EG at H = 0, as jackknife.se() does them
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R"))))); source(here("helper/objective.R"))
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
all <- read.csv(here("empirical/data_prepared_2026_c/preqin_cashflows_FW_VYP_North America_NAV.csv")); all$Date <- as.Date(all$Date)
v <- all[all$type == "VC" & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
U <- prepare.units(merge(pad.all.units(v, "first.cf", max(all$Date)), q, by = "Date", all.x = TRUE), c("RF", "MKT", "EG"))$units
crit <- function(p, U) { e <- vapply(U, function(u) unit.horizon.error(u, 0, c(RF = 1, MKT = p[1], EG = p[2])), 0); if (any(!is.finite(e))) 1e10 else sqrt(sum(e^2)) / length(e) }
p.hat <- c(0.1295389, 3.3057900); n <- length(U)
est <- t(sapply(seq_len(n), function(i) optim(p.hat, crit, U = U[-i], method = "Nelder-Mead")$par))
se <- sqrt((n - 1) / n * colSums(sweep(est, 2, colMeans(est))^2))
cat(sprintf("n %d | jackknife SE MKT %.4f EG %.4f (row: 2.6340, 9.7797)\n", n, se[1], se[2]))
vint <- sapply(U, `[[`, "vintage"); o <- order(-abs(est[, 2] - p.hat[2]))
print(data.frame(left.out = vint[o], MKT = round(est[o, 1], 3), EG = round(est[o, 2], 3))[1:4, ], row.names = FALSE)
i <- o[1]; s2 <- est[-i, ]; se2 <- sqrt((n - 2) / (n - 1) * colSums(sweep(s2, 2, colMeans(s2))^2))
cat(sprintf("without that one refit: jackknife SE MKT %.4f EG %.4f\n", se2[1], se2[2]))

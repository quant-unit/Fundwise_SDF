# Check of the external reviews (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/return_floor_count.R" <stall.rds from stall_check.R>
# Reviewer claim: 23 of 260 reported Multi-table points have a monthly gross benchmark return
# below 0.5 somewhere in the estimation period (all VC, October 1987); minimum 0.3945.
res <- readRDS(commandArgs(TRUE)[1])
q <- read.csv("empirical/data_prepared_2026/q_factors.csv"); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1
span <- list()
for (w in c("EW", "FW")) { a <- read.csv(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_VYP_North America_NAV.csv", w)); a$Date <- as.Date(a$Date)
  for (ty in c("BO", "VC")) span[[paste(w, ty)]] <- c(min(a$Date[a$type == ty & a$Vintage <= 2010]), max(a$Date)) }
res$minG <- NA; res$month <- NA
for (k in seq_len(nrow(res))) { x <- res[k, ]; sp <- span[[paste(x$w, x$Type)]]; s <- q[q$Date >= sp[1] & q$Date <= sp[2], ]
  g <- 1 + s$RF + x$MKT.row * s$MKT + (if (x$Factor == "Alpha") x$Coef.row else x$Coef.row * s[[x$Factor]])
  res$minG[k] <- min(g); res$month[k] <- format(s$Date[which.min(g)], "%Y-%m") }
cat("window BO EW:", format(span[["EW BO"]]), "| VC EW:", format(span[["EW VC"]]), "\n")
b <- res[res$minG < 0.5, ]
cat(sprintf("points with min gross return < 0.5: %d of %d (excluding H=1 duplicates: %d); types: %s; months: %s; overall min %.4f\n",
  nrow(b), nrow(res), sum(b$H != 1), paste(unique(b$Type), collapse = ","), paste(unique(b$month), collapse = ","), min(res$minG)))
print(table(b$Factor)); cat("among the 6 edge rows:", sum(b$rel.gap > 1e-6), "\n")
b <- b[order(b$minG), ]; print(data.frame(b$Type, b$w, b$Factor, H = b$H / 12, minG = round(b$minG, 4)), row.names = FALSE)

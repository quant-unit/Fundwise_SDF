# r5 implementation check (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r5/compare_r4_r5.R" <r5 cache folder of EW North America q-factors>
# r5 (scratch run, EW North America, q-factors) vs r4 rows: what moved, and is it expected?
suppressMessages(source("helper/result_cache.R")); invisible(capture.output(suppressMessages(source("helper/getNPVs.R"))))
d.new <- commandArgs(TRUE)[1]; d.old <- "results/data_out_2026-emp-c/cache_q_factors_preqin_EW_VYP_North America"
rd <- function(d, ver) { x <- Reduce(rbind.all.columns, lapply(Sys.glob(file.path(d, "*_cached_res.csv")), read.csv)); dedupe.result.rows(x[x$pipeline.version %in% ver, ]) }
r5 <- rd(d.new, PIPELINE.VERSION); r4 <- rd(d.old, "2026-11-16-r3"); r4 <- r4[as.character(r4$kernel.bandwidth) == "llsw", ]
k <- c("Type", "Factor", "max.month", "CV.key")
m <- merge(r5, r4, by = k, suffixes = c(".n", ".o"))
cat("rows r5", nrow(r5), "| r4", nrow(r4), "| matched", nrow(m), "\n")
pc <- c("MKT", "ME", "IA", "ROE", "EG", "Alpha")
m$dmax <- apply(sapply(pc, function(v) abs(m[[paste0(v, ".n")]] - m[[paste0(v, ".o")]])), 1, max, na.rm = TRUE)
m$dval <- (m$value.n - m$value.o) / m$value.o
cat("objective: r5 lower than r4 by > 1e-6 rel:", sum(m$dval < -1e-6), "| higher by > 1e-6 rel:", sum(m$dval > 1e-6), "| max higher", signif(max(m$dval), 3), "\n")
cat("parameter change: <= 1e-3:", sum(m$dmax <= 1e-3), "| > 1e-3:", sum(m$dmax > 1e-3), "\n")
cat("r5 rows on a bound:", sum(r5$at.bound %in% TRUE), "| by constraint:\n"); print(table(r5$active.constraints[r5$at.bound %in% TRUE], r5$CV.key[r5$at.bound %in% TRUE] == "ALL", dnn = c("active", "full sample")))
big <- m[m$dmax > 1e-3, ]
big <- big[order(-big$dmax), c(k, "MKT.o", "MKT.n", "value.o", "value.n", "dval", "dmax", "at.bound.n", "active.constraints")]
cat("\nlargest moves (full sample and folds):\n"); print(head(big, 25), row.names = FALSE)
cat("\nmoves by CV.key == ALL:", sum(big$CV.key == "ALL"), "; types:", paste(names(table(big$Type)), table(big$Type), collapse = ", "), "\n")
a <- r5[r5$CV.key == "ALL", ]
cat("\nfull-sample rows: Wald.stat present", sum(is.finite(a$Wald.stat.MKT_1)), "of", nrow(a), "| jk rows", sum(is.finite(a$jk.selfcheck)), "| jk.n.bound > 0:", sum(a$jk.n.bound > 0, na.rm = TRUE), "| jk.n.switch > 0:", sum(a$jk.n.switch > 0, na.rm = TRUE), "| max selfcheck", signif(max(a$jk.selfcheck, na.rm = TRUE), 3), "\n")
cat("git.dirty in r5 rows:", paste(unique(r5$git.dirty), collapse = ","), "| gross.floor:", paste(unique(r5$gross.floor), collapse = ","), "\n")
x <- a[a$Type %in% c("BO", "VC") & (a$jk.n.bound > 0 | a$jk.n.switch > 0) & is.finite(a$jk.n.bound), c("Type", "Factor", "max.month", "jk.n.bound", "jk.n.switch", "at.bound")]
cat("\nBO/VC full-sample rows with flagged jackknife:\n"); print(x, row.names = FALSE)

# Was the r4 point inside Theta? (domain rebuilt from the same input and factors as the estimator)
source("helper/admissible.R")
q <- read.csv("empirical/data_prepared_2026/q_factors.csv"); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1; q$Alpha <- 1
inp <- read.csv("empirical/data_prepared_2026_c/preqin_cashflows_EW_VYP_North America_NAV.csv"); inp$Date <- as.Date(inp$Date); Tmax <- max(inp$Date)
inp <- inp[inp$Vintage <= 2010, ]
months.of <- function(ty) { s <- q[q$Date >= min(inp$Date[inp$type == ty]) & q$Date <= Tmax, ]; s }
m$feasible.o <- NA
for (i in which(m$dval > 1e-6)) {
  nm <- if (m$Factor[i] == "MKT") "MKT" else c("MKT", m$Factor[i])
  p <- setNames(vapply(nm, function(v) m[[paste0(v, ".o")]][i], 0), nm)
  m$feasible.o[i] <- in.domain(p, admissible.constraints(months.of(m$Type[i]), nm))
}
h <- m[m$dval > 1e-6, ]
cat("\nrows with a higher r5 objective:", nrow(h), "| r4 point inside Theta:", sum(h$feasible.o), "\n")
if (any(h$feasible.o)) print(h[h$feasible.o, c(k, "MKT.o", "MKT.n", "value.o", "value.n", "dval", "active.constraints")], row.names = FALSE)
cat("\nfull-sample rows that moved (> 1e-3) or are on a bound:\n")
fa <- m[m$CV.key == "ALL" & (m$dmax > 1e-3 | m$at.bound.n %in% TRUE), ]
print(fa[order(fa$Type, fa$Factor, fa$max.month), c("Type", "Factor", "max.month", "MKT.o", "MKT.n", "value.o", "value.n", "dmax", "active.constraints", "min.G")], row.names = FALSE)

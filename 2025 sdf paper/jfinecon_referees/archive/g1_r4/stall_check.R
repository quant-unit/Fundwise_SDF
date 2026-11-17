# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/stall_check.R" <stall.rds>
# Are the full-sample two-parameter estimates in the G1 Multi tables minima?
# Exact pipeline criterion sqrt(sum e^2)/n; checked by reproducing each row's `value`.
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R")))))
suppressMessages(source(here("helper/result_cache.R"))); source(here("helper/objective.R")); source(here("helper/parallel_workers.R"))
tab <- do.call(rbind, lapply(c("BO", "VC"), function(ty) read.csv(here(sprintf("results/figures_bovc_us_g1/empirical_%s_US_Multi_data.csv", ty)))))
tab <- tab[tab$method == "Asymptotic" & tab$Factor != "MKT", ]
dirs <- list.dirs("results/data_out_2026-emp-c", recursive = FALSE)
cache <- Reduce(rbind.all.columns, lapply(dirs[grepl("q_factors_preqin_(alpha_ALL_)?(EW|FW)_VYP_North America$", dirs)], function(d) { x <- read.result.cache(d, verbose = FALSE); x$folder <- basename(d); x }))
cache <- cache[cache$CV.key == "ALL" & cache$lambda == 0, ]
m <- merge(tab[, c("Type", "Factor", "max.month", "weighting", "cache.identity", "MKT", "Coef", "hac.V", "jk.selfcheck")],
  cache[, c("Type", "Factor", "max.month", "cache.identity", "value", "X", "folder", "Alpha", "ME", "IA", "ROE", "EG")], by = c("Type", "Factor", "max.month", "cache.identity"))
cat("table rows", nrow(tab), "matched to cache", nrow(m), "\n")
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1; q$Alpha <- 1
UU <- list()
for (w in c("EW", "FW")) {
  all <- read.csv(here(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_VYP_North America_NAV.csv", w))); all$Date <- as.Date(all$Date); T.max <- max(all$Date)
  for (ty in c("BO", "VC")) {
    v <- all[all$type == ty & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
    UU[[paste(w, ty)]] <- prepare.units(merge(pad.all.units(v, "first.cf", T.max), q, by = "Date", all.x = TRUE), c("RF", "MKT", "ME", "IA", "ROE", "EG", "Alpha"))$units
  }
}
crit <- function(p, nm, H, U) { e <- vapply(U, function(u) unit.horizon.error(u, H, c(RF = 1, setNames(p, nm))), 0); if (any(!is.finite(e))) 1e10 else sqrt(sum(e^2)) / length(e) }
one <- function(k) {
  x <- m[k, ]; U <- UU[[paste(x$weighting, x$Type)]]; nm <- c("MKT", x$Factor); p0 <- c(x$MKT, x$Coef)
  c.rep <- crit(p0, nm, x$max.month, U)
  set.seed(k); st <- rbind(p0, c(1, 0), cbind(runif(30, -1, 4), if (x$Factor == "Alpha") runif(30, -0.02, 0.02) else runif(30, -4, 4)))
  fits <- t(apply(st, 1, function(s) {
    o <- if (x$Factor == "Alpha") optim(s, crit, nm = nm, H = x$max.month, U = U, method = "L-BFGS-B", lower = c(-Inf, -0.02), upper = c(Inf, 0.02), control = list(parscale = c(1, 0.001), factr = 1e3))
         else optim(s, crit, nm = nm, H = x$max.month, U = U, method = "Nelder-Mead", control = list(reltol = 1e-12, maxit = 5000))
    c(o$par, o$value) }))
  b <- fits[which.min(fits[, 3]), ]
  data.frame(Type = x$Type, w = x$weighting, Factor = x$Factor, H = x$max.month, V = length(U), V.row = x$hac.V,
    value.row = x$value, value.recomp = c.rep, value.best = b[3], rel.gap = (x$value - b[3]) / x$value,
    MKT.row = x$MKT, MKT.best = b[1], Coef.row = x$Coef, Coef.best = b[2], jk.selfcheck = x$jk.selfcheck)
}
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(m)), one, mc.cores = n.parallel.workers(nrow(m))))
cat(sprintf("units match hac.V: %s | max |recomputed - row value| / value: %.1e\n", all(res$V == res$V.row), max(abs(res$value.recomp - res$value.row) / res$value.row)))
cat(sprintf("rows where best found is lower by > 1e-6 (relative): %d of %d\n", sum(res$rel.gap > 1e-6), nrow(res)))
res <- res[order(-res$rel.gap), ]
num <- sapply(res, is.numeric); r2 <- res; r2[num] <- lapply(r2[num], signif, 4)
print(head(r2[, c("Type","w","Factor","H","value.row","value.best","rel.gap","MKT.row","MKT.best","Coef.row","Coef.best","jk.selfcheck")], 30), row.names = FALSE)
saveRDS(res, commandArgs(TRUE)[1])

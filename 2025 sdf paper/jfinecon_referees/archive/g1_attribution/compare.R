# G1 attribution table from the cells written by reproduce.sh.
# Run from r_project/:  Rscript <this file> <OUT>   (OUT = reproduce.sh output root)
# Cells: oldB/oldC = code d2bdf6e on data_prepared_2026_b/_c; newB = current code on _b;
# ownB/ownC = current code with the old window endpoint (each unit ends at its own last
# cash flow) on _b/_c; newC = the G1 results in results/data_out_2026-emp-c.
suppressMessages(source("helper/result_cache.R"))
out <- commandArgs(TRUE)[1]
keep <- c("Type", "Factor", "max.month", "MKT", "ME", "IA", "ROE", "EG", "Alpha")
coef.of <- function(d) { d$Coef <- mapply(function(f, i) d[[f]][i], d$Factor, seq_len(nrow(d))); d }
load.raw <- function(cell) do.call(rbind, lapply(c("EW", "FW"), function(w) {
  f <- Sys.glob(file.path(out, cell, sprintf("cache_q_factors_preqin_%s_VYP_North America", w), "*.csv"))
  if (!length(f)) stop("no results in cell ", cell, " (", w, "): did the run fail?")
  d <- do.call(rbind, lapply(f, function(x) { r <- read.csv(x); if ("CV.key" %in% names(r)) r <- r[r$CV.key == "ALL", ]; r[, keep] }))
  stopifnot(!anyDuplicated(d[, c("Type", "Factor", "max.month")])) # one fresh run per cell
  d$w <- w; coef.of(d) }))
load.g1 <- function() do.call(rbind, lapply(c("EW", "FW"), function(w) {
  d <- read.result.cache(file.path("results/data_out_2026-emp-c", sprintf("cache_q_factors_preqin_%s_VYP_North America", w)), verbose = FALSE)
  d <- d[d$CV.key == "ALL", keep]; stopifnot(!anyDuplicated(d[, c("Type", "Factor", "max.month")])); d$w <- w; coef.of(d) }))
cells <- list(oldB = load.raw("old_b"), oldC = load.raw("old_c"), newB = load.raw("new_b"), newC = load.g1(),
              ownB = load.raw("ownend_b"), ownC = load.raw("ownend_c"))
paper <- unique(do.call(rbind, lapply(c("BO", "VC"), function(ty) do.call(rbind, lapply(c("MKT", "Alpha", "Multi"), function(tab) {
  d <- read.csv(sprintf("../2025 sdf paper/tex/empirical_bovc_us/empirical_%s_US_%s_data.csv", ty, tab)); d <- d[d$method == "Asymptotic", ]
  data.frame(Type = d$Type, Factor = d$Factor, max.month = d$max.month, w = d$weighting, MKT = d$MKT, Coef = d$Coef) })))))
k <- c("Type", "Factor", "max.month", "w")
tab <- paper[, c(k, "MKT", "Coef")]; names(tab)[5:6] <- c("MKT.paper", "Coef.paper")
for (n in names(cells)) { x <- cells[[n]][, c(k, "MKT", "Coef")]; names(x)[5:6] <- paste0(c("MKT.", "Coef."), n); tab <- merge(tab, x, by = k, all = TRUE) }
tab <- tab[tab$Type %in% c("BO", "VC") & tab$max.month %in% c(0, 60, 120, 180, 240, 300, 360), ]
tab <- tab[order(tab$Type, tab$Factor, tab$w, tab$max.month), ]
num <- sapply(tab, is.numeric) & names(tab) != "max.month"; tab[num] <- lapply(tab[num], round, 5)
write.csv(tab, file.path(out, "g1_attribution_north_america.csv"), row.names = FALSE)
p <- !is.na(tab$MKT.paper) & !is.na(tab$MKT.oldB)
cat(sprintf("rows %d; paper vs oldB: %d rows, max |diff| %.1e\n", nrow(tab), sum(p), max(abs(c(tab$MKT.paper - tab$MKT.oldB, tab$Coef.paper - tab$Coef.oldB)[c(p, p)]))))

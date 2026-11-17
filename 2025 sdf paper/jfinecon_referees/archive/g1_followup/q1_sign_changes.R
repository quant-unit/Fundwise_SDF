# G1 follow-up (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_followup/q1_sign_changes.R" <OUT of archive/g1_attribution/reproduce.sh>
# Q1: which step flips the sign of multi-factor loadings (paper -> G1)?
suppressMessages(source("helper/result_cache.R"))
out <- commandArgs(TRUE)[1]; facs <- c("ME", "IA", "ROE", "EG")
raw <- function(cell, w) do.call(rbind, lapply(Sys.glob(file.path(out, cell, sprintf("cache_q_factors_preqin_%s_VYP_North America", w), "*.csv")),
  function(f) { r <- read.csv(f); if ("CV.key" %in% names(r)) r <- r[r$CV.key == "ALL", ]; r }))
pick <- function(d, w) { d <- d[d$Type %in% c("BO", "VC") & d$Factor %in% facs, ]
  data.frame(Type = d$Type, Factor = d$Factor, H = d$max.month, w = w, coef = mapply(function(f, i) d[[f]][i], d$Factor, seq_len(nrow(d)))) }
cell <- function(name) do.call(rbind, lapply(c("EW", "FW"), function(w) pick(raw(name, w), w)))
g1 <- do.call(rbind, lapply(c("EW", "FW"), function(w) { d <- read.result.cache(file.path("results/data_out_2026-emp-c", sprintf("cache_q_factors_preqin_%s_VYP_North America", w)), verbose = FALSE); pick(d[d$CV.key == "ALL", ], w) }))
k <- c("Type", "Factor", "H", "w")
m <- Reduce(function(a, b) merge(a, b, by = k), list(
  setNames(cell("old_b"), c(k, "oldB")), setNames(cell("old_c"), c(k, "oldC")), setNames(cell("ownend_b"), c(k, "ownB")),
  setNames(cell("ownend_c"), c(k, "ownC")), setNames(g1, c(k, "newC"))))
pap <- do.call(rbind, lapply(c("BO", "VC"), function(ty) { d <- read.csv(sprintf("../2025 sdf paper/tex/empirical_bovc_us/empirical_%s_US_Multi_data.csv", ty)); d <- d[d$method == "Asymptotic" & d$Factor %in% facs, ]
  data.frame(Type = d$Type, Factor = d$Factor, H = d$max.month, w = d$weighting, paper = d$Coef) }))
tab <- do.call(rbind, lapply(c("BO", "VC"), function(ty) { d <- read.csv(sprintf("results/figures_bovc_us_g1/empirical_%s_US_Multi_data.csv", ty)); d <- d[d$method == "Asymptotic" & d$Factor %in% facs, ]
  data.frame(Type = d$Type, Factor = d$Factor, H = d$max.month, w = d$weighting, g1tab = d$Coef, se = d$SE.Coef) }))
m <- merge(merge(m, pap, by = k), tab, by = k)
cat(sprintf("rows %d | checks: max|paper - oldB| %.1e, max|G1 table - newC| %.1e\n", nrow(m), max(abs(m$paper - m$oldB)), max(abs(m$g1tab - m$newC))))
flip <- function(a, b) sign(a) != sign(b)
m$flip <- flip(m$paper, m$newC); m$t <- m$newC / m$se; m$sig <- abs(m$t) > 1.96
# order A: code (B2, optimizer) -> data -> endpoint;  order B: data -> code -> endpoint
m$A.code <- flip(m$oldB, m$ownB); m$A.data <- flip(m$ownB, m$ownC); m$A.end <- flip(m$ownC, m$newC)
m$B.data <- flip(m$oldB, m$oldC); m$B.code <- flip(m$oldC, m$ownC)
for (ty in c("BO", "VC")) { x <- m[m$Type == ty & m$flip, ]
  cat(sprintf("\n%s: %d sign changes of %d (H = 0 and 1 both counted); significant in G1 (|t| > 1.96): %d\n", ty, nrow(x), sum(m$Type == ty), sum(x$sig)))
  cat(sprintf("  order code->data->endpoint: code %d, data %d, endpoint %d (a row can flip in more than one step)\n", sum(x$A.code), sum(x$A.data), sum(x$A.end)))
  cat(sprintf("  order data->code->endpoint: data %d, code %d, endpoint %d\n", sum(x$B.data), sum(x$B.code), sum(x$A.end)))
  x$step <- ifelse(x$A.end & !x$A.data & !x$A.code, "endpoint", ifelse(x$A.data & !x$A.end & !x$A.code, "data", ifelse(x$A.code & !x$A.data & !x$A.end, "code", "several")))
  print(table(x$step, x$H >= 240, dnn = c("step", "H >= 20y")))
  x <- x[order(x$Factor, x$w, x$H), c("Factor", "w", "H", "paper", "oldC", "ownC", "newC", "se", "t", "step")]; x[, 4:9] <- round(x[, 4:9], 3); print(x, row.names = FALSE) }
saveRDS(m, file.path(out, "q1_signs.rds"))

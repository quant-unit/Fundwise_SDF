# Fixed-b critical values of the paper (helper/fixedb.R) against the cubic approximation of
# Kiefer and Vogelsang (2005, Econometric Theory 21, Table 1; Bartlett kernel, 97.5% quantile):
#   cv(b) = 1.9600 + 2.9694 b + 0.4160 b^2 - 0.5324 b^3.
# Part 1 (public, seconds): our simulated grid b = 0.02, ..., 1 against the KV05 cubic and against
#   a cubic refitted to our grid with the intercept fixed at 1.96, as in KV05.
# Part 2 (needs the licensed empirical results; skipped if absent): Buyout and VC U_{v,H} tests
#   (beta = 0, beta = 1, second coefficient = 0; every model and horizon) whose |t| lies between
#   the paper's 5% critical value and the KV05 cubic value, i.e. whose 5% decision would change.
# Run from r_project/:  LANG=en_US.UTF-8 Rscript results/fixedb_kv05_comparison.R
suppressPackageStartupMessages(library(here))
source(here("helper/fixedb.R"))
kv05 <- function(b) 1.9600 + 2.9694 * b + 0.4160 * b^2 - 0.5324 * b^3

# ---- Part 1: critical values ----
grid <- seq(0.02, 1, by = 0.02)
cv.grid <- vapply(grid, function(b) fixedb.crit.t.table(0.05, b), 0) # quantile table, 1e5 draws, T = 1000
fit <- lm(I(cv.grid - 1.96) ~ 0 + grid + I(grid^2) + I(grid^3))
own <- function(b) 1.96 + sum(coef(fit) * c(b, b^2, b^3))
cat(sprintf("Cubic refitted to our grid: %.4f b %+.4f b^2 %+.4f b^3; max |residual| %.3f\n",
  coef(fit)[1], coef(fit)[2], coef(fit)[3], max(abs(resid(fit)))))
gap <- kv05(grid) - cv.grid
cat(sprintf("KV05 cubic minus our grid: mean %.3f, range %.3f to %.3f\n", mean(gap), min(gap), max(gap)))
B <- c(`7/25` = 7 / 25, `7/23` = 7 / 23, `7/27` = 7 / 27, `7/26` = 7 / 26, `7/29` = 7 / 29)
tab <- data.frame(b = names(B),
  conditional = vapply(B, function(b) fixedb.crit.t(0.05, b), 0),       # values used in the paper
  quantile.table = vapply(B, function(b) fixedb.crit.t.table(0.05, b), 0),
  own.cubic = vapply(B, own, 0), kv05.cubic = kv05(B))
tab$kv05.minus.paper <- tab$kv05.cubic - tab$conditional
print(format(tab, digits = 4), row.names = FALSE)

# ---- Part 2: 5% decisions that would change (licensed results) ----
dirs <- c("results/data_out_2026-emp-c", "results/data_out_2026-emp-c-allfunds", "results/data_out_2026-emp-c-navtrim20")
dirs <- dirs[dir.exists(dirs)]
if (!length(dirs)) { cat("Part 2 skipped: no empirical result folders\n"); quit(save = "no") }
u.dirs <- list.dirs(dirs, recursive = FALSE)
u.dirs <- u.dirs[grepl("U_pair$", u.dirs)]
hits <- list(); n.tests <- 0L
for (d in u.dirs) {
  f <- file.path(d, "0_asymptotic_inference_summary.csv")
  if (!file.exists(f)) next
  x <- read.csv(f, stringsAsFactors = FALSE)
  x <- x[x$Type %in% c("BO", "VC") & !x$at.bound & is.finite(x$SE.MKT) & x$SE.MKT > 0, ]
  for (i in seq_len(nrow(x))) {
    c.paper <- fixedb.crit.t(0.05, x$hac.b[i]); c.kv <- kv05(x$hac.b[i])
    t <- c(beta.0 = x$t.MKT[i], beta.1 = (x$MKT[i] - 1) / x$SE.MKT[i], coef.0 = x$t.Coef[i])
    for (k in names(t)[is.finite(t)]) {
      n.tests <- n.tests + 1L
      if (abs(t[[k]]) >= c.paper && abs(t[[k]]) < c.kv)
        hits[[length(hits) + 1]] <- data.frame(folder = d, Type = x$Type[i], H = x$max.month[i] / 12,
          model = x$Factor[i], test = k, t = t[[k]], b = x$hac.b[i], cv.paper = c.paper, cv.kv05 = c.kv)
    }
  }
}
cat(sprintf("Part 2: %d Buyout/VC U_{v,H} t-tests in %d result folders\n", n.tests, length(u.dirs)))
if (length(hits)) print(format(do.call(rbind, hits), digits = 4), row.names = FALSE) else cat("no 5% decision would change\n")

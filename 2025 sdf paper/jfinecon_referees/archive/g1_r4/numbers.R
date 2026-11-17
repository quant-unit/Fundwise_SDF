# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/numbers.R" 
# Numbers for the State-of-evidence update (r4: S = ceil(1.3 sqrt V), fixed-b)
source("helper/fixedb.R")
rd <- function(ty, f = "MKT") { d <- read.csv(sprintf("results/figures_bovc_us_g1/empirical_%s_US_%s_data.csv", ty, f)); d[d$method == "Asymptotic" & d$max.month != 1, ] }
old <- function(ty) { d <- read.csv(text = system(sprintf("git show 90e5421:r_project/results/figures_bovc_us_g1/empirical_%s_US_MKT_data.csv", ty), intern = TRUE)); d[d$method == "Asymptotic" & d$max.month != 1, ] }
H7 <- c(0, 60, 120, 180, 240, 300, 360)
cat("== single-factor L2: SE r3 -> r4, jackknife, N_eff, p(beta=0), Wald p(beta=1) fixed-b [chi2]\n")
for (ty in c("BO", "VC")) { d <- rd(ty); o <- old(ty)
  d <- merge(d[d$Factor == "MKT", ], o[o$Factor == "MKT", c("weighting", "max.month", "SE.MKT")], by = c("weighting", "max.month"), suffixes = c("", ".r3"))
  d <- d[d$max.month %in% H7, ]; d <- d[order(d$weighting, d$max.month), ]
  print(data.frame(ty, w = d$weighting, H = d$max.month / 12, beta = round(d$MKT, 3), SE.r3 = round(d$SE.MKT.r3, 3), SE.A = round(d$SE.MKT, 3), SE.ind = round(d$SE.MKT.indep, 3), SE.jk = round(d$SE.MKT.jk, 3),
    jk.A = round(d$SE.MKT.jk / d$SE.MKT, 2), Neff = round(d$neff.units, 2), S = d$hac.S, V = d$hac.V, p0 = signif(d$p.MKT, 2), pW1 = signif(d$Wald.p.value.MKT_1, 2), pW1.chi = signif(d$Wald.p.value.MKT_1.chisq, 2)), row.names = FALSE) }
cat("\n== ratio SE_JK / SE_A, single-factor L2, all 13 horizons (H=1 excluded):\n")
for (ty in c("BO", "VC")) { d <- rd(ty); d <- d[d$Factor == "MKT", ]; for (w in c("EW", "FW")) { x <- d[d$weighting == w, ]; r <- x$SE.MKT.jk / x$SE.MKT
  cat(sprintf("  %s %s: %.2f-%.2f | Wald beta=1 rejected at 5%% (fixed-b): %d of %d, H = %s\n", ty, w, min(r), max(r), sum(x$Wald.p.value.MKT_1 < 0.05), nrow(x), paste(x$max.month[x$Wald.p.value.MKT_1 < 0.05] / 12, collapse = ","))) } }
cat("\n== (alpha, beta) model L2 (Coef = Alpha per month)\n")
for (ty in c("BO", "VC")) { d <- rd(ty); d <- d[d$Factor == "Alpha" & d$max.month %in% c(H7, 270, 330), ]; d <- d[order(d$weighting, d$max.month), ]
  print(data.frame(ty, w = d$weighting, H = d$max.month / 12, beta = round(d$MKT, 3), SE.b = round(d$SE.MKT, 3), alpha.pct = round(100 * d$Coef, 3), SE.a.pct = round(100 * d$SE.Coef, 3), SE.a.jk = round(100 * d$SE.Coef.jk, 3), p.alpha = signif(d$p.Coef, 2), bound = d$at.bound), row.names = FALSE) }
cat("\n== multi-factor loadings that changed sign vs the paper: |t| and fixed-b p under r4\n")
facs <- c("ME", "IA", "ROE", "EG")
pap <- do.call(rbind, lapply(c("BO", "VC"), function(ty) { d <- read.csv(sprintf("../2025 sdf paper/tex/empirical_bovc_us/empirical_%s_US_Multi_data.csv", ty)); d <- d[d$method == "Asymptotic" & d$Factor %in% facs, ]
  data.frame(Type = d$Type, Factor = d$Factor, H = d$max.month, w = d$weighting, paper = d$Coef) }))
new <- do.call(rbind, lapply(c("BO", "VC"), function(ty) { d <- read.csv(sprintf("results/figures_bovc_us_g1/empirical_%s_US_Multi_data.csv", ty)); d <- d[d$method == "Asymptotic" & d$Factor %in% facs, ]
  data.frame(Type = d$Type, Factor = d$Factor, H = d$max.month, w = d$weighting, coef = d$Coef, t = d$Coef / d$SE.Coef, p = d$p.Coef) }))
m <- merge(pap, new, by = c("Type", "Factor", "H", "w")); f <- m[sign(m$paper) != sign(m$coef), ]
cat(sprintf("rows %d, sign changes %d (BO %d, VC %d); max |t| %.2f; min fixed-b p %.3f; any p < 0.05: %s\n", nrow(m), nrow(f), sum(f$Type == "BO"), sum(f$Type == "VC"), max(abs(f$t)), min(f$p), any(f$p < 0.05)))
cat("\n== VC EW ME at 20 years (D example)\n"); x <- rd("VC", "Multi"); print(x[x$Factor == "ME" & x$weighting == "EW" & x$max.month == 240, c("Coef", "SE.Coef", "SE.Coef.indep", "SE.Coef.jk", "p.Coef")])
cat("\n== fixed-b critical values used: BO b =", round(7/25, 3), "5%:", round(fixedb.crit.t(0.05, 7/25), 3), "| VC b =", round(7/27, 3), "5%:", round(fixedb.crit.t(0.05, 7/27), 3), "\n")

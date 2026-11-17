# G1 follow-up (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_followup/q3_alpha_beta_ridge.R" <output folder for the contour PDF>
# Q3: how well are alpha and beta separated? Buyout, VYP, North America, V <= 2010,
# G1 criterion Q = sqrt(sum e^2) / n with G_u = 1 + RF + beta MKT + alpha.
suppressPackageStartupMessages(library(here))
invisible(capture.output(suppressMessages(source(here("helper/getNPVs.R")))))
source(here("helper/objective.R"))
q <- read.csv(here("empirical/data_prepared_2026/q_factors.csv")); q$Date <- as.Date(format(as.Date(q$Date) + 1, "%Y-%m-01")) - 1; q$Alpha <- 1
outdir <- commandArgs(TRUE)[1]; res <- list()
pdf(file.path(outdir, "bo_alpha_beta_contours.pdf"), width = 11, height = 7.5); par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (w in c("EW", "FW")) {
  all <- read.csv(here(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_VYP_North America_NAV.csv", w))); all$Date <- as.Date(all$Date)
  v <- all[all$type == "BO" & all$Vintage <= 2010, c("Fund.ID", "Date", "CF", "type", "Vintage")]; v$Fund.ID <- paste(v$Fund.ID, v$type, sep = "_")
  U <- prepare.units(merge(pad.all.units(v, "first.cf", max(all$Date)), q, by = "Date", all.x = TRUE), c("RF", "MKT", "Alpha"))$units
  g1 <- read.csv(here("results/figures_bovc_us_g1/empirical_BO_US_Alpha_data.csv")); g1 <- g1[g1$method == "Asymptotic" & g1$Factor == "Alpha" & g1$weighting == w, ]
  for (m in c(0, 120, 240)) {
    Q <- function(p) { e <- vapply(U, function(u) horizon.error(u$cf, exp(cumsum(log(1 + u$X %*% c(1, p[1], p[2])))), m), 0)
      if (any(!is.finite(e))) NA else sqrt(sum(e^2)) / length(e) }
    o <- optim(c(1, 0), function(p) { z <- Q(p); if (is.na(z)) 1e10 else z }, method = "L-BFGS-B", lower = c(-Inf, -0.02), upper = c(Inf, 0.02),
               control = list(parscale = c(1, 0.001), factr = 1e3))
    gb <- seq(-0.5, 2, length.out = 101); ga <- seq(-0.01, 0.02, length.out = 101)
    Z <- outer(gb, ga, Vectorize(function(b, a) Q(c(b, a))))
    prof <- sapply(gb, function(b) { z <- optimize(function(a) { y <- Q(c(b, a)); if (is.na(y)) 1e10 else y }, c(-0.02, 0.02), tol = 1e-10); z$objective })
    flat <- gb[prof <= o$value * 1.001]
    p1 <- optimize(function(a) { y <- Q(c(1, a)); if (is.na(y)) 1e10 else y }, c(-0.02, 0.02), tol = 1e-10) # beta = 1, alpha profiled
    # Hessian of Q in (beta, alpha in %/month) at the optimum -> correlation of the estimates
    h <- optimHess(o$par, function(p) Q(p), control = list(parscale = c(1, 0.001)))
    cr <- -h[1, 2] / sqrt(h[1, 1] * h[2, 2]) # corr of inverse Hessian (2x2)
    gg <- g1[g1$max.month == m, ]
    res[[length(res) + 1]] <- data.frame(w = w, H = m, beta = o$par[1], alpha = o$par[2], Q = o$value, beta.g1 = gg$MKT, alpha.g1 = gg$Coef,
      corr = cr, flat.lo = min(flat), flat.hi = max(flat), se.beta.g1 = gg$SE.MKT, Q.beta1.rel = p1$objective / o$value, alpha.at.beta1 = p1$minimum)
    contour(gb, ga * 100, Z / o$value, levels = c(1.001, 1.01, 1.05, 1.2, 1.5, 2), xlab = "beta", ylab = "alpha (% per month)",
            main = sprintf("Buyout %s, H = %d y: Q / min Q", w, m / 12))
    points(o$par[1], o$par[2] * 100, pch = 19); lines(gb, sapply(gb, function(b) 100 * optimize(function(a) { y <- Q(c(b, a)); if (is.na(y)) 1e10 else y }, c(-0.02, 0.02))$minimum), lty = 2)
  }
}
dev.off()
r <- do.call(rbind, res); cat(sprintf("check: max |beta - G1| %.4f, max |alpha - G1| %.5f\n", max(abs(r$beta - r$beta.g1)), max(abs(r$alpha - r$alpha.g1))))
r[, -(1:2)] <- signif(r[, -(1:2)], 4); print(r, row.names = FALSE)

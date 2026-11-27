# Specification tests of the U_{v,H} identifying conditions (Theorem T7) on the Preqin
# Buyout and VC panels (helper/horizon_invariance.R):
#   1. horizon-invariance test beta_0 = beta_H for the single-factor U_{v,H} (q_v = 1),
#      primary horizon H = 15 yr (the zero-bias horizon range of Proposition 2), with
#      H = 5, 10, 12.5, 20 yr reported as a profile;
#   2. two-moment shared-shock test at H = 0 (vintage mean pricing error and pair product).
# Run from r_project/ after empirical/prep_preqin.R:
#   LANG=en_US.UTF-8 Rscript results/horizon_invariance_test.R
# Writes results/figures_bovc_us_g1/horizon_invariance_test.csv and
# results/figures_bovc_us_g1/two_moment_test.csv. Needs the licensed Preqin data.
suppressMessages({
  source("simulation/mc_driver.R")
  if (!exists("fixedb.crit.t")) source("helper/fixedb.R")
  source("helper/horizon_invariance.R")
})
git.provenance <- function() list(
  commit = tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1], error = function(e) NA_character_),
  dirty = tryCatch(length(system("git status --porcelain --untracked-files=no", intern = TRUE, ignore.stderr = TRUE)) > 0, error = function(e) NA))

pub <- mc.public.factors()
fc <- unique(c("RF", setdiff(colnames(pub), "Date")))
load.fund <- function(w, type) {
  d <- read.csv(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s_North America_NAV.csv", w), stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date)
  d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  d <- pad.all.units(d, "first.cf", max(d$Date))
  d <- merge(d, pub, by = "Date", all.x = TRUE)
  d <- d[d$type == type & d$Vintage <= 2010, ]
  d$Fund.ID <- droplevels(as.factor(as.character(d$Fund.ID)))
  prepare.fast(d, fc)
}

H.PRIMARY <- 180
H.PROFILE <- c(60, 120, 150, 180, 240)
hi <- list(); tm <- list()
for (ty in c("BO", "VC")) for (w in c("EW", "FW")) {
  fp <- load.fund(w, ty)
  est <- if (w == "FW") MC.ESTIMATORS$U.eq.fw else MC.ESTIMATORS$U.eq
  f0 <- mc.fit(est, fp, "MKT", 0)
  for (m in H.PROFILE) {
    r <- u.horizon.contrast(est, fp, 0, m, f0 = f0)
    hi[[length(hi) + 1]] <- data.frame(type = ty, weighting = w, H = m / 12, primary = m == H.PRIMARY,
      se.fit0 = f0$SE.MKT, r, stringsAsFactors = FALSE)
  }
  tm[[length(tm) + 1]] <- data.frame(type = ty, weighting = w, V = sum(fp$group.n >= 2), u.two.moment.test(est, fp), stringsAsFactors = FALSE)
}
hi <- do.call(rbind, hi); tm <- do.call(rbind, tm)
# The H = 0 fit must reproduce the production estimate and SE (same code path)
stopifnot(all(abs(hi$se0 - hi$se.fit0) < 1e-8 | hi$bound))
prov <- git.provenance()
hi$git.commit <- prov$commit; hi$git.dirty <- prov$dirty
tm$git.commit <- prov$commit; tm$git.dirty <- prov$dirty
write.csv(hi, "results/figures_bovc_us_g1/horizon_invariance_test.csv", row.names = FALSE)
write.csv(tm, "results/figures_bovc_us_g1/two_moment_test.csv", row.names = FALSE)
cat("wrote results/figures_bovc_us_g1/horizon_invariance_test.csv and two_moment_test.csv\n")
print(format(hi[, c("type", "weighting", "H", "beta0", "beta1", "diff", "se0", "se1", "corr", "se.diff", "t", "p", "V", "b", "bound")], digits = 3), row.names = FALSE)
print(format(tm, digits = 3), row.names = FALSE)

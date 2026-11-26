# Step 4 Monte Carlo, design D0 (revision plan, Section I, "Step 4: scope and
# compute plan"): size of the fixed-b t-test in a Gaussian location model on
# the actual vintage calendars, with gaps and finite V. Run from r_project/:
#   LANG=en_US.UTF-8 Rscript simulation/mc_d0.R
# x_v = mu + u_v at the vintage years of the calendar; u is an AR(1) series in
# calendar years (rho in {0, 0.3, 0.6}, unit variance) sampled at those years,
# so gaps shorten the effective dependence. Test of mu = 0 at 5%: Bartlett
# sandwich over vintage-year distances with S = ceil(1.3 sqrt V) (the
# estimator's kernel), Kiefer-Vogelsang fixed-b critical value at b = S / V
# (conditional method, helper/fixedb.R) and, for comparison, the normal value.
# 10^5 draws per cell; one .rds in simulation/data_out_mc/d0/.
source("helper/inference.R") # bartlett.weight, hac.bandwidth.llsw, fixed-b functions
source("simulation/mc_generate.R")

d0.size <- function(years, rho, reps = 1e5, seed = 1, alpha = 0.05) {
  set.seed(seed)
  V <- length(years); S <- hac.bandwidth.llsw(V); b <- S / V
  W <- outer(years, years, bartlett.weight, bw = S)
  span <- min(years):max(years)
  # AR(1) with unit variance on the calendar span, sampled at the vintage years
  u <- matrix(rnorm(length(span) * reps), length(span))
  if (rho != 0) for (t in 2:length(span)) u[t, ] <- rho * u[t - 1, ] + sqrt(1 - rho^2) * u[t, ]
  x <- u[match(years, span), , drop = FALSE]
  m <- colMeans(x); e <- sweep(x, 2, m)
  se <- sqrt(colSums(e * (W %*% e)) / V^2)
  t <- m / se
  c.fb <- fixedb.crit.t(alpha, b)
  data.frame(V = V, S = S, b = b, rho = rho, crit.fixedb = c.fb, size.fixedb = mean(abs(t) > c.fb),
             size.normal = mean(abs(t) > qnorm(1 - alpha / 2)), mcse = sqrt(alpha * (1 - alpha) / reps), reps = reps)
}

shapes <- list(BO = mc.data.shape("BO"), VC = mc.data.shape("VC"))
cal <- list(
  "Buyout, all vintages (L2)" = as.integer(names(shapes$BO)),
  "Buyout, n_v >= 2 (U)" = as.integer(names(shapes$BO))[shapes$BO >= 2],
  "VC, all vintages (L2)" = as.integer(names(shapes$VC)),
  "VC, n_v >= 2 (U)" = as.integer(names(shapes$VC))[shapes$VC >= 2],
  "25 consecutive vintages (no gaps)" = 1986:2010
)
res <- do.call(rbind, lapply(names(cal), function(k) do.call(rbind, lapply(c(0, 0.3, 0.6), function(r) cbind(calendar = k, d0.size(cal[[k]], r))))))
dir.create("simulation/data_out_mc/d0", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(results = res, commit = system("git rev-parse --short HEAD", intern = TRUE),
             dirty = length(system("git status --porcelain --untracked-files=no", intern = TRUE)) > 0), "simulation/data_out_mc/d0/d0.rds")
print(res, digits = 3, row.names = FALSE)

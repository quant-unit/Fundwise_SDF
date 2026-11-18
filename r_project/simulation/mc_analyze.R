######################
### Step 4 Monte Carlo: analysis
### Revision plan, Section I, "Step 4: scope and compute plan"
######################
#
# mc.summary(): per estimator, model, parameter and horizon over the
# replications of a design: mean, bias against the true value, its Monte
# Carlo SE, sd, RMSE, share on a bound, and the size of the two-sided
# fixed-b t-test of the true value at 5% (helper/fixedb.R; estimates on a bound
# have no SE and count as not tested, their share is reported separately).
# Criterion targets theta_H^dagger and the error decomposition follow in the
# next stage.
if (!exists("fixedb.pvalue")) source("helper/fixedb.R")

mc.summary <- function(rows, theta0, level = 0.05) {
  rows <- rows[is.na(rows$error %||% NA), , drop = FALSE]
  key <- interaction(rows$estimator, rows$model, rows$max.month, drop = TRUE)
  out <- lapply(split(rows, key), function(d) {
    pars <- if (d$model[1] == "MKT") "MKT" else c("MKT", d$model[1])
    do.call(rbind, lapply(pars, function(p) {
      est <- if (p == "MKT") d$MKT else d$Coef
      se <- if (p == "MKT") d$SE.MKT else d$SE.Coef
      truth <- if (p %in% names(theta0)) theta0[[p]] else 0
      t <- (est - truth) / se
      pv <- ifelse(is.finite(t), fixedb.pvalue(t^2, d$hac.b, 1), NA_real_)
      tested <- is.finite(pv)
      data.frame(estimator = d$estimator[1], model = d$model[1], parameter = p, H = d$max.month[1] / 12, reps = nrow(d),
        truth = truth, mean = mean(est), bias = mean(est) - truth, mcse.bias = sd(est) / sqrt(nrow(d)), sd = sd(est),
        rmse = sqrt(mean((est - truth)^2)), median = median(est), on.bound = mean(d$at.bound),
        size = if (any(tested)) mean(pv[tested] < level) else NA_real_, mcse.size = if (any(tested)) sqrt(level * (1 - level) / sum(tested)) else NA_real_,
        tested = sum(tested), mean.se = mean(se, na.rm = TRUE), seconds = mean(d$seconds), stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, out)
  out[order(out$estimator, out$model, out$parameter, out$H), ]
}
`%||%` <- function(a, b) if (is.null(a)) b else a

######################
### Step 4 Monte Carlo: analysis
### Revision plan, Section I, "Step 4: scope and compute plan"
######################
#
# Criterion targets (mc.targets): theta_H^dagger of an estimator, model and
# horizon is the minimizer of the expected criterion, approximated by the
# criterion averaged over replications of the design (the same draws as the
# estimates, seeds 1, ..., R):
#   L2:      Q(theta) = mean over replications of mean_i e_i(theta)^2
#            (the loss whose minimizer the estimator targets; the reported
#            objective sqrt(sum e^2) / N is a monotone transform per draw, not
#            under the expectation)
#   U_{v,H}: Q(theta) = mean over replications of the U criterion.
# Like the estimator, the target is the minimizer over the admissible domain
# Theta, built from the replication whose grids start earliest (the most
# months, so the most restrictive floor): one parameter by Brent on the
# feasible interval; two by fit.admissible() started at theta_0 (wall outside
# Theta, polish along the boundary). A target on a bound is flagged: the
# averaged criterion then falls towards the edge (e.g. the exploding alpha at
# H = 0, T5). An unconstrained search gave a nonsense target there (the first
# version; corrected before any use). Without a shared vintage shock the U target is
# theta_0 exactly (E U_v = A^2 (mu_v^2 + c_v), zero at theta_0; T7); it is
# computed only as a check.
#
# mc.summary(): per estimator, model, parameter and horizon: mean, bias
# against the true value and its Monte Carlo SE, sd, RMSE, share on a bound,
# rejection rate of the two-sided fixed-b t-test (5%) of the true value and,
# with targets, the decomposition bias = target shift (target - truth) +
# finite-sample bias (mean - target) and the rejection rate at the target.
# Estimates on a bound have no SE: counted as not tested (column `tested`).
if (!exists("fixedb.pvalue")) source("helper/fixedb.R")
`%||%` <- function(a, b) if (is.null(a)) b else a

# Average criterion of an estimator over prepared replications (fast structures)
mc.mean.criterion <- function(est, fps, max.month) {
  function(p) {
    v <- vapply(fps, function(fp) {
      e <- fast.unit.errors(p, fp, max.month)
      if (any(!is.finite(e))) return(Inf)
      if (est$criterion == "L2") mean(e^2) else fast.u.criterion(e, fp, vintage.weights = est$vintage.weights)
    }, 0)
    mean(v)
  }
}

# Targets of every estimator, model and horizon of a design (spec from
# mc.load.designs()), from replications `reps`; parallel over cells.
mc.targets <- function(spec, public, reps = 1:200, workers = 1, check.u = FALSE) {
  levels <- unique(vapply(MC.ESTIMATORS[spec$estimators], `[[`, "", "level"))
  draws <- lapply(reps, function(s) mc.draw(spec$design, s, public))
  fps <- lapply(setNames(levels, levels), function(l) lapply(draws, function(d) mc.as.fast(d, public, l)))
  th0 <- spec$design$theta0
  first <- lapply(fps, function(l) l[[which.min(vapply(l, function(fp) min(fp$grids$k0), 0))]])
  cells <- expand.grid(estimator = spec$estimators, model = spec$models, max.month = spec$horizons, stringsAsFactors = FALSE)
  one <- function(i) {
    est <- MC.ESTIMATORS[[cells$estimator[i]]]; model <- cells$model[i]; m <- cells$max.month[i]
    par <- if (model == "MKT") c(MKT = unname(th0["MKT"])) else setNames(c(th0[["MKT"]], if (model %in% names(th0)) th0[[model]] else 0), c("MKT", model))
    exact <- est$criterion == "U" && spec$design$shared.sd == 0 && !check.u
    if (exact) return(data.frame(cells[i, ], MKT.target = par[["MKT"]], Coef.target = if (length(par) > 1) par[[2]] else NA_real_,
                                 method = "theta_0 (T7)", target.bound = ""))
    f <- mc.mean.criterion(est, fps[[est$level]], m)
    con <- mc.domain(first[[est$level]], names(par))
    if (length(par) == 1) {
      o <- optimize(function(b) f(c(MKT = b)), feasible.interval(con), tol = 1e-6)
      out <- setNames(o$minimum, "MKT")
    } else {
      fit <- fit.admissible(f, list(par), con, optim.parscale(par), smooth = TRUE, reltol = 1e-12)
      out <- fit$par
    }
    data.frame(cells[i, ], MKT.target = unname(out[1]), Coef.target = if (length(out) > 1) unname(out[2]) else NA_real_,
               method = "averaged criterion", target.bound = paste(active.constraints(out, con), collapse = ","))
  }
  res <- if (workers > 1) parallel::mclapply(seq_len(nrow(cells)), one, mc.cores = workers) else lapply(seq_len(nrow(cells)), one)
  bad <- !vapply(res, is.data.frame, TRUE)
  if (any(bad)) stop("target computation failed: ", as.character(res[[which(bad)[1]]])[1])
  out <- do.call(rbind, res); rownames(out) <- NULL
  out$reps.target <- length(reps)
  out
}

mc.summary <- function(rows, theta0, targets = NULL, level = 0.05) {
  rows <- rows[is.na(rows$error %||% NA), , drop = FALSE]
  key <- interaction(rows$estimator, rows$model, rows$max.month, drop = TRUE)
  p.of <- function(est, se, null, b) {
    t <- (est - null) / se
    ok <- is.finite(t) & is.finite(b)
    p <- rep(NA_real_, length(t))
    if (any(ok)) p[ok] <- fixedb.pvalue(t[ok]^2, b[ok], 1)
    p
  }
  out <- lapply(split(rows, key), function(d) {
    pars <- if (d$model[1] == "MKT") "MKT" else c("MKT", d$model[1])
    tg <- if (is.null(targets)) NULL else targets[targets$estimator == d$estimator[1] & targets$model == d$model[1] & targets$max.month == d$max.month[1], ]
    do.call(rbind, lapply(pars, function(p) {
      est <- if (p == "MKT") d$MKT else d$Coef
      se <- if (p == "MKT") d$SE.MKT else d$SE.Coef
      truth <- if (p %in% names(theta0)) theta0[[p]] else 0
      pv <- p.of(est, se, truth, d$hac.b)
      tested <- is.finite(pv)
      target <- if (is.null(tg) || !nrow(tg)) NA_real_ else if (p == "MKT") tg$MKT.target[1] else tg$Coef.target[1]
      pv.t <- if (is.finite(target)) p.of(est, se, target, d$hac.b) else rep(NA_real_, nrow(d))
      data.frame(estimator = d$estimator[1], model = d$model[1], parameter = p, H = d$max.month[1] / 12, reps = nrow(d),
        truth = truth, target = target, mean = mean(est), bias = mean(est) - truth, mcse.bias = sd(est) / sqrt(nrow(d)),
        shift = target - truth, fs.bias = mean(est) - target, sd = sd(est), rmse = sqrt(mean((est - truth)^2)),
        on.bound = mean(d$at.bound), size = if (any(tested)) mean(pv[tested] < level) else NA_real_,
        size.target = if (any(is.finite(pv.t))) mean(pv.t[is.finite(pv.t)] < level) else NA_real_,
        mcse.size = if (any(tested)) sqrt(level * (1 - level) / sum(tested)) else NA_real_, tested = sum(tested),
        mean.se = mean(se, na.rm = TRUE), se.ratio = mean(se, na.rm = TRUE) / sd(est), neff = mean(d$neff.MKT, na.rm = TRUE),
        stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, out)
  out[order(out$estimator, out$model, out$parameter, out$H), ]
}

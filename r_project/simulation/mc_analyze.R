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
mc.mean.criterion <- function(est, fps, max.month, fixed = NULL) {
  full.par <- if (length(fixed)) function(p) c(fixed, p) else identity
  function(p) {
    v <- vapply(fps, function(fp) {
      e <- fast.unit.errors(full.par(p), fp, max.month)
      if (any(!is.finite(e))) return(Inf)
      if (est$criterion == "L2") mean(e^2)
      else if (est$criterion == "U.cross") fast.u.cross.criterion(e, fp, est$band)
      else fast.u.criterion(e, fp, mc.u.sizes(est, fp), vintage.weights = est$vintage.weights)
    }, 0)
    mean(v)
  }
}

# Smallest root of the IV moment averaged over replications, on the interval
# where every replication's moment is defined (grid of 401 points, uniroot)
mc.mean.iv.root <- function(est, fps, max.month) {
  zs <- lapply(fps, function(fp) fast.horizon.scale(est$reference, fp, max.month))
  g <- function(b) mean(vapply(seq_along(fps), function(r) { e <- fast.unit.errors(c(MKT = b), fps[[r]], 0); if (any(!is.finite(e))) NA_real_ else mean(zs[[r]] * e) }, 0))
  iv <- vapply(fps, function(fp) feasible.interval(mc.domain(fp, "MKT")), numeric(2)); br <- c(max(iv[1, ]), min(iv[2, ]))
  grid <- seq(br[1], br[2], length.out = 401); gv <- vapply(grid, g, 0)
  ch <- which(is.finite(gv[-1]) & is.finite(gv[-length(gv)]) & sign(gv[-1]) != sign(gv[-length(gv)]))
  if (!length(ch)) return(list(root = NA_real_, n.roots = 0L))
  list(root = uniroot(g, grid[c(ch[1], ch[1] + 1)], tol = 1e-10)$root, n.roots = length(ch))
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
    fixed <- if (model == "Alpha.only") c(MKT = unname(th0["MKT"])) else NULL
    par <- if (model == "MKT") c(MKT = unname(th0["MKT"]))
           else if (model == "Alpha.only") c(Alpha = if ("Alpha" %in% names(th0)) unname(th0["Alpha"]) else 0)
           else setNames(c(th0[["MKT"]], if (model %in% names(th0)) th0[[model]] else 0), c("MKT", model))
    # U: theta_0 by T7 only without vintage-level effects (no shared shock, no path shock, no calendar shock)
    exact <- est$criterion == "U" && spec$design$shared.sd == 0 && is.null(spec$design$path.shock) && is.null(spec$design$cal.shock) && !check.u
    if (exact) return(data.frame(cells[i, ],
                                 MKT.target = if (model == "Alpha.only") unname(th0["MKT"]) else par[["MKT"]],
                                 Coef.target = if (model == "Alpha.only") par[["Alpha"]] else if (length(par) > 1) par[[2]] else NA_real_,
                                 method = "theta_0 (T7)", target.bound = ""))
    # IV: the instrument is a function of the factor path. Without a path
    # shock every design prices each fund at theta_0 given the path (shared
    # shocks have mean one and are independent of it), so
    # E[z_v PV_v(theta_0) | path] = 0 (T2) and the target is theta_0. With a
    # path shock (X) fails: the target is the smallest root of the moment
    # averaged over the replications (mc.mean.iv.root), as in mc.fit.iv.
    if (est$criterion == "IV" && is.null(spec$design$path.shock)) return(data.frame(cells[i, ], MKT.target = par[["MKT"]], Coef.target = NA_real_,
                                                 method = "theta_0 (moment)", target.bound = ""))
    if (est$criterion == "IV") { r <- mc.mean.iv.root(est, fps[[est$level]], m)
      return(data.frame(cells[i, ], MKT.target = r$root, Coef.target = NA_real_, method = sprintf("averaged moment, smallest of %d root(s)", r$n.roots), target.bound = "")) }
    f <- mc.mean.criterion(est, fps[[est$level]], m, fixed = fixed)
    con <- mc.domain(first[[est$level]], names(par), fixed = fixed)
    if (length(par) == 1) {
      # with resampled paths each replication has its own domain: the target
      # is searched where every replication's criterion is finite
      br <- if (is.null(spec$design$path)) feasible.interval(con) else {
        iv <- vapply(fps[[est$level]], function(fp) feasible.interval(mc.domain(fp, names(par), fixed = fixed)), numeric(2))
        c(max(iv[1, ]), min(iv[2, ]))
      }
      if (est$criterion == "U.cross") {
        o <- mc.opt.1d.global(f, br, names(par))
        out <- setNames(o$par, names(par))
      } else {
        o <- optimize(function(x) f(setNames(x, names(par))), br, tol = 1e-6)
        out <- setNames(o$minimum, names(par))
      }
    } else {
      if (!is.null(spec$design$path)) stop("targets of two-parameter models with resampled paths are not implemented")
      ps <- optim.parscale(par)
      fit <- fit.admissible(f, admissible.starts(par, con, ps), con, ps, smooth = TRUE, reltol = 1e-12)
      out <- fit$par
    }
    data.frame(cells[i, ],
               MKT.target = if (model == "Alpha.only") unname(th0["MKT"]) else unname(out[1]),
               Coef.target = if (model == "Alpha.only") unname(out[1]) else if (length(out) > 1) unname(out[2]) else NA_real_,
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
    pars <- if (d$model[1] == "MKT") "MKT" else if (d$model[1] == "Alpha.only") "Alpha" else c("MKT", d$model[1])
    tg <- if (is.null(targets)) NULL else targets[targets$estimator == d$estimator[1] & targets$model == d$model[1] & targets$max.month == d$max.month[1], ]
    do.call(rbind, lapply(pars, function(p) {
      est <- if (p == "MKT") d$MKT else d$Coef
      se <- if (p == "MKT") d$SE.MKT else d$SE.Coef
      neff.vec <- if (p == "MKT") d$neff.MKT else if ("neff.Coef" %in% names(d)) d$neff.Coef else rep(NA_real_, nrow(d))
      truth <- if (p %in% names(theta0)) theta0[[p]] else 0
      pv <- p.of(est, se, truth, d$hac.b)
      tested <- is.finite(pv)
      target <- if (is.null(tg) || !nrow(tg)) NA_real_ else if (p == "MKT") tg$MKT.target[1] else tg$Coef.target[1]
      pv.t <- if (is.finite(target)) p.of(est, se, target, d$hac.b) else rep(NA_real_, nrow(d))
      sd.tested <- if (sum(tested) > 1L) sd(est[tested]) else NA_real_
      mean.se.tested <- if (any(tested)) mean(se[tested]) else NA_real_
      data.frame(estimator = d$estimator[1], model = d$model[1], parameter = p, H = d$max.month[1] / 12, reps = nrow(d),
        truth = truth, target = target, mean = mean(est), bias = mean(est) - truth, mcse.bias = sd(est) / sqrt(nrow(d)),
        shift = target - truth, fs.bias = mean(est) - target, sd = sd(est), sd.tested = sd.tested,
        rmse = sqrt(mean((est - truth)^2)),
        on.bound = mean(d$at.bound), size = if (any(tested)) mean(pv[tested] < level) else NA_real_,
        size.target = if (any(is.finite(pv.t))) mean(pv.t[is.finite(pv.t)] < level) else NA_real_,
        mcse.size = if (any(tested)) sqrt(level * (1 - level) / sum(tested)) else NA_real_, tested = sum(tested),
        mean.se = mean(se, na.rm = TRUE), se.ratio = mean(se, na.rm = TRUE) / sd(est),
        se.ratio.tested = if (is.finite(sd.tested) && sd.tested > 0) mean.se.tested / sd.tested else NA_real_,
        neff = if (any(is.finite(neff.vec))) mean(neff.vec, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, out)
  out[order(out$estimator, out$model, out$parameter, out$H), ]
}

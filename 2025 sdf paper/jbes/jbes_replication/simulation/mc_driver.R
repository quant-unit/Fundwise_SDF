######################
### Step 4 Monte Carlo: driver
### Revision plan, Section I, "Step 4: scope and compute plan"; tests in
### helper/test_mc.R (T3: estimates equal the pipeline's; T4: restart)
######################
#
# One replication: draw (simulation/mc_generate.R), then every estimator,
# model and horizon on the same draw (common random numbers). The estimator
# is the pipeline's (estim_model_optimized.R, r5):
#   - criteria: L2 sqrt(sum e^2) / N over vintage portfolios or funds;
#     U_{v,H} over the funds of each vintage, q_v = 1 or q_v = n_v;
#   - domain Theta (helper/admissible.R): beta_MKT in [-1, 4], alpha in
#     [-0.02, 0.02], monthly gross return >= 0.08 in every month of the window;
#   - fit: Brent on the feasible interval (one parameter); fit.admissible()
#     from admissible.starts() around (MKT = 1, others 0) (two parameters);
#   - SE: Bartlett sandwich, S = ceil(1.3 sqrt V), unit losses e^2 (L2) or
#     w_v U_v (U_{v,H}), numerical derivatives with delta = 1e-4.
# Only the criterion evaluation is the fast one (helper/objective.R, T1).
# No jackknife, no CV (scope principle (e)). Each row keeps what size and
# coverage need at any null: estimates, SEs, their covariance and b; fixed-b
# p-values are computed in the analysis (helper/fixedb.R).
#
# Checkpoints: <out.dir>/<design>/batch_<k>.rds, one per batch of replications
# (seeds = replication indices); a restart skips existing batches.

for (f in c("helper/result_cache.R", "helper/parallel_workers.R", "helper/getNPVs.R", "helper/objective.R", "helper/admissible.R",
           "helper/inference.R", "simulation/mc_generate.R")) invisible(capture.output(suppressMessages(source(f))))

MC.ESTIMATORS <- list(
  L2.vyp  = list(level = "vyp",  criterion = "L2"),
  L2.fund = list(level = "fund", criterion = "L2"),
  U.eq    = list(level = "fund", criterion = "U", vintage.weights = "equal"),
  U.n     = list(level = "fund", criterion = "U", vintage.weights = "n.funds"),
  # benchmark-path IV comparator (plan W9, Gate G3; one instrument, single
  # factor only): the root of V^-1 sum_v z_v PV_v(beta) = 0 with
  # z_v = A_{v,H}(reference) / A_{v,0}(reference) (fast.horizon.scale), a
  # function of the factor path only; z_v = 1 at H = 0 (the unweighted mean
  # pricing moment). Fitted by mc.fit.iv().
  IV      = list(level = "vyp",  criterion = "IV", reference = c(MKT = 1)),
  # fund-size weighting (designs with sizes; log hardening, design B): L2 on
  # FW vintage portfolios and U with sizes in the pair normalization, as the
  # pipeline's FW files (use.size)
  L2.vyp.fw = list(level = "vyp.fw", criterion = "L2"),
  U.eq.fw   = list(level = "fund.fw", criterion = "U", vintage.weights = "equal", size = TRUE),
  # leave-band-out cross-vintage pair-product diagnostics (Proposition 2(iv))
  U.cross.m0 = list(level = "vyp", criterion = "U.cross", band = 0L),
  U.cross.m5 = list(level = "vyp", criterion = "U.cross", band = 5L),
  U.cross.m7 = list(level = "vyp", criterion = "U.cross", band = 7L)
)
# sizes passed to the U criterion: the units' sizes for FW U, else 1
mc.u.sizes <- function(est, fp) if (isTRUE(est$size)) fp$size else rep(1, fp$n.units)
MC.SETTINGS <- list(gross.floor = 0.08, mkt.bounds = c(-1, 4), alpha.bounds = c(-0.02, 0.02),
                    bandwidth = "llsw", delta = 1e-4, wall = 1e10)

# Criterion of one estimator on a fast structure (same values as err.sqr.calc)
mc.criterion <- function(est, fp, max.month, fixed = NULL) {
  full.par <- if (length(fixed)) function(p) c(fixed, p) else identity
  if (est$criterion == "L2") {
    function(p) {
      e <- fast.unit.errors(full.par(p), fp, max.month)
      if (any(!is.finite(e))) return(MC.SETTINGS$wall)
      sqrt(sum(e^2)) / length(e)
    }
  } else if (est$criterion == "U.cross") {
    function(p) {
      e <- fast.unit.errors(full.par(p), fp, max.month)
      if (any(!is.finite(e))) return(MC.SETTINGS$wall)
      u <- fast.u.cross.criterion(e, fp, est$band)
      if (!is.finite(u)) MC.SETTINGS$wall else u
    }
  } else {
    function(p) {
      e <- fast.unit.errors(full.par(p), fp, max.month)
      if (any(!is.finite(e))) return(MC.SETTINGS$wall)
      u <- fast.u.criterion(e, fp, mc.u.sizes(est, fp), vintage.weights = est$vintage.weights)
      if (!is.finite(u)) MC.SETTINGS$wall else u
    }
  }
}

# Sandwich SEs with the pipeline's units and unit losses (estim.CoVM)
mc.sandwich <- function(est, fp, par, max.month, fixed = NULL) {
  full.par <- if (length(fixed)) function(p) c(fixed, p) else identity
  if (est$criterion == "L2") {
    units <- lapply(seq_len(fp$n.units), function(j) fast.subset(fp, seq_len(fp$n.units) == j))
    loss <- function(p, u) fast.unit.errors(full.par(p), u, max.month)^2
    vint <- fp$vintage
  } else {
    w <- u.unit.weights(fp$group.n, est$vintage.weights)
    keep <- names(w)[!is.na(w)]
    units <- lapply(keep, function(k) fast.subset(fp, fp$key == k))
    names(units) <- keep
    loss <- function(p, u) w[[u$key[1]]] * fast.u.criterion(fast.unit.errors(full.par(p), u, max.month), u, mc.u.sizes(est, u), vintage.weights = est$vintage.weights)
    vint <- vapply(units, function(u) u$vintage[1], 0)
  }
  sandwich.covariance(units, par, loss, vint, bw = MC.SETTINGS$bandwidth, delta = MC.SETTINGS$delta)
}

# Domain months of a replication: every calendar month from the earliest
# vintage origin to T.max (the pipeline's admissible.months of the type)
mc.domain <- function(fp, par.names, fixed = NULL) {
  rows <- min(fp$grids$k0):nrow(fp$X)
  months <- data.frame(Date = rows, fp$X[rows, , drop = FALSE])
  if (length(fixed)) {
    for (nm in names(fixed)) months$RF <- months$RF + fixed[[nm]] * months[[nm]]
  }
  admissible.constraints(months, par.names, MC.SETTINGS$gross.floor, MC.SETTINGS$mkt.bounds, MC.SETTINGS$alpha.bounds)
}

# Global 1D minimizer on [br[1], br[2]] for potentially non-convex / indefinite
# criteria (U.cross): evaluates a 101-point grid including both endpoints and
# refines up to 5 lowest interior local brackets with optimize().
mc.opt.1d.global <- function(f, br, par.name) {
  grid <- seq(br[1], br[2], length.out = 101L)
  vals <- vapply(grid, function(x) f(setNames(x, par.name)), 0)
  best.x <- br[1]; best.v <- vals[1]
  if (vals[length(vals)] < best.v) { best.x <- br[2]; best.v <- vals[length(vals)] }
  n <- length(vals)
  loc <- which(vals[2:(n - 1L)] <= vals[1:(n - 2L)] & vals[2:(n - 1L)] <= vals[3:n]) + 1L
  idx <- unique(c(which.min(vals), loc))
  idx <- idx[idx > 1L & idx < n]
  if (length(idx) > 5L) idx <- idx[order(vals[idx])[1:5]]
  for (i in idx) {
    o <- stats::optimize(function(x) f(setNames(x, par.name)), c(grid[i - 1L], grid[i + 1L]), tol = 1e-10)
    if (is.finite(o$objective) && o$objective < best.v) { best.x <- o$minimum; best.v <- o$objective }
  }
  list(par = best.x, value = best.v, convergence = 0L)
}

# IV comparator, single factor: the root of gbar(beta) = V^-1 sum_v z_v PV_v(beta)
# on the feasible interval of Theta. gbar is evaluated on a grid of 201 points;
# each sign change brackets a root (n.roots); the smallest root is taken (the
# economic root lies below the growth-optimal one, T8) and refined by uniroot.
# Without a sign change: the minimizer of gbar^2 (n.roots = 0), reported with
# no SE (the moment has no zero in Theta, so there is no estimate to test;
# like an estimate on a bound, it counts as not tested). SE: sandwich
# of the just-identified moment, Var = Omega / (V G^2), Omega the Bartlett sum
# of the vintage scores s_v = z_v PV_v over vintage-year distances with the
# estimator's S = ceil(1.3 sqrt V), G = gbar'(beta) by central differences.
mc.fit.iv <- function(est, fp, model, max.month) {
  if (model != "MKT") stop("IV comparator: single factor only")
  t0 <- proc.time()[["elapsed"]]
  z <- fast.horizon.scale(est$reference, fp, max.month)
  pv <- function(b) fast.unit.errors(c(MKT = b), fp, 0)
  gbar <- function(b) { e <- pv(b); if (any(!is.finite(e))) NA_real_ else mean(z * e) }
  con <- mc.domain(fp, "MKT"); br <- feasible.interval(con)
  grid <- seq(br[1], br[2], length.out = 201); gv <- vapply(grid, gbar, 0)
  ch <- which(is.finite(gv[-1]) & is.finite(gv[-length(gv)]) & sign(gv[-1]) != sign(gv[-length(gv)]))
  if (length(ch)) {
    i <- ch[1]
    b.hat <- if (gv[i] == 0) grid[i] else if (gv[i + 1] == 0) grid[i + 1] else uniroot(gbar, grid[c(i, i + 1)], tol = 1e-12)$root
  } else {
    b.hat <- optimize(function(b) { g <- gbar(b); if (is.finite(g)) g^2 else MC.SETTINGS$wall }, br, tol = 1e-10)$minimum
  }
  par <- c(MKT = b.hat)
  act <- active.constraints(par, con)
  se <- NA_real_; S <- NA_real_; V <- fp$n.units; b <- NA_real_
  if (!length(act) && length(ch)) {
    d <- MC.SETTINGS$delta
    G <- (gbar(b.hat + d) - gbar(b.hat - d)) / (2 * d)
    s <- z * pv(b.hat)
    S <- if (identical(MC.SETTINGS$bandwidth, "llsw")) hac.bandwidth.llsw(V) else MC.SETTINGS$bandwidth
    K <- outer(fp$vintage, fp$vintage, bartlett.weight, bw = S)
    omega <- sum(K * outer(s, s)) / V
    se <- sqrt(omega / (V * G^2)); b <- S / V
  }
  data.frame(model = model, max.month = max.month, MKT = b.hat, Coef = NA_real_, SE.MKT = se, SE.Coef = NA_real_,
    cov.MKT.Coef = NA_real_, value = gbar(b.hat), convergence = 0L, fit.method = if (length(ch)) "root" else "min.square",
    n.endpoints = 1L, at.bound = length(act) > 0, active = paste(act, collapse = ","), min.G = min.gross.return(par, con)$value,
    hac.S = S, hac.V = V, hac.b = b, Wald.stat = NA_real_, q = 1L, neff.MKT = NA_real_, neff.Coef = NA_real_, se.singular = NA,
    n.roots = length(ch), seconds = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE)
}

# One fit: estimator est, model ("MKT", "Alpha.only", or a second parameter), horizon (months)
mc.fit <- function(est, fp, model, max.month) {
  if (est$criterion == "IV") return(mc.fit.iv(est, fp, model, max.month))
  t0 <- proc.time()[["elapsed"]]
  fixed <- if (model == "Alpha.only") c(MKT = 1) else NULL
  par <- if (model == "MKT") c(MKT = 1) else if (model == "Alpha.only") c(Alpha = 0) else setNames(c(1, 0), c("MKT", model))
  con <- mc.domain(fp, names(par), fixed = fixed)
  f <- mc.criterion(est, fp, max.month, fixed = fixed)
  if (length(par) == 1) {
    br <- feasible.interval(con)
    if (est$criterion == "U.cross") {
      ans <- mc.opt.1d.global(f, br, names(par))
      par[] <- ans$par; value <- ans$value; conv <- ans$convergence; method <- "grid+opt"; n.end <- 1L
    } else {
      ans <- stats::optim(par, function(p) f(setNames(p, names(par))), method = "Brent", lower = br[1], upper = br[2])
      par[] <- ans$par; value <- ans$value; conv <- ans$convergence; method <- "Brent"; n.end <- 1L
    }
  } else {
    ps <- optim.parscale(par)
    fit <- fit.admissible(f, admissible.starts(par, con, ps), con, ps, smooth = TRUE)
    par <- fit$par; value <- fit$value; conv <- fit$convergence; method <- fit$method; n.end <- length(fit$endpoints)
  }
  act <- active.constraints(par, con)
  cov <- if (length(act) == 0 && est$criterion != "U.cross") mc.sandwich(est, fp, par, max.month, fixed = fixed) else NULL # SEs not valid on a bound or indefinite U.cross
  se <- function(v) if (is.null(cov)) NA_real_ else unname(cov$SE[paste0("SE.", v)])
  if (model == "Alpha.only") {
    return(data.frame(model = model, max.month = max.month, MKT = 1, Coef = unname(par["Alpha"]),
      SE.MKT = NA_real_, SE.Coef = se("Alpha"), cov.MKT.Coef = NA_real_,
      value = value, convergence = conv, fit.method = method, n.endpoints = n.end,
      at.bound = length(act) > 0, active = paste(act, collapse = ","), min.G = min.gross.return(par, con)$value,
      hac.S = if (is.null(cov)) NA_real_ else cov$hac.S, hac.V = if (is.null(cov)) NA_real_ else cov$hac.V,
      hac.b = if (is.null(cov)) NA_real_ else cov$hac.b, Wald.stat = if (is.null(cov)) NA_real_ else cov$Wald.stat,
      q = 1L, neff.MKT = NA_real_, neff.Coef = if (is.null(cov)) NA_real_ else cov$neff[["Alpha"]],
      se.singular = if (is.null(cov)) NA else isTRUE(cov$se.singular),
      seconds = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE))
  }
  second <- if (length(par) > 1) names(par)[2] else NA_character_
  data.frame(model = model, max.month = max.month, MKT = unname(par["MKT"]),
    Coef = if (is.na(second)) NA_real_ else unname(par[second]),
    SE.MKT = se("MKT"), SE.Coef = if (is.na(second)) NA_real_ else se(second),
    cov.MKT.Coef = if (is.null(cov) || is.na(second)) NA_real_ else cov$CoVM["MKT", second],
    value = value, convergence = conv, fit.method = method, n.endpoints = n.end,
    at.bound = length(act) > 0, active = paste(act, collapse = ","), min.G = min.gross.return(par, con)$value,
    hac.S = if (is.null(cov)) NA_real_ else cov$hac.S, hac.V = if (is.null(cov)) NA_real_ else cov$hac.V,
    hac.b = if (is.null(cov)) NA_real_ else cov$hac.b, Wald.stat = if (is.null(cov)) NA_real_ else cov$Wald.stat,
    q = length(par), neff.MKT = if (is.null(cov)) NA_real_ else cov$neff[["MKT"]],
    neff.Coef = if (is.null(cov) || is.na(second)) NA_real_ else cov$neff[[second]],
    se.singular = if (is.null(cov)) NA else isTRUE(cov$se.singular),
    seconds = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE)
}

# Every estimator, model and horizon of one replication
# (default: the pipeline's four estimators; the IV comparator is not in the
# empirical pipeline and is requested by name)
mc.replication <- function(design, seed, public, estimators = c("L2.vyp", "L2.fund", "U.eq", "U.n"), models = c("MKT", "Alpha"),
                           horizons = c(0, 60, 120, 180, 240, 360)) {
  dr <- mc.draw(design, seed, public)
  fps <- list(fund = NULL, vyp = NULL)
  for (lev in unique(vapply(MC.ESTIMATORS[estimators], `[[`, "", "level"))) fps[[lev]] <- mc.as.fast(dr, public, lev)
  out <- list()
  for (e in estimators) for (m in models) for (h in horizons) {
    est <- MC.ESTIMATORS[[e]]
    r <- mc.fit(est, fps[[est$level]], m, h)
    out[[length(out) + 1]] <- data.frame(rep = seed, estimator = e, r, stringsAsFactors = FALSE)
  }
  do.call(rbind.all.columns.mc, out) # IV rows add n.roots
}

# Run replications of one design in batches with checkpoints. spec: list with
# design (mc.design()), estimators, models, horizons. Returns all rows.
mc.run <- function(name, spec, reps, public, out.dir = "simulation/data_out_mc", batch.size = 50,
                   workers = if (exists("n.parallel.workers")) n.parallel.workers(length(reps)) else 1, verbose = TRUE) {
  reps <- as.integer(reps)
  if (!length(reps) || anyNA(reps) || any(reps < 1L) || anyDuplicated(reps)) {
    stop("Design ", name, ": reps must be unique positive integers")
  }
  dir <- file.path(out.dir, name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  batches <- split(reps, ceiling(reps / batch.size))
  prov <- list(git.commit = tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1], error = function(e) NA),
               git.dirty = tryCatch(length(system("git status --porcelain --untracked-files=no", intern = TRUE, ignore.stderr = TRUE)) > 0, error = function(e) NA),
               spec.md5 = mc.spec.md5(spec), public.md5 = value.md5(public), settings = MC.SETTINGS)
  failed.batches <- character()
  for (k in names(batches)) {
    file <- file.path(dir, sprintf("batch_%04d.rds", as.integer(k)))
    if (file.exists(file)) {
      old <- readRDS(file)
      if (!identical(old$provenance$spec.md5, prov$spec.md5)) stop("Checkpoint ", file, " was written for another design specification")
      if (!identical(old$provenance$public.md5, prov$public.md5)) stop("Checkpoint ", file, " was written for different public factors")
      if (!identical(old$provenance$settings, prov$settings)) stop("Checkpoint ", file, " was written with different MC.SETTINGS")
      if (!identical(as.integer(old$reps), as.integer(batches[[k]]))) {
        stop("Checkpoint ", file, " has different replication IDs than requested")
      }
      if ("error" %in% names(old$rows) && any(!is.na(old$rows$error))) {
        stop("Checkpoint ", file, " contains failed replication rows")
      }
      if (!identical(sort(unique(as.integer(old$rows$rep))), sort(as.integer(batches[[k]])))) {
        stop("Checkpoint ", file, " rows do not match its recorded replication IDs")
      }
      if (verbose) cat("skip", basename(file), "\n")
      next
    }
    t0 <- Sys.time()
    run1 <- function(s) tryCatch(mc.replication(spec$design, s, public, spec$estimators, spec$models, spec$horizons),
                                 error = function(e) data.frame(rep = s, error = conditionMessage(e)))
    # dynamic scheduling: one replication per job, so no worker idles while
    # another works through a prescheduled share (slow designs); rows do not
    # depend on the schedule (each replication is seeded by its index)
    rows <- if (workers > 1) parallel::mclapply(batches[[k]], run1, mc.cores = workers, mc.preschedule = FALSE) else lapply(batches[[k]], run1)
    # a worker that died (e.g. out of memory) returns a try-error, not rows
    rows <- lapply(seq_along(rows), function(j) {
      x <- rows[[j]]
      if (is.data.frame(x)) x else data.frame(rep = batches[[k]][j], error = paste("worker failed:", as.character(x)[1]))
    })
    rows <- do.call(rbind.all.columns.mc, rows)
    failed <- if ("error" %in% names(rows)) unique(rows$rep[!is.na(rows$error)]) else integer()
    out <- list(provenance = prov, name = name, reps = batches[[k]], rows = rows)
    if (length(failed)) {
      # not a finished batch: kept for inspection, recomputed by a restart
      saveRDS(out, sub("\\.rds$", "_failed.rds", file))
      failed.batches <- c(failed.batches, k)
      warning(name, " batch ", k, ": ", length(failed), " replication(s) failed, e.g. ", rows$error[!is.na(rows$error)][1], "; batch not marked finished")
    } else {
      saveRDS(out, file)
      unlink(sub("\\.rds$", "_failed.rds", file)) # a failed earlier attempt is superseded
    }
    if (verbose) cat(sprintf("%s batch %s: %d replications in %.0f s%s\n", name, k, length(batches[[k]]), as.numeric(difftime(Sys.time(), t0, units = "secs")),
                             if (length(failed)) sprintf(", %d FAILED", length(failed)) else ""))
  }
  if (length(failed.batches)) {
    stop("Design ", name, ": failed replications in batch(es) ", paste(failed.batches, collapse = ", "))
  }
  mc.read(name, out.dir, expected.spec = spec, expected.public = public,
          expected.settings = MC.SETTINGS, expected.reps = reps, allow.extra.reps = TRUE)
}

# Hash of what determines the rows: the design, estimators, models and
# horizons, not the number of replications (extending a design keeps its
# finished batches)
mc.spec.md5 <- function(spec) value.md5(spec[c("design", "estimators", "models", "horizons")])

rbind.all.columns.mc <- function(...) {
  l <- list(...)
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[, cols] }))
}

# All rows of a design, with each file's provenance checked for consistency
mc.read <- function(name, out.dir = "simulation/data_out_mc",
                    expected.spec = NULL, expected.public = NULL, expected.settings = NULL,
                    expected.reps = NULL, allow.extra.reps = FALSE) {
  files <- sort(list.files(file.path(out.dir, name), "^batch_\\d+\\.rds$", full.names = TRUE))
  if (!length(files)) {
    if (!is.null(expected.reps) && length(expected.reps) > 0L) {
      stop("Design ", name, ": no finished batches found (expected ", length(expected.reps), " replications)")
    }
    return(data.frame())
  }
  b <- lapply(files, readRDS)
  md5 <- unique(vapply(b, function(x) x$provenance$spec.md5, ""))
  if (length(md5) > 1L) stop("Design ", name, ": batches from different specifications")
  pub.md5 <- unique(vapply(b, function(x) x$provenance$public.md5, ""))
  if (length(pub.md5) > 1L) stop("Design ", name, ": batches from different public factors")
  st0 <- b[[1]]$provenance$settings
  if (!all(vapply(b, function(x) identical(x$provenance$settings, st0), TRUE))) {
    stop("Design ", name, ": batches from different MC.SETTINGS")
  }
  if (!is.null(expected.spec) && !identical(md5, mc.spec.md5(expected.spec))) {
    stop("Design ", name, ": batch specification hash does not match current spec")
  }
  if (!is.null(expected.public)) {
    exp.pub.md5 <- if (is.character(expected.public) && length(expected.public) == 1L) expected.public else value.md5(expected.public)
    if (!identical(pub.md5, exp.pub.md5)) {
      stop("Design ", name, ": batch public-factor hash does not match current public factors")
    }
  }
  if (!is.null(expected.settings) && !identical(st0, expected.settings)) {
    stop("Design ", name, ": batch MC.SETTINGS do not match current settings")
  }
  batch.reps <- lapply(b, function(x) as.integer(x$reps))
  all.reps <- unlist(batch.reps, use.names = FALSE)
  if (anyNA(all.reps) || anyDuplicated(all.reps)) {
    stop("Design ", name, ": duplicate or missing replication IDs across batches")
  }
  rows <- do.call(rbind.all.columns.mc, lapply(b, `[[`, "rows"))
  if ("error" %in% names(rows) && any(!is.na(rows$error))) {
    stop("Design ", name, ": finished batch contains failed replication rows")
  }
  row.reps <- sort(unique(as.integer(rows$rep)))
  if (!identical(row.reps, sort(all.reps))) {
    stop("Design ", name, ": row replication IDs do not match batch metadata")
  }
  if (!is.null(expected.reps)) {
    exp.reps <- sort(as.integer(expected.reps))
    if (allow.extra.reps) {
      missing.reps <- setdiff(exp.reps, row.reps)
      if (length(missing.reps)) {
        stop("Design ", name, ": missing ", length(missing.reps), " requested replication(s), e.g. ", missing.reps[1])
      }
    } else if (!identical(row.reps, exp.reps)) {
      missing.reps <- setdiff(exp.reps, row.reps)
      extra.reps <- setdiff(row.reps, exp.reps)
      stop(sprintf("Design %s: replication mismatch (missing %d, extra %d)", name, length(missing.reps), length(extra.reps)))
    }
  }
  attr(rows, "provenance") <- lapply(b, `[[`, "provenance")
  rows
}

# Designs from config/mc_designs.yaml: one spec per design (design, estimators,
# models, horizons, reps), defaults merged in. Empirical fund counts are read
# once per shape from the checked-in sample funnel table (aggregates only).
mc.load.designs <- function(file = "config/mc_designs.yaml", shape.file = "empirical/data_prepared_2026_c/preqin_sample_funnel_North America.csv") {
  cfg <- yaml::read_yaml(file)
  shapes <- list()
  lapply(setNames(names(cfg$designs), names(cfg$designs)), function(id) {
    d <- modifyList(cfg$defaults, cfg$designs[[id]])
    if (!is.null(d$consecutive)) {
      k <- d$consecutive # keys funds, count, last ("n" would be read as a YAML boolean)
      n <- setNames(rep(as.integer(k$funds), k$count), as.character(seq(k$last - k$count + 1, k$last)))
    } else {
      if (is.null(shapes[[d$shape]])) shapes[[d$shape]] <<- mc.data.shape(d$shape, file = shape.file)
      n <- shapes[[d$shape]]
      if (!is.null(d$vintages)) { r <- as.integer(strsplit(d$vintages, ":")[[1]]); n <- n[as.integer(names(n)) >= r[1] & as.integer(names(n)) <= r[2]] }
      if (!is.null(d$shape_scale)) n <- setNames(pmax(2L, as.integer(round(n * d$shape_scale))), names(n))
    }
    list(design = mc.design(n, theta0 = unlist(d$theta0), sigma = d$sigma, shared.sd = d$shared_sd, no.deals = d$no_deals,
                            investment.period = d$investment_period, max.holding.period = d$max_holding_period,
                            t.max = as.Date(d$t_max), quarterly = isTRUE(d$quarterly),
                            shared.rho = if (is.null(d$shared_rho)) 0 else d$shared_rho, path = d$path,
                            path.shock = if (is.null(d$path_shock)) NULL else list(a = d$path_shock$a, window.years = d$path_shock$window_years, level = d$path_shock$level),
                            sizes = d$sizes,
                            cal.shock = if (is.null(d$cal_shock)) NULL else list(sigma_W = d$cal_shock$sigma_W, phi_m = if (is.null(d$cal_shock$phi_m)) 0 else d$cal_shock$phi_m)),
         estimators = unlist(d$estimators), models = unlist(d$models), horizons = as.numeric(unlist(d$horizons)), reps = d$reps)
  })
}

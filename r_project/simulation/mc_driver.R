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
  IV      = list(level = "vyp",  criterion = "IV", reference = c(MKT = 1))
)
MC.SETTINGS <- list(gross.floor = 0.08, mkt.bounds = c(-1, 4), alpha.bounds = c(-0.02, 0.02),
                    bandwidth = "llsw", delta = 1e-4, wall = 1e10)

# Criterion of one estimator on a fast structure (same values as err.sqr.calc)
mc.criterion <- function(est, fp, max.month) {
  if (est$criterion == "L2") {
    function(p) {
      e <- fast.unit.errors(p, fp, max.month)
      if (any(!is.finite(e))) return(MC.SETTINGS$wall)
      sqrt(sum(e^2)) / length(e)
    }
  } else {
    function(p) {
      e <- fast.unit.errors(p, fp, max.month)
      if (any(!is.finite(e))) return(MC.SETTINGS$wall)
      u <- fast.u.criterion(e, fp, vintage.weights = est$vintage.weights)
      if (!is.finite(u)) MC.SETTINGS$wall else u
    }
  }
}

# Sandwich SEs with the pipeline's units and unit losses (estim.CoVM)
mc.sandwich <- function(est, fp, par, max.month) {
  if (est$criterion == "L2") {
    units <- lapply(seq_len(fp$n.units), function(j) fast.subset(fp, seq_len(fp$n.units) == j))
    loss <- function(p, u) fast.unit.errors(p, u, max.month)^2
    vint <- fp$vintage
  } else {
    w <- u.unit.weights(fp$group.n, est$vintage.weights)
    keep <- names(w)[!is.na(w)]
    units <- lapply(keep, function(k) fast.subset(fp, fp$key == k))
    names(units) <- keep
    loss <- function(p, u) w[[u$key[1]]] * fast.u.criterion(fast.unit.errors(p, u, max.month), u, vintage.weights = est$vintage.weights)
    vint <- vapply(units, function(u) u$vintage[1], 0)
  }
  sandwich.covariance(units, par, loss, vint, bw = MC.SETTINGS$bandwidth, delta = MC.SETTINGS$delta)
}

# Domain months of a replication: every calendar month from the earliest
# vintage origin to T.max (the pipeline's admissible.months of the type)
mc.domain <- function(fp, par.names) {
  rows <- min(fp$grids$k0):nrow(fp$X)
  months <- data.frame(Date = rows, fp$X[rows, , drop = FALSE])
  admissible.constraints(months, par.names, MC.SETTINGS$gross.floor, MC.SETTINGS$mkt.bounds, MC.SETTINGS$alpha.bounds)
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
    hac.S = S, hac.V = V, hac.b = b, Wald.stat = NA_real_, q = 1L, neff.MKT = NA_real_, se.singular = NA,
    n.roots = length(ch), seconds = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE)
}

# One fit: estimator est, model ("MKT" or a second parameter), horizon (months)
mc.fit <- function(est, fp, model, max.month) {
  if (est$criterion == "IV") return(mc.fit.iv(est, fp, model, max.month))
  t0 <- proc.time()[["elapsed"]]
  par <- if (model == "MKT") c(MKT = 1) else setNames(c(1, 0), c("MKT", model))
  con <- mc.domain(fp, names(par))
  f <- mc.criterion(est, fp, max.month)
  if (length(par) == 1) {
    br <- feasible.interval(con)
    ans <- stats::optim(par, function(p) f(setNames(p, "MKT")), method = "Brent", lower = br[1], upper = br[2])
    par[] <- ans$par; value <- ans$value; conv <- ans$convergence; method <- "Brent"; n.end <- 1L
  } else {
    ps <- optim.parscale(par)
    fit <- fit.admissible(f, admissible.starts(par, con, ps), con, ps, smooth = TRUE)
    par <- fit$par; value <- fit$value; conv <- fit$convergence; method <- fit$method; n.end <- length(fit$endpoints)
  }
  act <- active.constraints(par, con)
  cov <- if (length(act) == 0) mc.sandwich(est, fp, par, max.month) else NULL # SEs not valid on a bound
  se <- function(v) if (is.null(cov)) NA_real_ else unname(cov$SE[paste0("SE.", v)])
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
  dir <- file.path(out.dir, name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  batches <- split(reps, ceiling(reps / batch.size))
  prov <- list(git.commit = tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1], error = function(e) NA),
               git.dirty = tryCatch(length(system("git status --porcelain --untracked-files=no", intern = TRUE, ignore.stderr = TRUE)) > 0, error = function(e) NA),
               spec.md5 = mc.spec.md5(spec), public.md5 = value.md5(public), settings = MC.SETTINGS)
  for (k in names(batches)) {
    file <- file.path(dir, sprintf("batch_%04d.rds", as.integer(k)))
    if (file.exists(file)) {
      old <- readRDS(file)
      if (!identical(old$provenance$spec.md5, prov$spec.md5)) stop("Checkpoint ", file, " was written for another design specification")
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
      warning(name, " batch ", k, ": ", length(failed), " replication(s) failed, e.g. ", rows$error[!is.na(rows$error)][1], "; batch not marked finished")
    } else {
      saveRDS(out, file)
      unlink(sub("\\.rds$", "_failed.rds", file)) # a failed earlier attempt is superseded
    }
    if (verbose) cat(sprintf("%s batch %s: %d replications in %.0f s%s\n", name, k, length(batches[[k]]), as.numeric(difftime(Sys.time(), t0, units = "secs")),
                             if (length(failed)) sprintf(", %d FAILED", length(failed)) else ""))
  }
  mc.read(name, out.dir)
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
mc.read <- function(name, out.dir = "simulation/data_out_mc") {
  files <- sort(list.files(file.path(out.dir, name), "^batch_\\d+\\.rds$", full.names = TRUE))
  if (!length(files)) return(data.frame())
  b <- lapply(files, readRDS)
  md5 <- unique(vapply(b, function(x) x$provenance$spec.md5, ""))
  if (length(md5) > 1) stop("Design ", name, ": batches from different specifications")
  rows <- do.call(rbind.all.columns.mc, lapply(b, `[[`, "rows"))
  attr(rows, "provenance") <- lapply(b, `[[`, "provenance")
  rows
}

# Designs from config/mc_designs.yaml: one spec per design (design, estimators,
# models, horizons, reps), defaults merged in. Empirical fund counts are read
# once per shape from the local fund file (aggregates only).
mc.load.designs <- function(file = "config/mc_designs.yaml", shape.file = "empirical/data_prepared_2026_c/preqin_cashflows_EW_North America_NAV.csv") {
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
                            shared.rho = if (is.null(d$shared_rho)) 0 else d$shared_rho, path = d$path),
         estimators = unlist(d$estimators), models = unlist(d$models), horizons = as.numeric(unlist(d$horizons)), reps = d$reps)
  })
}

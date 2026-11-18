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
  U.n     = list(level = "fund", criterion = "U", vintage.weights = "n.funds")
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

# One fit: estimator est, model ("MKT" or a second parameter), horizon (months)
mc.fit <- function(est, fp, model, max.month) {
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
mc.replication <- function(design, seed, public, estimators = names(MC.ESTIMATORS), models = c("MKT", "Alpha"),
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
  do.call(rbind, out)
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
               spec.md5 = value.md5(spec), public.md5 = value.md5(public), settings = MC.SETTINGS)
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
    rows <- if (workers > 1) parallel::mclapply(batches[[k]], run1, mc.cores = workers) else lapply(batches[[k]], run1)
    rows <- do.call(rbind.all.columns.mc, rows)
    saveRDS(list(provenance = prov, name = name, reps = batches[[k]], rows = rows), file)
    if (verbose) cat(sprintf("%s batch %s: %d replications in %.0f s\n", name, k, length(batches[[k]]), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  mc.read(name, out.dir)
}

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

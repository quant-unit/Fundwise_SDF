######################
### Admissible parameter domain (Assumption A7) and constrained fitting
### Decision of 17 November 2026 (revision plan, State of evidence, section I, D2)
######################
#
# Theta: beta_MKT in brent.interval ([-1, 4], the single-factor bracket),
# alpha in [alpha.lower, alpha.upper], other loadings unrestricted, and the
# monthly gross benchmark return
#   G_t(theta) = 1 + r_f,t + theta' X_t >= GROSS.FLOOR
# in every month t of the type's full-sample window (the same months for the
# full-sample, CV and jackknife fits). GROSS.FLOOR = 0.08 is the lowest monthly
# gross return the single-factor domain already admits (beta = 4 in October
# 1987 gives 0.0786); it is not chosen from the multi-factor results.
# All constraints are linear in theta: ui %*% theta - ci >= 0.
#
# Fitting (fit.admissible):
#   1. Nelder-Mead from a fixed start set (pipeline start; points at 90% of the
#      distance to the boundary of Theta along the coordinate axes and, for
#      up to three parameters, the diagonals; directions scaled by parscale);
#      points outside Theta get a large value. For interior optima this is the
#      former estimator, started more widely.
#   2. If the best point lies within reach of a constraint, it is polished by
#      constrOptim (log barrier; BFGS for smooth criteria, Nelder-Mead for L1),
#      which can move along the boundary; the better of the two is kept.
#   An estimate is on a bound when a constraint's slack is below BOUND.TOL
#   (parameter units for the box, gross-return units for the floor).

GROSS.FLOOR <- 0.08
BOUND.TOL <- 1e-4
WALL.VALUE <- 1e10

# months: data frame with one row per month (columns: Date, RF and the factor
# columns; the Alpha column is 1). Returns the constraint set of Theta.
admissible.constraints <- function(months, par.names, floor = GROSS.FLOOR,
                                   mkt.bounds = c(-1, 4), alpha.bounds = c(-0.02, 0.02)) {
  months <- months[stats::complete.cases(months[, c("RF", par.names), drop = FALSE]), , drop = FALSE]
  ui <- as.matrix(months[, par.names, drop = FALSE])
  ci <- floor - 1 - months$RF
  kind <- rep("floor", nrow(ui))
  box <- function(name, lo, hi) {
    if (!(name %in% par.names)) return(NULL)
    e <- as.numeric(par.names == name)
    keep <- is.finite(c(lo, hi)) # an infinite bound is no constraint
    list(ui = rbind(e, -e)[keep, , drop = FALSE], ci = c(lo, -hi)[keep], kind = paste0(name, c(".lower", ".upper"))[keep])
  }
  for (b in list(box("MKT", mkt.bounds[1], mkt.bounds[2]), box("Alpha", alpha.bounds[1], alpha.bounds[2]))) {
    if (is.null(b)) next
    ui <- rbind(ui, b$ui); ci <- c(ci, b$ci); kind <- c(kind, b$kind)
  }
  dimnames(ui) <- list(NULL, par.names)
  list(ui = ui, ci = ci, kind = kind, dates = months$Date, rf = months$RF,
       X = as.matrix(months[, par.names, drop = FALSE]), floor = floor)
}

slack <- function(theta, con) as.numeric(con$ui %*% theta - con$ci)
in.domain <- function(theta, con) all(slack(theta, con) >= 0)

# lowest monthly gross benchmark return at theta, and its month
min.gross.return <- function(theta, con) {
  g <- 1 + con$rf + as.numeric(con$X %*% theta)
  i <- which.min(g)
  list(value = g[i], date = con$dates[i])
}

# constraints active at theta (slack below tol), by kind
active.constraints <- function(theta, con, tol = BOUND.TOL) unique(con$kind[slack(theta, con) < tol])

# feasible interval of a one-parameter model
feasible.interval <- function(con) {
  a <- con$ui[, 1]; c0 <- con$ci
  lo <- max(c(-Inf, (c0 / a)[a > 0])); hi <- min(c(Inf, (c0 / a)[a < 0]))
  c(lo, hi)
}

# Fixed start set: the centre (pipeline start) and points at frac of the distance
# to the boundary along +-axes and (k <= 3) +-diagonals, in parscale units.
admissible.starts <- function(center, con, parscale = rep(1, length(center)), frac = 0.9, max.step = 10) {
  k <- length(center)
  dirs <- rbind(diag(k), -diag(k))
  if (k >= 2 && k <= 3) {
    for (i in 1:(k - 1)) for (j in (i + 1):k) for (s in list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))) {
      d <- numeric(k); d[i] <- s[1]; d[j] <- s[2]; dirs <- rbind(dirs, d / sqrt(2))
    }
  }
  s0 <- slack(center, con)
  stopifnot(all(s0 > 0)) # the centre must be strictly inside Theta
  starts <- list(center)
  for (r in seq_len(nrow(dirs))) {
    d <- dirs[r, ] * parscale
    rate <- as.numeric(con$ui %*% d)
    step <- min(c(max.step, (s0 / -rate)[rate < 0]))
    starts[[length(starts) + 1]] <- setNames(center + frac * step * d, names(center))
  }
  starts
}

# f: criterion of theta (named vector). starts: list of named start vectors.
# Returns the best fit with diagnostics; $stage1 holds the value reached from each start.
fit.admissible <- function(f, starts, con, parscale = rep(1, length(starts[[1]])), smooth = TRUE,
                           reltol = 1e-10, maxit = 2000, polish.margin = 1e-2) {
  nm <- names(starts[[1]])
  fw <- function(p) {
    p <- setNames(p, nm)
    if (!in.domain(p, con)) return(WALL.VALUE)
    v <- f(p)
    if (!is.finite(v)) WALL.VALUE else v
  }
  runs <- lapply(starts, function(s) {
    if (!in.domain(s, con)) return(list(par = s, value = Inf, convergence = NA, counts = NA))
    o <- stats::optim(s, fw, method = "Nelder-Mead", control = list(parscale = parscale, reltol = reltol, maxit = maxit))
    # one restart from the end point (Nelder-Mead can stop early on a ridge)
    o2 <- stats::optim(o$par, fw, method = "Nelder-Mead", control = list(parscale = parscale, reltol = reltol, maxit = maxit))
    if (o2$value <= o$value) o2 else o
  })
  values <- vapply(runs, `[[`, 0, "value")
  best <- runs[[which.min(values)]]
  par <- setNames(best$par, nm); value <- best$value; method <- "Nelder-Mead"; convergence <- best$convergence
  # polish along the boundary when a constraint is within reach
  # slack in comparable units: gross return for the floor, parscale units for the box
  row.scale <- ifelse(con$kind == "floor", 1, parscale[match(sub("\\.(lower|upper)$", "", con$kind), nm)])
  if (min(slack(par, con) / row.scale) < polish.margin) {
    # the barrier needs a strictly interior start: move par slightly towards a
    # strictly interior start if it sits exactly on the boundary
    inner <- par
    if (any(slack(inner, con) <= 0)) {
      anchor <- Find(function(s) all(slack(s, con) > 0), starts)
      if (!is.null(anchor)) inner <- anchor + 0.999 * (par - anchor)
    }
    if (all(slack(inner, con) > 0)) {
      grad <- if (smooth) function(p) {
        h <- 1e-6 * parscale
        vapply(seq_along(p), function(k) { e <- numeric(length(p)); e[k] <- h[k]; (f(setNames(p + e, nm)) - f(setNames(p - e, nm))) / (2 * h[k]) }, 0)
      } else NULL
      co <- tryCatch(stats::constrOptim(inner, function(p) { v <- f(setNames(p, nm)); if (is.finite(v)) v else WALL.VALUE },
        grad = grad, ui = con$ui, ci = con$ci, method = if (smooth) "BFGS" else "Nelder-Mead",
        mu = 1e-8, control = list(parscale = parscale, reltol = reltol, maxit = maxit)), error = function(e) NULL)
      if (!is.null(co) && in.domain(co$par, con) && co$value < value) {
        par <- setNames(co$par, nm); value <- co$value; method <- "constrOptim"; convergence <- co$convergence
      }
    }
  }
  act <- active.constraints(par, con)
  mg <- min.gross.return(par, con)
  ends <- lapply(runs[is.finite(values) & values < WALL.VALUE], function(r) setNames(r$par, nm))
  list(par = par, value = value, convergence = convergence, method = method,
       at.bound = length(act) > 0, active = paste(act, collapse = ","),
       min.G = mg$value, min.G.date = mg$date, stage1 = values, stage1.par = lapply(runs, function(r) setNames(r$par, nm)), n.starts = length(starts),
       endpoints = distinct.points(c(list(par), ends), parscale))
}

# distinct points (distance in parscale units above tol), first occurrence kept
distinct.points <- function(points, parscale, tol = 1e-3) {
  keep <- list()
  for (p in points) {
    if (!any(vapply(keep, function(q) max(abs((p - q) / parscale)) < tol, TRUE))) keep[[length(keep) + 1]] <- p
  }
  keep
}

# One refit of the jackknife (or any refit that should report basin switches):
# starts at the full-sample estimate and at the distinct local optima found in
# the full-sample search. switched = the refit's estimate is a materially
# different point (more than switch.tol apart in parscale units) with a lower
# criterion (by more than 1e-6 relative) than the point reached from the
# full-sample estimate: the refit left that basin. A polish that only refines
# the same point is not a switch.
refit.admissible <- function(f, par.hat, endpoints, con, parscale, smooth = TRUE, switch.tol = 1e-2) {
  # par.hat stays first (it lies in Theta: the domain is the same in every refit)
  starts <- distinct.points(c(list(par.hat), endpoints), parscale)
  r <- fit.admissible(f, starts, con, parscale, smooth)
  warm <- r$stage1.par[[1]]
  r$switched <- is.finite(r$stage1[1]) && (r$stage1[1] - r$value) > 1e-6 * abs(r$value) &&
    max(abs((r$par - warm) / parscale)) > switch.tol
  r
}

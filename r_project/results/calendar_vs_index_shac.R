# WP1.2 Diagnostic: Calendar-lag Bartlett SHAC K(|y_v - y_w|/S) vs
# index-lag Bartlett HAC K(|v - w|/S) on empirical Preqin BO and VC panels
# (given missing early vintage years: 1989 in BO; 1983-84 in VC L2; 1988, 1989, 1991 in VC U).
# Run from r_project/:  LANG=en_US.UTF-8 Rscript results/calendar_vs_index_shac.R
# Writes results/figures_bovc_us_g1/calendar_vs_index_shac.csv.
suppressMessages({
  source("simulation/mc_driver.R")
  if (!exists("fixedb.crit.t")) source("helper/fixedb.R")
})

pub <- mc.public.factors()
fc <- unique(c("RF", setdiff(colnames(pub), "Date")))

load.units <- function(w, vyp, type) {
  d <- read.csv(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s%s_North America_NAV.csv", w, if (vyp) "_VYP" else ""), stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date)
  d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  d <- pad.all.units(d, "first.cf", max(d$Date))
  d <- merge(d, pub, by = "Date", all.x = TRUE)
  d <- d[d$type == type & d$Vintage <= 2010, ]
  d$Fund.ID <- droplevels(as.factor(as.character(d$Fund.ID)))
  list(rows = d, fp = prepare.fast(d, fc))
}

se_both_lags <- function(est, fp, par, max.month) {
  if (est$criterion == "L2") {
    units <- lapply(seq_len(fp$n.units), function(j) fast.subset(fp, seq_len(fp$n.units) == j))
    loss <- function(p, u) fast.unit.errors(p, u, max.month)^2
    vint_cal <- fp$vintage
  } else {
    w <- u.unit.weights(fp$group.n, est$vintage.weights)
    keep <- names(w)[!is.na(w)]
    units <- lapply(keep, function(k) fast.subset(fp, fp$key == k))
    names(units) <- keep
    loss <- function(p, u) w[[u$key[1]]] * fast.u.criterion(fast.unit.errors(p, u, max.month), u, mc.u.sizes(est, u), vintage.weights = est$vintage.weights)
    vint_cal <- vapply(units, function(u) u$vintage[1], 0)
  }
  ord <- order(vint_cal)
  units <- units[ord]; vint_cal <- vint_cal[ord]; vint_idx <- seq_along(vint_cal)
  list(
    cal = sandwich.covariance(units, par, loss, vint_cal, bw = "llsw", delta = 1e-4),
    idx = sandwich.covariance(units, par, loss, vint_idx, bw = "llsw", delta = 1e-4)
  )
}

ar_both_lags <- function(fp, m, alpha = 0.05) {
  br <- feasible.interval(mc.domain(fp, "MKT"))
  z <- fast.horizon.scale(MC.ESTIMATORS$IV$reference, fp, m)
  V <- fp$n.units; S <- hac.bandwidth.llsw(V); b.rat <- S / V; crit <- fixedb.crit.t(alpha, b.rat)
  v_cal <- fp$vintage; v_idx <- match(v_cal, sort(unique(v_cal)))
  run_ar <- function(vvec) {
    K <- outer(vvec, vvec, bartlett.weight, bw = S)
    margin <- function(beta) {
      e <- fast.unit.errors(c(MKT = beta), fp, 0)
      if (any(!is.finite(e))) return(NA_real_)
      s <- z * e; gbar <- mean(s); sc <- s - gbar
      omega <- sum(K * outer(sc, sc)) / V
      if (!is.finite(omega) || omega <= 0) return(if (abs(gbar) < 1e-12) -crit else Inf)
      abs(gbar) / sqrt(omega / V) - crit
    }
    grid <- seq(br[1], br[2], length.out = 401L); mv <- vapply(grid, margin, 0)
    in.set <- is.finite(mv) & mv <= 0
    runs <- rle(in.set); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1L
    ivs <- lapply(which(runs$values), function(r) {
      i1 <- starts[r]; i2 <- ends[r]
      lo <- if (i1 == 1L) br[1] else uniroot(margin, c(grid[i1 - 1L], grid[i1]), tol = 1e-9)$root
      hi <- if (i2 == length(grid)) br[2] else uniroot(margin, c(grid[i2], grid[i2 + 1L]), tol = 1e-9)$root
      c(lo, hi)
    })
    list(lo = ivs[[1]][1], hi = ivs[[length(ivs)]][2], p1 = fixedb.q1.pvalue((margin(1) + crit)^2, b.rat))
  }
  list(cal = run_ar(v_cal), idx = run_ar(v_idx), V = V, S = S, b = b.rat)
}

rows <- list()
for (ty in c("BO", "VC")) {
  for (w in c("EW", "FW")) {
    vyp <- load.units(w, TRUE, ty)$fp
    fund <- load.units(w, FALSE, ty)$fp
    estU <- if (w == "FW") MC.ESTIMATORS$U.eq.fw else MC.ESTIMATORS$U.eq
    estL <- MC.ESTIMATORS$L2.vyp
    for (m in c(0, 60, 120, 180, 360)) {
      fU <- mc.fit(estU, fund, "MKT", m); sU <- se_both_lags(estU, fund, c(MKT = fU$MKT), m)
      fL <- mc.fit(estL, vyp, "MKT", m);  sL <- se_both_lags(estL, vyp, c(MKT = fL$MKT), m)
      bU <- sU$cal$hac.b
      rows[[length(rows) + 1L]] <- data.frame(
        type = ty, weighting = w, estimator = "U.eq", H = m / 12, V = sU$cal$hac.V, S = sU$cal$hac.S, b = bU,
        est = fU$MKT, se_cal = sU$cal$SE[[1]], se_idx = sU$idx$SE[[1]], se_diff = sU$idx$SE[[1]] - sU$cal$SE[[1]],
        p1_cal = fixedb.q1.pvalue(((fU$MKT - 1)/sU$cal$SE[[1]])^2, bU),
        p1_idx = fixedb.q1.pvalue(((fU$MKT - 1)/sU$idx$SE[[1]])^2, bU),
        p0_cal = fixedb.q1.pvalue((fU$MKT/sU$cal$SE[[1]])^2, bU),
        p0_idx = fixedb.q1.pvalue((fU$MKT/sU$idx$SE[[1]])^2, bU),
        neff = sU$cal$neff.units, ar_lo_cal = NA_real_, ar_lo_idx = NA_real_, stringsAsFactors = FALSE
      )
      bL <- sL$cal$hac.b
      rows[[length(rows) + 1L]] <- data.frame(
        type = ty, weighting = w, estimator = "L2.vyp", H = m / 12, V = sL$cal$hac.V, S = sL$cal$hac.S, b = bL,
        est = fL$MKT, se_cal = sL$cal$SE[[1]], se_idx = sL$idx$SE[[1]], se_diff = sL$idx$SE[[1]] - sL$cal$SE[[1]],
        p1_cal = fixedb.q1.pvalue(((fL$MKT - 1)/sL$cal$SE[[1]])^2, bL),
        p1_idx = fixedb.q1.pvalue(((fL$MKT - 1)/sL$idx$SE[[1]])^2, bL),
        p0_cal = fixedb.q1.pvalue((fL$MKT/sL$cal$SE[[1]])^2, bL),
        p0_idx = fixedb.q1.pvalue((fL$MKT/sL$idx$SE[[1]])^2, bL),
        neff = sL$cal$neff.units, ar_lo_cal = NA_real_, ar_lo_idx = NA_real_, stringsAsFactors = FALSE
      )
    }
    ar0 <- ar_both_lags(vyp, 0)
    rows[[length(rows) + 1L]] <- data.frame(
      type = ty, weighting = w, estimator = "AR.set", H = 0, V = ar0$V, S = ar0$S, b = ar0$b,
      est = NA_real_, se_cal = NA_real_, se_idx = NA_real_, se_diff = NA_real_,
      p1_cal = ar0$cal$p1, p1_idx = ar0$idx$p1, p0_cal = NA_real_, p0_idx = NA_real_,
      neff = NA_real_, ar_lo_cal = ar0$cal$lo, ar_lo_idx = ar0$idx$lo, stringsAsFactors = FALSE
    )
  }
}
out <- do.call(rbind, rows)
write.csv(out, "results/figures_bovc_us_g1/calendar_vs_index_shac.csv", row.names = FALSE)
cat("wrote results/figures_bovc_us_g1/calendar_vs_index_shac.csv\n")
print(out[out$H == 0, c("type","weighting","estimator","V","est","se_cal","se_idx","se_diff","p1_cal","p1_idx","neff","ar_lo_cal","ar_lo_idx")], row.names = FALSE)

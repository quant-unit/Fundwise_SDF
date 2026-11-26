# Empirical IV moment comparison, identification-robust (Anderson-Rubin-style)
# confidence sets, matched-vintage L2 support, and observed vs. predicted L2-U
# gaps (revision plan, Section L; reconciled Nov 22 fixes).
# Run from r_project/:  LANG=en_US.UTF-8 Rscript results/empirical_iv.R
#
# IV: H = 0 mean pricing moment V^-1 sum_v PV_v(beta) = 0 (and H > 0 path
# instrument z_v = A_{v,H}(1)/A_{v,0}(1)) on the r5 admissible domain.
# Reports:
#   1. All detected roots on the fixed full-sample domain, plus leave-one-vintage-out
#      (LOO) zero-root failure counts (never silently substituting min.square into
#      jackknife SEs).
#   2. Nominal 95% confidence sets obtained by inverting the centered Bartlett
#      fixed-b moment test at every candidate beta in the admissible interval
#      (with refined interval endpoints via uniroot; records centering convention,
#      S, V, b, critical value, and all acceptance intervals). Explicitly qualified
#      as identification-robust / calibrated-by-design (still over-rejects under
#      persistent finite-V shocks; rates read from results/mc_step4/inference_checks.csv,
#      written by simulation/mc_inference_checks.R).
#   3. Both full-sample L2/IV and dynamically derived matched-vintage (n_v >= 2)
#      L2/IV so L2 - U_{v,H} gaps are compared on identical vintage support.
# Writes to results/figures_bovc_us_g1/:
#   iv_comparison.tex, iv_comparison_data.csv, iv_comparison_rows.csv,
#   iv_estimates_all_horizons.csv, iv_ar_confidence_sets.csv
source("simulation/mc_driver.R")
if (!exists("fixedb.q1.pvalue")) source("helper/fixedb.R")
out.dir <- "results/figures_bovc_us_g1"
H.ALL <- c(0, 1, 30, 60, 120, 150, 180, 210, 240, 270, 300, 330, 360)
H.TAB <- c(0, 60, 180, 360)

pub <- local({
  p <- read.csv("empirical/data_prepared_2026/q_factors.csv")
  p$Date <- as.Date(p$Date)
  p$Date <- as.Date(format(p$Date, "%Y-%m-01"))
  p$Date <- seq.Date(p$Date[1], by = "month", length.out = nrow(p))
  p$Date <- as.Date(format(p$Date + 32, "%Y-%m-01")) - 1
  p$Alpha <- 1
  p
})
fc <- unique(c("RF", setdiff(colnames(pub), "Date")))

# padded rows and fast structure of one file, type (as estim_model_optimized.R)
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

cells <- expand.grid(type = c("BO", "VC"), w = c("EW", "FW"), stringsAsFactors = FALSE)
U <- list()
for (i in seq_len(nrow(cells))) {
  k <- paste(cells$type[i], cells$w[i])
  vyp <- load.units(cells$w[i], TRUE, cells$type[i])
  fund <- load.units(cells$w[i], FALSE, cells$type[i])
  # Dynamically derive eligible U_{v,H} vintages (n_v >= 2 funds) from fund file
  v.cnt <- table(fund$fp$vintage)
  elig.v <- sort(as.integer(names(v.cnt)[v.cnt >= 2L]))
  vyp.matched.fp <- fast.subset(vyp$fp, vyp$fp$vintage %in% elig.v)
  U[[k]] <- list(vyp = vyp, fund = fund, elig.v = elig.v, vyp.matched = vyp.matched.fp)
}
est.L2 <- MC.ESTIMATORS$L2.vyp
est.U <- function(w) if (w == "FW") MC.ESTIMATORS$U.eq.fw else MC.ESTIMATORS$U.eq

# All roots of the IV moment on a bracket br (default: full-sample domain)
iv.roots <- function(fp, m, br = NULL) {
  z <- fast.horizon.scale(MC.ESTIMATORS$IV$reference, fp, m)
  g <- function(b) { e <- fast.unit.errors(c(MKT = b), fp, 0); if (any(!is.finite(e))) NA_real_ else mean(z * e) }
  if (is.null(br)) br <- feasible.interval(mc.domain(fp, "MKT"))
  grid <- seq(br[1], br[2], length.out = 401)
  gv <- vapply(grid, g, 0)
  ch <- which(is.finite(gv[-1]) & is.finite(gv[-length(gv)]) & sign(gv[-1]) != sign(gv[-length(gv)]))
  vapply(ch, function(i) {
    if (gv[i] == 0) grid[i] else if (gv[i + 1] == 0) grid[i + 1] else uniroot(g, grid[c(i, i + 1)], tol = 1e-12)$root
  }, 0)
}

# Centered Bartlett fixed-b moment-test inversion (Anderson-Rubin-style 95% set).
# Tests H0: E[z_v PV_v(beta)] = 0 using centered scores s_tilde_v = s_v - mean(s),
# Bartlett kernel with S = ceil(1.3 * sqrt(V)) over vintage-year distances, and
# two-sided fixed-b critical value c_0.05(b = S/V). Refines all acceptance
# interval endpoints via uniroot.
iv.ar.set <- function(fp, m, alpha = 0.05, br = NULL, n.grid = 401L) {
  if (is.null(br)) br <- feasible.interval(mc.domain(fp, "MKT"))
  z <- fast.horizon.scale(MC.ESTIMATORS$IV$reference, fp, m)
  V <- fp$n.units
  S <- hac.bandwidth.llsw(V)
  b.rat <- S / V
  K <- outer(fp$vintage, fp$vintage, bartlett.weight, bw = S)
  crit <- fixedb.crit.t(alpha, b.rat)

  # Signed margin: <= 0 means inside the acceptance region (|t| <= crit)
  margin <- function(beta) {
    e <- fast.unit.errors(c(MKT = beta), fp, 0)
    if (any(!is.finite(e))) return(NA_real_)
    s <- z * e
    gbar <- mean(s)
    sc <- s - gbar # centered scores
    omega <- sum(K * outer(sc, sc)) / V
    if (!is.finite(omega) || omega <= 0) return(if (abs(gbar) < 1e-12) -crit else Inf)
    abs(gbar) / sqrt(omega / V) - crit
  }

  grid <- seq(br[1], br[2], length.out = n.grid)
  mv <- vapply(grid, margin, 0)
  in.set <- is.finite(mv) & mv <= 0

  # Find every contiguous acceptance segment and refine interior boundary crossings
  intervals <- list()
  if (any(in.set)) {
    runs <- rle(in.set)
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    for (r in which(runs$values)) {
      i1 <- starts[r]; i2 <- ends[r]
      lo <- if (i1 == 1L) br[1] else {
        if (is.finite(mv[i1 - 1L]) && is.finite(mv[i1]) && mv[i1 - 1L] > 0 && mv[i1] <= 0) {
          uniroot(margin, c(grid[i1 - 1L], grid[i1]), tol = 1e-9)$root
        } else grid[i1]
      }
      hi <- if (i2 == length(grid)) br[2] else {
        if (is.finite(mv[i2]) && is.finite(mv[i2 + 1L]) && mv[i2] <= 0 && mv[i2 + 1L] > 0) {
          uniroot(margin, c(grid[i2], grid[i2 + 1L]), tol = 1e-9)$root
        } else grid[i2]
      }
      intervals[[length(intervals) + 1L]] <- c(lo, hi)
    }
  }
  # Also check if beta = 1 is covered directly
  m.at.1 <- margin(1)
  fmt.iv <- if (!length(intervals)) "empty" else paste(vapply(intervals, function(iv) sprintf("[%.3f, %.3f]", iv[1], iv[2]), ""), collapse = " U ")
  list(
    V = V, S = S, b = b.rat, crit = crit, centering = "centered",
    br.lo = br[1], br.hi = br[2], n.intervals = length(intervals),
    lo = if (length(intervals)) intervals[[1]][1] else NA_real_,
    hi = if (length(intervals)) intervals[[length(intervals)]][2] else NA_real_,
    intervals.str = fmt.iv, covers.beta1 = is.finite(m.at.1) && m.at.1 <= 0,
    t.at.beta1 = if (is.finite(m.at.1)) m.at.1 + crit else NA_real_,
    p.at.beta1 = if (is.finite(m.at.1)) fixedb.q1.pvalue((m.at.1 + crit)^2, b.rat) else NA_real_
  )
}

pfb <- function(b, se, hb) if (is.finite(se) && is.finite(hb)) fixedb.q1.pvalue(((b - 1) / se)^2, hb) else NA_real_

# ---- IV, all horizons (full sample and matched-vintage support) -------------------
iv <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  k <- paste(cells$type[i], cells$w[i]); fp <- U[[k]]$vyp$fp
  do.call(rbind, lapply(H.ALL, function(m) {
    r <- mc.fit.iv(MC.ESTIMATORS$IV, fp, "MKT", m)
    rt <- iv.roots(fp, m)
    data.frame(type = cells$type[i], weighting = cells$w[i], max.month = m, H = m / 12,
      beta = r$MKT, se = r$SE.MKT, p.beta1 = pfb(r$MKT, r$SE.MKT, r$hac.b),
      n.roots = r$n.roots, roots = paste(sprintf("%.3f", rt), collapse = "; "),
      at.bound = r$at.bound, active = r$active, V = r$hac.V, S = r$hac.S, b = r$hac.b,
      stringsAsFactors = FALSE)
  }))
}))
write.csv(iv, file.path(out.dir, "iv_estimates_all_horizons.csv"), row.names = FALSE)

# ---- check 1: IV against an independent computation from the pipeline's unit errors ----
chk1 <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  k <- paste(cells$type[i], cells$w[i]); prep <- prepare.units(U[[k]]$vyp$rows, fc)
  do.call(rbind, lapply(c(0, 60, 180), function(m) {
    r <- iv[iv$type == cells$type[i] & iv$weighting == cells$w[i] & iv$max.month == m, ]
    pe <- function(b, mm) vapply(prep$units, unit.horizon.error, 0, max.month = mm, par0 = c(RF = 1, MKT = b), sdf.model = "linear")
    zz <- pe(1, m) / pe(1, 0); gb <- function(b) mean(zz * pe(b, 0))
    yrs <- vapply(prep$units, `[[`, 0, "vintage"); V <- length(yrs); S <- ceiling(1.3 * sqrt(V)); s <- zz * pe(r$beta, 0)
    om <- 0; for (a in seq_len(V)) for (j in seq_len(V)) om <- om + max(0, 1 - abs(yrs[a] - yrs[j]) / S) * s[a] * s[j]
    Gd <- (gb(r$beta + 1e-4) - gb(r$beta - 1e-4)) / 2e-4
    data.frame(cell = k, H = m / 12, moment.at.estimate = gb(r$beta), se.ratio = r$se / sqrt(om / V / (V * Gd^2)))
  }))
}))
cat("check 1 (IV against the pipeline's unit errors, explicit kernel sum): max |moment| at the estimate", signif(max(abs(chk1$moment.at.estimate)), 2),
    "| max |SE ratio - 1|", signif(max(abs(chk1$se.ratio - 1), na.rm = TRUE), 2), "\n")
stopifnot(max(abs(chk1$moment.at.estimate)) < 1e-8, max(abs(chk1$se.ratio - 1), na.rm = TRUE) < 1e-6)

# ---- L2 and U from the table data, refitted as a check; matched-vintage support; AR sets ----
tab.l2 <- rbind(read.csv(file.path(out.dir, "empirical_BO_US_MKT_data.csv")), read.csv(file.path(out.dir, "empirical_VC_US_MKT_data.csv")))
tab.l2 <- tab.l2[tab.l2$Factor == "MKT" & tab.l2$method == "Asymptotic", ]
tab.u <- rbind(read.csv(file.path(out.dir, "upair_BO_US_data.csv")), read.csv(file.path(out.dir, "upair_VC_US_data.csv")))
tab.u <- tab.u[tab.u$model == "MKT" & tab.u$parameter == "MKT" & tab.u$vintage.weights == "equal", ]

# Point estimate only on fixed full-sample bracket br; IV returns NA if no root
fit.b <- function(est, fp, m, br = NULL) {
  if (is.null(br)) br <- feasible.interval(mc.domain(fp, "MKT"))
  if (est$criterion == "IV") {
    rt <- iv.roots(fp, m, br = br)
    return(if (length(rt)) rt[1] else NA_real_)
  }
  f <- mc.criterion(est, fp, m)
  stats::optim(c(MKT = 1), function(p) f(setNames(p, "MKT")), method = "Brent", lower = br[1], upper = br[2])$par[[1]]
}

res <- list(); ar.rows <- list(); repro <- 0
for (i in seq_len(nrow(cells))) {
  ty <- cells$type[i]; w <- cells$w[i]; k <- paste(ty, w)
  fv <- U[[k]]$vyp$fp; ff <- U[[k]]$fund$fp; fvm <- U[[k]]$vyp.matched
  br.v <- feasible.interval(mc.domain(fv, "MKT"))
  br.f <- feasible.interval(mc.domain(ff, "MKT"))
  br.vm <- feasible.interval(mc.domain(fvm, "MKT"))
  for (m in H.TAB) {
    l2 <- tab.l2[tab.l2$Type == ty & tab.l2$weighting == w & tab.l2$max.month == m, ]
    u <- tab.u[tab.u$Type == ty & tab.u$weighting == w & tab.u$max.month == m, ]
    v <- iv[iv$type == ty & iv$weighting == w & iv$max.month == m, ]
    stopifnot(nrow(l2) == 1, nrow(u) == 1, nrow(v) == 1)
    b.l2 <- fit.b(est.L2, fv, m, br = br.v)
    b.u <- fit.b(est.U(w), ff, m, br = br.f)
    repro <- max(repro, abs(b.l2 - l2$MKT), abs(b.u - u$estimate))

    # Matched-vintage L2 and IV (identical vintage set as U_{v,H}, n_v >= 2)
    l2.m <- mc.fit(est.L2, fvm, "MKT", m)
    iv.m <- mc.fit.iv(MC.ESTIMATORS$IV, fvm, "MKT", m)
    rt.m <- iv.roots(fvm, m, br = br.vm)

    # Centered fixed-b Bartlett AR 95% confidence sets (full and matched)
    ar.full <- iv.ar.set(fv, m, alpha = 0.05, br = br.v)
    ar.match <- iv.ar.set(fvm, m, alpha = 0.05, br = br.vm)
    ar.rows[[length(ar.rows) + 1L]] <- data.frame(type = ty, weighting = w, support = "full", max.month = m, H = m / 12,
      n.roots = v$n.roots, roots = v$roots, ar.full, stringsAsFactors = FALSE)
    ar.rows[[length(ar.rows) + 1L]] <- data.frame(type = ty, weighting = w, support = "matched", max.month = m, H = m / 12,
      n.roots = iv.m$n.roots, roots = paste(sprintf("%.3f", rt.m), collapse = "; "), ar.match, stringsAsFactors = FALSE)

    # Leave-one-vintage-out jackknife on fixed full-sample domains
    vs <- sort(unique(fv$vintage))
    jk <- t(vapply(vs, function(x) c(
      fit.b(est.L2, fast.subset(fv, fv$vintage != x), m, br = br.v),
      fit.b(est.U(w), fast.subset(ff, ff$vintage != x), m, br = br.f),
      fit.b(MC.ESTIMATORS$IV, fast.subset(fv, fv$vintage != x), m, br = br.v)), numeric(3)))
    vs.m <- sort(unique(fvm$vintage))
    jk.m <- t(vapply(vs.m, function(x) c(
      fit.b(est.L2, fast.subset(fvm, fvm$vintage != x), m, br = br.vm),
      fit.b(est.U(w), fast.subset(ff, ff$vintage != x), m, br = br.f),
      fit.b(MC.ESTIMATORS$IV, fast.subset(fvm, fvm$vintage != x), m, br = br.vm)), numeric(3)))
    jse <- function(x) if (anyNA(x)) NA_real_ else sqrt((length(x) - 1) / length(x) * sum((x - mean(x))^2))
    n.no.root <- sum(is.na(jk[, 3]))
    n.no.root.m <- sum(is.na(jk.m[, 3]))

    res[[length(res) + 1L]] <- data.frame(
      type = ty, weighting = w, max.month = m, H = m / 12,
      V.full = fv$n.units, V.matched = fvm$n.units,
      l2 = l2$MKT, l2.se = if (isTRUE(l2$at.bound)) NA else l2$SE.MKT, l2.bound = isTRUE(l2$at.bound),
      l2.matched = l2.m$MKT, l2.matched.se = if (isTRUE(l2.m$at.bound)) NA else l2.m$SE.MKT, l2.matched.bound = isTRUE(l2.m$at.bound),
      u = u$estimate, u.se = if (isTRUE(u$at.bound)) NA else u$SE.A, u.bound = isTRUE(u$at.bound),
      iv = v$beta, iv.se = v$se, iv.bound = v$at.bound, iv.roots = v$n.roots, iv.root.list = v$roots,
      iv.loo.no.root = n.no.root,
      iv.matched = iv.m$MKT, iv.matched.se = iv.m$SE.MKT, iv.matched.roots = iv.m$n.roots,
      iv.matched.root.list = paste(sprintf("%.3f", rt.m), collapse = "; "), iv.matched.loo.no.root = n.no.root.m,
      ar.lo = ar.full$lo, ar.hi = ar.full$hi, ar.set = ar.full$intervals.str, ar.crit = ar.full$crit, ar.p1 = ar.full$p.at.beta1,
      ar.matched.lo = ar.match$lo, ar.matched.hi = ar.match$hi, ar.matched.set = ar.match$intervals.str, ar.matched.crit = ar.match$crit,
      gap.l2u = l2$MKT - u$estimate, gap.l2u.jse = jse(jk[, 1] - jk[, 2]),
      gap.l2u.matched = l2.m$MKT - u$estimate, gap.l2u.matched.jse = jse(jk.m[, 1] - jk.m[, 2]),
      gap.uiv = u$estimate - v$beta, gap.uiv.jse = jse(jk[, 2] - jk[, 3]),
      jk.V = length(vs), stringsAsFactors = FALSE
    )
  }
}
res <- do.call(rbind, res)
ar.df <- do.call(rbind, ar.rows)
write.csv(ar.df, file.path(out.dir, "iv_ar_confidence_sets.csv"), row.names = FALSE)
cat("check 2 (L2 and U refitted with the fast criterion = the G1/G2 table data): max |difference|", signif(repro, 2), "\n")
stopifnot(repro < 1e-5)

# ---- predicted L2 - U gap from the data-shaped Monte Carlo (paired over replications) ----
pred <- function(type, w, m) {
  d <- if (w == "EW") paste0("core_", type) else paste0("s_core_", type)
  r <- mc.read(d); r <- r[r$model == "MKT" & r$max.month == m, ]
  e1 <- if (w == "EW") "L2.vyp" else "L2.vyp.fw"; e2 <- if (w == "EW") "U.eq" else "U.eq.fw"
  a <- r[r$estimator == e1, c("rep", "MKT")]; b <- r[r$estimator == e2, c("rep", "MKT")]
  if (!nrow(a) || !nrow(b)) return(c(NA, NA))
  x <- merge(a, b, by = "rep"); g <- x$MKT.x - x$MKT.y
  c(mean(g), sd(g) / sqrt(length(g)))
}
pg <- t(vapply(seq_len(nrow(res)), function(j) pred(res$type[j], res$weighting[j], res$max.month[j]), numeric(2)))
res$pred.l2u <- pg[, 1]; res$pred.l2u.mcse <- pg[, 2]
write.csv(res, file.path(out.dir, "iv_comparison_rows.csv"), row.names = FALSE)

# ---- table ---------------------------------------------------------------------------
cells.rec <- list()
rec <- function(key, column, value, printed) {
  cells.rec[[length(cells.rec) + 1L]] <<- data.frame(key, column = column, value = value, printed = printed, stringsAsFactors = FALSE)
  printed
}
f3 <- function(x) if (is.na(x)) "--" else formatC(round(x, 3) + 0, format = "f", digits = 3)
es <- function(key, col, b, se, bound) {
  sprintf("%s%s (%s)", rec(key, col, b, f3(b)), if (isTRUE(bound)) "$^\\dagger$" else "",
          rec(key, paste0(col, ".se"), if (isTRUE(bound)) NA else se, f3(if (isTRUE(bound)) NA else se)))
}
lines <- character(0)
for (ty in c("BO", "VC")) for (w in c("EW", "FW")) {
  r0 <- res[res$type == ty & res$weighting == w & res$max.month == 0, ]
  sub.note <- if (r0$V.full == r0$V.matched) {
    sprintf("$V=%d$", r0$V.full)
  } else {
    sprintf("$V=%d$ full / %d matched ($n_v\\ge2$)", r0$V.full, r0$V.matched)
  }
  lines <- c(lines, sprintf("\\multicolumn{8}{l}{\\textit{%s, %s (%s)}} \\\\",
    c(BO = "Buyout", VC = "Venture capital")[[ty]], c(EW = "equally weighted", FW = "fund-size weighted")[[w]], sub.note))
  for (m in H.TAB) {
    r <- res[res$type == ty & res$weighting == w & res$max.month == m, ]
    key <- data.frame(type = ty, weighting = w, H = m / 12)
    l2.str <- es(key, "l2", r$l2, r$l2.se, r$l2.bound)
    u.str <- es(key, "u", r$u, r$u.se, r$u.bound)
    gap.str <- if (r$V.full == r$V.matched) {
      sprintf("%s (%s)", rec(key, "gap.l2u", r$gap.l2u, f3(r$gap.l2u)), rec(key, "gap.l2u.jse", r$gap.l2u.jse, f3(r$gap.l2u.jse)))
    } else {
      sprintf("%s [%s]",
        rec(key, "gap.l2u", r$gap.l2u, f3(r$gap.l2u)),
        rec(key, "gap.l2u.matched", r$gap.l2u.matched, f3(r$gap.l2u.matched)))
    }
    pred.str <- if (is.finite(r$pred.l2u)) {
      sprintf("%s (%s)", rec(key, "pred.l2u", r$pred.l2u, f3(r$pred.l2u)), rec(key, "pred.l2u.mcse", r$pred.l2u.mcse, f3(r$pred.l2u.mcse)))
    } else "--"
    rt.str <- rec(key, "iv.root.list", NA_real_, if (nzchar(r$iv.root.list)) r$iv.root.list else "none")
    loo.str <- rec(key, "iv.loo.no.root", r$iv.loo.no.root, sprintf("%d/%d", r$iv.loo.no.root, r$V.full))
    ar.str <- rec(key, "ar.set", NA_real_, r$ar.set)
    lines <- c(lines, sprintf("\\quad %s yr & %s & %s & %s & %s & %s & %s & %s \\\\",
      format(m / 12), u.str, l2.str, gap.str, pred.str, rt.str, loo.str, ar.str))
  }
  if (!(ty == "VC" && w == "FW")) lines <- c(lines, "\\addlinespace")
}
# rejection rates of the AR moment test at beta = 1 in the Monte Carlo (simulation/mc_inference_checks.R)
ic <- read.csv("results/mc_step4/inference_checks.csv", stringsAsFactors = FALSE); ic <- ic[ic$check == "ar_replay", ]
arp <- function(d) { x <- ic$reject[ic$design == d]; stopifnot(length(x) == 1); sprintf("%.1f", 100 * x) }
tex <- c(
  sprintf("%% Auto-generated by r_project/results/empirical_iv.R at commit %s%s",
          system("git rev-parse --short HEAD", intern = TRUE),
          if (length(system("git status --porcelain --untracked-files=no", intern = TRUE))) " (uncommitted changes)" else ""),
  "\\begin{table}[htbp]", "\\centering", "\\footnotesize",
  "\\caption{Single-factor $\\hat\\beta_{MKT}$: $U_{v,H}$, L2, L2$-U_{v,H}$ against the expected core gap, and identification-robust IV moment sets}",
  "\\label{tab:iv_comparison}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{l rr rr r c c}", "\\toprule",
  "Horizon & $U_{v,H}$, $q_v=1$ & L2 (full) & L2$-U_{v,H}$ & Expected gap & IV roots & LOO no-root & Nominal 95\\% AR set \\\\", "\\midrule",
  lines, "\\bottomrule", "\\end{tabular}%", "}", "\\par\\vspace{2pt}",
  paste0("{\\scriptsize\\textit{Notes:} North America, vintages up to 2010, single-factor $\\hat\\beta_{MKT}$ (Bartlett sandwich SE, $S=\\lceil1.3\\sqrt V\\rceil$). ",
    "L2$-U_{v,H}$: for Buyout ($V=25$ for both estimators), full-sample difference with leave-one-vintage-out (LOO) jackknife SE in parentheses (fixed full-sample domain); for VC ($V=27$ full vs.\\ $V=23$ eligible $n_v\\ge2$ vintages), full-sample difference followed in brackets by the matched-vintage difference on the same 23 vintages (matched L2 estimates, jackknife SEs, and matched AR sets are recorded in \\texttt{iv\\_comparison\\_rows.csv}). ",
    "Expected gap: paired mean L2 minus $U_{v,H}$ difference in the data-shaped core Monte Carlo without vintage-level effects (EW: \\texttt{core\\_BO/VC}; FW: \\texttt{s\\_core\\_BO/VC} at 0 and 15 yr), MC SE in parentheses. ",
    "IV roots: all zeros of $V^{-1}\\sum_v z_v PV_v(\\beta)=0$ ($z_v=A_{v,H}(1)/A_{v,0}(1)$, $z_v=1$ at $H=0$) on the admissible domain; LOO no-root: share of delete-one-vintage refits on the fixed full-sample domain that lose every root (zero-root refits are never replaced by squared-moment minimizers). ",
    "Nominal 95\\% AR set: acceptance region from inverting the centered Bartlett fixed-$b$ moment test (scores $s_v(\\beta)-\\bar s(\\beta)$, $S=\\lceil1.3\\sqrt V\\rceil$, $b=S/V$, critical values 2.790 for $V=25$ and 2.722 for $V=27$), avoiding root selection and weak-slope division. ",
    sprintf("Calibration caveat: robustness to roots and coverage under finite-$V$ persistence are distinct properties---at true $\\beta=1$, this moment test rejects %s\\%%/%s\\%% in independent BO/VC cores, %s\\%%/%s\\%% under persistent M4 ($c/2c$, $\\rho_W=0.6$), and %s\\%%/%s\\%% ($c$) and %s\\%%/%s\\%% ($2c$) in Design~A$'$ ($V=25/100$), so the sets are reported as a nominal, calibrated-by-design diagnostic. $^\\dagger$: estimate on a bound of $\\Theta$, no SE.}",
      arp("core_BO"), arp("core_VC"), arp("m4p_shared"), arp("m4p_shared2"), arp("pd_V25_n20_c"), arp("pd_V100_n20_c"), arp("pd_V25_n20_2c"), arp("pd_V100_n20_2c"))),
  "\\end{table}"
)
writeLines(tex, file.path(out.dir, "iv_comparison.tex"))
cr <- do.call(rbind, cells.rec)
write.csv(cr, file.path(out.dir, "iv_comparison_data.csv"), row.names = FALSE)
cat("wrote", file.path(out.dir, "iv_comparison.tex"), "with", nrow(cr), "recorded numbers\n")

# Empirical IV estimates next to L2 and U_{v,H}, and the observed against the
# predicted L2-U gap (revision plan, log entries "empirical IV estimates and the
# predicted L2--U gap: scope" and "...: results"; Section L, weak link 2).
# Run from r_project/:  LANG=en_US.UTF-8 Rscript results/empirical_iv.R
#
# IV: mc.fit.iv() of the Step 4 Monte Carlo (simulation/mc_driver.R), the
# smallest root in beta of V^-1 sum_v z_v PV_v(beta) = 0, z_v = A_vH(1)/A_v0(1)
# (z_v = 1 at H = 0), on the r5 domain, Bartlett sandwich SE, fixed-b p at
# analysis; every root is listed. Input: the prepared main VYP files (North
# America, EW and FW), prepared as estim_model_optimized.R does.
# L2 and U (q_v = 1): the G1/G2 table data files; refitted here with the same
# fast criterion as a reproduction check and for leave-one-vintage-out
# jackknife SEs of the gaps L2 - U and U - IV (descriptive, no test).
# Predicted L2 - U gap: L2 bias minus U bias in the data-shaped Monte Carlo
# (EW: core_BO/core_VC at 0, 5, 15, 30 years; FW: s_core_BO/s_core_VC at 0 and
# 15 years), paired over replications.
# Writes results/figures_bovc_us_g1/iv_comparison.tex and iv_comparison_data.csv
# (every printed number) and iv_estimates_all_horizons.csv (the 52 IV fits).
source("simulation/mc_driver.R")
if (!exists("fixedb.q1.pvalue")) source("helper/fixedb.R")
out.dir <- "results/figures_bovc_us_g1"
H.ALL <- c(0, 1, 30, 60, 120, 150, 180, 210, 240, 270, 300, 330, 360) # the empirical table grid (1 = alias of 0)
H.TAB <- c(0, 60, 180, 360)

pub <- local({ p <- read.csv("empirical/data_prepared_2026/q_factors.csv"); p$Date <- as.Date(p$Date)
  p$Date <- as.Date(format(p$Date, "%Y-%m-01")); p$Date <- seq.Date(p$Date[1], by = "month", length.out = nrow(p))
  p$Date <- as.Date(format(p$Date + 32, "%Y-%m-01")) - 1; p$Alpha <- 1; p })
fc <- unique(c("RF", setdiff(colnames(pub), "Date")))
# padded rows and fast structure of one file, type (as estim_model_optimized.R)
load.units <- function(w, vyp, type) {
  d <- read.csv(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s%s_North America_NAV.csv", w, if (vyp) "_VYP" else ""), stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date); d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  d <- pad.all.units(d, "first.cf", max(d$Date))
  d <- merge(d, pub, by = "Date", all.x = TRUE); d <- d[d$type == type & d$Vintage <= 2010, ]
  d$Fund.ID <- droplevels(as.factor(as.character(d$Fund.ID)))
  list(rows = d, fp = prepare.fast(d, fc))
}
cells <- expand.grid(type = c("BO", "VC"), w = c("EW", "FW"), stringsAsFactors = FALSE)
U <- list(); for (i in seq_len(nrow(cells))) { k <- paste(cells$type[i], cells$w[i])
  U[[k]] <- list(vyp = load.units(cells$w[i], TRUE, cells$type[i]), fund = load.units(cells$w[i], FALSE, cells$type[i])) }
est.L2 <- MC.ESTIMATORS$L2.vyp
est.U <- function(w) if (w == "FW") MC.ESTIMATORS$U.eq.fw else MC.ESTIMATORS$U.eq

# every root of the IV moment on the domain (same grid and moment as mc.fit.iv)
iv.roots <- function(fp, m) {
  z <- fast.horizon.scale(MC.ESTIMATORS$IV$reference, fp, m)
  g <- function(b) { e <- fast.unit.errors(c(MKT = b), fp, 0); if (any(!is.finite(e))) NA_real_ else mean(z * e) }
  br <- feasible.interval(mc.domain(fp, "MKT")); grid <- seq(br[1], br[2], length.out = 201); gv <- vapply(grid, g, 0)
  ch <- which(is.finite(gv[-1]) & is.finite(gv[-length(gv)]) & sign(gv[-1]) != sign(gv[-length(gv)]))
  vapply(ch, function(i) if (gv[i] == 0) grid[i] else uniroot(g, grid[c(i, i + 1)], tol = 1e-12)$root, 0)
}
pfb <- function(b, se, hb) if (is.finite(se) && is.finite(hb)) fixedb.q1.pvalue(((b - 1) / se)^2, hb) else NA_real_

# ---- IV, all horizons ------------------------------------------------------------
iv <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) { k <- paste(cells$type[i], cells$w[i]); fp <- U[[k]]$vyp$fp
  do.call(rbind, lapply(H.ALL, function(m) { r <- mc.fit.iv(MC.ESTIMATORS$IV, fp, "MKT", m); rt <- iv.roots(fp, m)
    data.frame(type = cells$type[i], weighting = cells$w[i], max.month = m, H = m / 12, beta = r$MKT, se = r$SE.MKT, p.beta1 = pfb(r$MKT, r$SE.MKT, r$hac.b),
      n.roots = r$n.roots, roots = paste(sprintf("%.3f", rt), collapse = "; "), at.bound = r$at.bound, active = r$active, V = r$hac.V, S = r$hac.S, b = r$hac.b) })) }))
write.csv(iv, file.path(out.dir, "iv_estimates_all_horizons.csv"), row.names = FALSE)

# ---- check 1: IV against an independent computation from the pipeline's unit errors ----
chk1 <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) { k <- paste(cells$type[i], cells$w[i]); prep <- prepare.units(U[[k]]$vyp$rows, fc)
  do.call(rbind, lapply(c(0, 60, 180), function(m) { r <- iv[iv$type == cells$type[i] & iv$weighting == cells$w[i] & iv$max.month == m, ]
    pe <- function(b, mm) vapply(prep$units, unit.horizon.error, 0, max.month = mm, par0 = c(RF = 1, MKT = b), sdf.model = "linear")
    zz <- pe(1, m) / pe(1, 0); gb <- function(b) mean(zz * pe(b, 0))
    yrs <- vapply(prep$units, `[[`, 0, "vintage"); V <- length(yrs); S <- ceiling(1.3 * sqrt(V)); s <- zz * pe(r$beta, 0)
    om <- 0; for (a in seq_len(V)) for (j in seq_len(V)) om <- om + max(0, 1 - abs(yrs[a] - yrs[j]) / S) * s[a] * s[j]
    Gd <- (gb(r$beta + 1e-4) - gb(r$beta - 1e-4)) / 2e-4
    data.frame(cell = k, H = m / 12, moment.at.estimate = gb(r$beta), se.ratio = r$se / sqrt(om / V / (V * Gd^2))) })) }))
cat("check 1 (IV against the pipeline's unit errors, explicit kernel sum): max |moment| at the estimate", signif(max(abs(chk1$moment.at.estimate)), 2),
    "| max |SE ratio - 1|", signif(max(abs(chk1$se.ratio - 1), na.rm = TRUE), 2), "\n")
stopifnot(max(abs(chk1$moment.at.estimate)) < 1e-8, max(abs(chk1$se.ratio - 1), na.rm = TRUE) < 1e-6)

# ---- L2 and U from the table data, refitted as a check; jackknife SEs of the gaps ----
tab.l2 <- rbind(read.csv(file.path(out.dir, "empirical_BO_US_MKT_data.csv")), read.csv(file.path(out.dir, "empirical_VC_US_MKT_data.csv")))
tab.l2 <- tab.l2[tab.l2$Factor == "MKT" & tab.l2$method == "Asymptotic", ]
tab.u <- rbind(read.csv(file.path(out.dir, "upair_BO_US_data.csv")), read.csv(file.path(out.dir, "upair_VC_US_data.csv")))
tab.u <- tab.u[tab.u$model == "MKT" & tab.u$parameter == "MKT" & tab.u$vintage.weights == "equal", ]
# point estimate only (the Brent fit of mc.fit(), without its sandwich; IV as mc.fit.iv)
fit.b <- function(est, fp, m) {
  if (est$criterion == "IV") return(mc.fit.iv(est, fp, "MKT", m)$MKT)
  f <- mc.criterion(est, fp, m); br <- feasible.interval(mc.domain(fp, "MKT"))
  stats::optim(c(MKT = 1), function(p) f(setNames(p, "MKT")), method = "Brent", lower = br[1], upper = br[2])$par[[1]]
}
res <- list(); repro <- 0
for (i in seq_len(nrow(cells))) { ty <- cells$type[i]; w <- cells$w[i]; k <- paste(ty, w); fv <- U[[k]]$vyp$fp; ff <- U[[k]]$fund$fp
  for (m in H.TAB) {
    l2 <- tab.l2[tab.l2$Type == ty & tab.l2$weighting == w & tab.l2$max.month == m, ]; u <- tab.u[tab.u$Type == ty & tab.u$weighting == w & tab.u$max.month == m, ]
    v <- iv[iv$type == ty & iv$weighting == w & iv$max.month == m, ]
    stopifnot(nrow(l2) == 1, nrow(u) == 1, nrow(v) == 1)
    b.l2 <- fit.b(est.L2, fv, m); b.u <- fit.b(est.U(w), ff, m)
    repro <- max(repro, abs(b.l2 - l2$MKT), abs(b.u - u$estimate))
    vs <- sort(unique(fv$vintage)) # leave-one-vintage-out (vintages with funds in the fund file)
    jk <- t(vapply(vs, function(x) c(fit.b(est.L2, fast.subset(fv, fv$vintage != x), m), fit.b(est.U(w), fast.subset(ff, ff$vintage != x), m),
      fit.b(MC.ESTIMATORS$IV, fast.subset(fv, fv$vintage != x), m)), numeric(3)))
    jse <- function(x) sqrt((length(x) - 1) / length(x) * sum((x - mean(x))^2))
    res[[length(res) + 1]] <- data.frame(type = ty, weighting = w, max.month = m, H = m / 12,
      l2 = l2$MKT, l2.se = if (isTRUE(l2$at.bound)) NA else l2$SE.MKT, l2.bound = isTRUE(l2$at.bound),
      u = u$estimate, u.se = if (isTRUE(u$at.bound)) NA else u$SE.A, u.bound = isTRUE(u$at.bound),
      iv = v$beta, iv.se = v$se, iv.bound = v$at.bound, iv.roots = v$n.roots, iv.root.list = v$roots,
      gap.l2u = l2$MKT - u$estimate, gap.l2u.jse = jse(jk[, 1] - jk[, 2]), gap.uiv = u$estimate - v$beta, gap.uiv.jse = jse(jk[, 2] - jk[, 3]),
      jk.V = length(vs)) } }
res <- do.call(rbind, res)
cat("check 2 (L2 and U refitted with the fast criterion = the G1/G2 table data): max |difference|", signif(repro, 2), "\n")
stopifnot(repro < 1e-5)

# ---- predicted L2 - U gap from the data-shaped Monte Carlo (paired over replications) ----
pred <- function(type, w, m) {
  d <- if (w == "EW") paste0("core_", type) else paste0("s_core_", type); r <- mc.read(d); r <- r[r$model == "MKT" & r$max.month == m, ]
  e1 <- if (w == "EW") "L2.vyp" else "L2.vyp.fw"; e2 <- if (w == "EW") "U.eq" else "U.eq.fw"
  a <- r[r$estimator == e1, c("rep", "MKT")]; b <- r[r$estimator == e2, c("rep", "MKT")]
  if (!nrow(a) || !nrow(b)) return(c(NA, NA)); x <- merge(a, b, by = "rep"); g <- x$MKT.x - x$MKT.y
  c(mean(g), sd(g) / sqrt(length(g))) }
pg <- t(vapply(seq_len(nrow(res)), function(j) pred(res$type[j], res$weighting[j], res$max.month[j]), numeric(2)))
res$pred.l2u <- pg[, 1]; res$pred.l2u.mcse <- pg[, 2]
write.csv(res, file.path(out.dir, "iv_comparison_rows.csv"), row.names = FALSE)

# ---- table ---------------------------------------------------------------------------
cells.rec <- list()
rec <- function(key, column, value, printed) { cells.rec[[length(cells.rec) + 1]] <<- data.frame(key, column = column, value = value, printed = printed, stringsAsFactors = FALSE); printed }
f3 <- function(x) if (is.na(x)) "--" else formatC(round(x, 3) + 0, format = "f", digits = 3)
es <- function(key, col, b, se, bound) sprintf("%s%s (%s)", rec(key, col, b, f3(b)), if (bound) "$^\\dagger$" else "", rec(key, paste0(col, ".se"), if (bound) NA else se, f3(if (bound) NA else se)))
lines <- character(0)
for (ty in c("BO", "VC")) for (w in c("EW", "FW")) {
  lines <- c(lines, sprintf("\\multicolumn{8}{l}{\\textit{%s, %s}} \\\\", c(BO = "Buyout", VC = "Venture capital")[[ty]], c(EW = "equally weighted", FW = "fund-size weighted")[[w]]))
  for (m in H.TAB) { r <- res[res$type == ty & res$weighting == w & res$max.month == m, ]; key <- data.frame(type = ty, weighting = w, H = m / 12)
    pr <- if (is.finite(r$pred.l2u)) sprintf("%s (%s)", rec(key, "pred.l2u", r$pred.l2u, f3(r$pred.l2u)), rec(key, "pred.l2u.mcse", r$pred.l2u.mcse, f3(r$pred.l2u.mcse))) else "--"
    lines <- c(lines, sprintf("\\quad %s yr & %s & %s & %s & %s & %s (%s) & %s & %s (%s) \\\\", format(m / 12),
      es(key, "l2", r$l2, r$l2.se, r$l2.bound), es(key, "u", r$u, r$u.se, r$u.bound), es(key, "iv", r$iv, r$iv.se, r$iv.bound),
      rec(key, "iv.roots", r$iv.roots, format(r$iv.roots)),
      rec(key, "gap.l2u", r$gap.l2u, f3(r$gap.l2u)), rec(key, "gap.l2u.jse", r$gap.l2u.jse, f3(r$gap.l2u.jse)), pr,
      rec(key, "gap.uiv", r$gap.uiv, f3(r$gap.uiv)), rec(key, "gap.uiv.jse", r$gap.uiv.jse, f3(r$gap.uiv.jse)))) } }
tex <- c(sprintf("%% Auto-generated by r_project/results/empirical_iv.R at commit %s%s", system("git rev-parse --short HEAD", intern = TRUE),
    if (length(system("git status --porcelain --untracked-files=no", intern = TRUE))) " (uncommitted changes)" else ""),
  "\\begin{table}[htbp]", "\\centering", "\\footnotesize",
  "\\caption{Single-factor $\\hat\\beta_{MKT}$ by L2, $U_{v,H}$ and the IV moment, and the gap between L2 and $U_{v,H}$ against the noise tilt}", "\\label{tab:iv_comparison}",
  "\\resizebox{\\linewidth}{!}{%", "\\begin{tabular}{l rrr r rrr}", "\\toprule",
  "Horizon & L2 & $U_{v,H}$, $q_v=1$ & IV & Roots & L2$-U_{v,H}$ & predicted & $U_{v,H}-$IV \\\\", "\\midrule",
  lines, "\\bottomrule", "\\end{tabular}%", "}", "\\par\\vspace{2pt}",
  paste0("{\\scriptsize\\textit{Notes:} North America, vintages up to 2010, single-factor $\\hat\\beta_{MKT}$ (SE). L2 on vintage portfolios and $U_{v,H}$ with $q_v=1$: the estimates of the main tables. ",
    "IV: the smallest root in $\\beta$ of $V^{-1}\\sum_vz_vPV_v(\\beta)=0$ with the benchmark-path instrument $z_v=A_{v,H}(1)/A_{v,0}(1)$ ($z_v=1$ at $H=0$: the mean pricing moment), on the admissible domain; Roots: the number of roots on the domain (all are listed in the data file). ",
    "SEs: Bartlett sandwich, $S=\\lceil1.3\\sqrt V\\rceil$. L2$-U_{v,H}$ and $U_{v,H}-$IV: differences of the estimates, with leave-one-vintage-out jackknife SEs in brackets (descriptive, no test). ",
    "Predicted: L2 bias minus $U_{v,H}$ bias in the data-shaped Monte Carlo without vintage-level effects (the noise tilt; EW: the core designs, FW: the core designs with the data's fund-size dispersion, simulated at 0 and 15 years only), Monte Carlo SE in brackets. ",
    "In the Monte Carlo, under market-correlated vintage shocks IV at $H=0$ lies below L2 and $U_{v,H}$, and IV at $H>0$ is not robust. $^\\dagger$: estimate on the boundary of the domain, no SE.}"),
  "\\end{table}")
writeLines(tex, file.path(out.dir, "iv_comparison.tex"))
cr <- do.call(rbind, cells.rec)
write.csv(cr, file.path(out.dir, "iv_comparison_data.csv"), row.names = FALSE)
cat("wrote", file.path(out.dir, "iv_comparison.tex"), "with", nrow(cr), "recorded numbers\n")

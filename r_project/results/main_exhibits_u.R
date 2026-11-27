# U-first main exhibits of the empirical section (revision plan, Section L;
# reconciled Nov 22 fixes).
# Run from r_project/:  LANG=en_US.UTF-8 Rscript results/main_exhibits_u.R
#
# No estimation: reads existing result files only.
#   L2 on vintage portfolios: results/figures_bovc_us_g1/empirical_{BO,VC}_US_MKT_data.csv
#     (G1 table data, r5; rows method == "Asymptotic"; Factor "MKT" = single
#     factor, "Alpha" = (alpha, beta) model with Coef = alpha).
#   U_{v,H}, q_v = 1: results/figures_bovc_us_g1/upair_{BO,VC}_US_data.csv (G2
#     table data; vintage.weights == "equal").
#   Expected core L2 - U gap and matched-vintage L2:
#     results/figures_bovc_us_g1/iv_comparison_rows.csv
#   VC (alpha, beta) without 1993-1997:
#     results/figures_bovc_us_g1/vc_alpha_without_1993_1997.csv
# Writes to results/figures_bovc_us_g1/main_u/:
#   figure_E1_beta_profiles.pdf/.tex  main text: single-factor beta across horizons (fixed-b 95% bands)
#   table_E1_beta.tex                 main text: single-factor beta with p(beta=1), H = 0, 5, 10, 15, 30
#   table_E1_beta_all.tex             Online Appendix: the same at all 12 horizons
#   table_E2_alpha_beta.tex           main text: (alpha, beta), estimates & SEs (no stars per calibration)
#   table_OA_vc_without_1993_1997.tex Online Appendix: VC (alpha, beta) without 1993-1997
#   exhibits_data.csv                 every printed number with its source
suppressPackageStartupMessages(library(ggplot2))
source("helper/fixedb.R")
src.dir <- "results/figures_bovc_us_g1"
out.dir <- file.path(src.dir, "main_u")
dir.create(out.dir, showWarnings = FALSE)
H.ALL <- c(0, 30, 60, 120, 150, 180, 210, 240, 270, 300, 330, 360)
H.MAIN <- c(0, 60, 120, 180, 360)
H.E2 <- c(0, 30, 60, 120, 180, 240, 360)
TYPES <- c(BO = "Buyout", VC = "Venture capital")
WTS <- c(EW = "equally weighted", FW = "fund-size weighted")

# ---- inputs -------------------------------------------------------------------------
rd <- function(f) read.csv(file.path(src.dir, f), stringsAsFactors = FALSE)
l2 <- rbind(rd("empirical_BO_US_MKT_data.csv"), rd("empirical_VC_US_MKT_data.csv"))
l2 <- l2[l2$method == "Asymptotic" & l2$Factor %in% c("MKT", "Alpha"), ]
u <- rbind(rd("upair_BO_US_data.csv"), rd("upair_VC_US_data.csv"))
u <- u[u$vintage.weights == "equal" & u$model %in% c("MKT", "Alpha"), ]
pred <- rd("iv_comparison_rows.csv")
vcx <- rd("vc_alpha_without_1993_1997.csv")
stopifnot(all(l2$pipeline.version == "2026-11-17-r5"), !anyNA(match(H.ALL, l2$max.month)), !anyNA(match(H.ALL, u$max.month)))
if (any(l2$git.dirty) || any(u$git.dirty)) warning("Estimation rows were generated from a dirty git tree (recorded in headers)")

# one cell: estimate, SE, fixed-b p (test of zero), fixed-b p (test of 1 for MKT), b, bound flag
cell.l2 <- function(type, w, m, model, par) {
  r <- l2[l2$Type == type & l2$weighting == w & l2$max.month == m & l2$Factor == model, ]
  if (nrow(r) != 1) stop("L2: not one row for ", type, w, m, model)
  hb <- r$hac.b
  if (par == "MKT") {
    p1 <- if (!isTRUE(r$at.bound) && is.finite(r$SE.MKT) && is.finite(hb)) fixedb.q1.pvalue(((r$MKT - 1) / r$SE.MKT)^2, hb) else NA_real_
    list(est = r$MKT, se = r$SE.MKT, p = r$p.MKT, p1 = p1, b = hb, bound = isTRUE(r$at.bound), src = "empirical_*_US_MKT_data.csv")
  } else {
    list(est = r$Coef, se = r$SE.Coef, p = r$p.Coef, p1 = NA_real_, b = hb, bound = isTRUE(r$at.bound), src = "empirical_*_US_MKT_data.csv")
  }
}
cell.u <- function(type, w, m, model, par) {
  r <- u[u$Type == type & u$weighting == w & u$max.month == m & u$model == model & u$parameter == par, ]
  if (nrow(r) != 1) stop("U: not one row for ", type, w, m, model, par)
  hb <- r$hac.b
  p1 <- if (par == "MKT" && !isTRUE(r$at.bound) && is.finite(r$SE.A) && is.finite(hb)) fixedb.q1.pvalue(((r$estimate - 1) / r$SE.A)^2, hb) else NA_real_
  list(est = r$estimate, se = r$SE.A, p = r$p, p1 = p1, b = hb, bound = isTRUE(r$at.bound), src = "upair_*_US_data.csv")
}

# ---- recording and formatting --------------------------------------------------------
cells <- list()
rec <- function(exhibit, key, column, value, printed, source) {
  cells[[length(cells) + 1]] <<- data.frame(exhibit = exhibit, key, column = column, value = value, printed = printed, source = source, stringsAsFactors = FALSE)
  printed
}
f3 <- function(x) if (is.na(x)) "--" else formatC(round(x, 3) + 0, format = "f", digits = 3)
fp <- function(p) if (is.na(p)) "--" else if (p < 0.001) "<0.001" else formatC(round(p, 3) + 0, format = "f", digits = 3)
stars <- function(p) if (is.na(p)) "" else if (p < 0.001) "$^{***}$" else if (p < 0.01) "$^{**}$" else if (p < 0.05) "$^{*}$" else ""

est.cell <- function(ex, key, c, col, scale = 1, marks = TRUE) {
  e <- rec(ex, key, col, c$est * scale, f3(c$est * scale), c$src)
  mk <- if (c$bound) "$^{\\dagger}$" else if (marks) stars(c$p) else ""
  if (marks && !c$bound) rec(ex, key, paste0(col, ".p"), c$p, gsub("[^*]", "", mk), c$src)
  s <- if (c$bound) rec(ex, key, paste0(col, ".se"), NA, "--", c$src) else rec(ex, key, paste0(col, ".se"), c$se * scale, f3(c$se * scale), c$src)
  c(paste0(e, mk), s)
}

GIT.STATE <- sprintf("%s%s", system("git rev-parse --short HEAD", intern = TRUE),
  if (length(system("git status --porcelain --untracked-files=no", intern = TRUE))) " (uncommitted changes)" else "")
header <- function() c(sprintf("%% Auto-generated by r_project/results/main_exhibits_u.R at commit %s", GIT.STATE),
  sprintf("%% Estimation results: L2 commit(s) %s%s; U_{v,H} commit(s) %s; no estimation in this script",
    paste(unique(l2$git.commit), collapse = ", "), if (isTRUE(any(l2$git.dirty))) " (run with uncommitted changes)" else "",
    paste(unique(u$git.commit), collapse = ", ")), "")
write.tex <- function(file, body) { writeLines(c(header(), body), file.path(out.dir, file)); message("wrote ", file.path(out.dir, file)) }
hlab <- function(m) paste(format(m / 12), "yr")

nV <- function(d, ty) { v <- unique(d$hac.V[d$Type == ty]); stopifnot(length(v) == 1); v }
bU.BO <- unique(u$hac.b[u$Type == "BO" & u$model == "MKT"]); stopifnot(length(bU.BO) == 1)
bU.VC <- unique(u$hac.b[u$Type == "VC" & u$model == "MKT"]); stopifnot(length(bU.VC) == 1)
crit.BO <- fixedb.crit.t(0.05, bU.BO)
crit.VC <- fixedb.crit.t(0.05, bU.VC)

NOTE.COMMON <- paste0(sprintf("North America, vintages up to 2010 (FW baseline: within-vintage median fund-size imputation; sensitivity in Table~\\ref{tab:fw_imputation_sensitivity}): Buyout %d vintages; VC %d, of which the %d with at least two funds enter $U_{v,H}$. ",
    nV(l2, "BO"), nV(l2, "VC"), nV(u, "VC")), if (nV(u, "BO") != nV(l2, "BO")) stop("Buyout vintage counts differ between L2 and U") else "",
  "$U_{v,H}$: the cross-fund pair-product criterion with equal vintage weights ($q_v=1$); under Theorem~T7 (conditional independence of funds given the factor path, eligibility $n_v\\ge2$, fund sizes, and vintage origin, zero sign-aligned conditional mean pricing error, no shared vintage effects $c_v=0$, and unique identification), it targets $\\theta_0$ at every horizon, $H=0$ included. ",
  "L2: per-commitment squared loss on vintage portfolios (matching the inverse-FoF-size normalization of \\citealp[p.~533]{DLP12} at $H=0$; unnormalized dollar sums scale each vintage squared error by $S_v^2$), whose target is tilted by within-vintage fund noise and moves with $H$ (matched-vintage 23-vintage VC L2 reported in Table~\\ref{tab:iv_comparison}). ",
  "SE: Bartlett sandwich over vintage years, $S=\\lceil1.3\\sqrt V\\rceil$. ")
NOTE.STARS.E1 <- "Markers on $U_{v,H}$: two-sided fixed-$b$ $p$-values of a test of $\\beta=0$, $^{*}$\\,$p<0.05$, $^{**}$\\,$p<0.01$, $^{***}$\\,$p<0.001$; column $p_1$ gives the two-sided fixed-$b$ $p$-value of $\\beta=1$ (under shared shocks T8, per-vintage score centering can fail at $H>0$, Remark~\\ref{rem:per_vintage_centering}). FW single-factor markers are qualified by mild finite-sample over-rejection (8.1\\% Buyout / 8.5\\% VC at $H=0$ in Design~B, Table~\\ref{tab:mc_hardening}). $^{\\dagger}$: estimate on a bound of the admissible domain ($\\beta\\in[-1,4]$, $\\alpha\\in[-2\\%,2\\%]$ per month, gross benchmark return at least 0.08), no SE. "

# ---- Table E1: single-factor beta (with p(beta=1) on U_{v,H}) --------------------------
table.e1 <- function(H, file, label, caption, ex) {
  lines <- character(0)
  for (ty in names(TYPES)) {
    lines <- c(lines, sprintf("\\multicolumn{11}{l}{\\textit{%s}} \\\\", TYPES[[ty]]))
    for (m in H) { key <- data.frame(type = ty, H = m / 12)
      v <- character(0)
      for (w in names(WTS)) { k <- cbind(key, weighting = w)
        cu <- cell.u(ty, w, m, "MKT", "MKT"); cl <- cell.l2(ty, w, m, "MKT", "MKT")
        p1.str <- if (cu$bound) rec(ex, k, "U.p1", NA_real_, "--", cu$src) else rec(ex, k, "U.p1", cu$p1, fp(cu$p1), cu$src)
        v <- c(v, est.cell(ex, k, cu, "U.beta"), p1.str, est.cell(ex, k, cl, "L2.beta", marks = FALSE)) }
      lines <- c(lines, sprintf("\\quad %s & %s \\\\", hlab(m), paste(v, collapse = " & "))) }
    if (ty == "BO") lines <- c(lines, "\\addlinespace")
  }
  write.tex(file, c("\\begin{table}[htbp]", "\\centering", "\\footnotesize", sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{l rrr rr rrr rr}", "\\toprule",
    " & \\multicolumn{5}{c}{Equally weighted funds (EW)} & \\multicolumn{5}{c}{Fund-size weighted funds (FW)} \\\\",
    "\\cmidrule(lr){2-6} \\cmidrule(lr){7-11}",
    " & \\multicolumn{3}{c}{$U_{v,H}$ (main)} & \\multicolumn{2}{c}{L2} & \\multicolumn{3}{c}{$U_{v,H}$ (main)} & \\multicolumn{2}{c}{L2} \\\\",
    "\\cmidrule(lr){2-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-9} \\cmidrule(lr){10-11}",
    "Horizon & $\\hat\\beta$ & SE & $p_1$ & $\\hat\\beta$ & SE & $\\hat\\beta$ & SE & $p_1$ & $\\hat\\beta$ & SE \\\\", "\\midrule",
    lines, "\\bottomrule", "\\end{tabular}%", "}", "\\par\\vspace{2pt}",
    paste0("{\\scriptsize\\textit{Notes:} Single-factor $\\hat\\beta_{MKT}$ ($\\alpha$ fixed at zero). ", NOTE.COMMON, NOTE.STARS.E1,
      "The headline estimate is $U_{v,H}$ at $H=0$: in our data-shaped core Monte Carlo, $H=0$ minimizes $U_{v,H}$ RMSE for single-factor $\\beta$ and $(\\alpha,\\beta)$ (whereas in two-factor models, Table~\\ref{tab:mc_two_factor}, $H=15$ yr lowers $U_{v,H}$'s second-loading RMSE for ME, IA, and ROE while market-beta RMSE rises in ROE).",
      if (length(H) < length(H.ALL)) " All horizons: Table~\\ref{tab:u_l2_beta_all}." else "", "}"),
    "\\end{table}"))
}
table.e1(H.MAIN, "table_E1_beta.tex", "tab:u_l2_beta", "Single-factor $\\hat\\beta_{MKT}$ ($\\alpha$ fixed at zero): $U_{v,H}$ and L2", "E1")
table.e1(H.ALL, "table_E1_beta_all.tex", "tab:u_l2_beta_all", "Single-factor $\\hat\\beta_{MKT}$ ($\\alpha$ fixed at zero): $U_{v,H}$ and L2, all horizons", "E1.all")

# ---- Table E1 (JBES main): single-factor beta + compact Panel B inception diagnostics ---
shac.cmp <- if (file.exists(file.path(src.dir, "calendar_vs_index_shac.csv"))) rd("calendar_vs_index_shac.csv") else NULL
panel.b.lines <- character(0)
for (ty in names(TYPES)) for (w in names(WTS)) {
  r0 <- pred[pred$type == ty & pred$weighting == w & pred$max.month == 0, ]
  stopifnot(nrow(r0) == 1)
  gap.s <- if (r0$V.full == r0$V.matched) sprintf("%s (%s)", f3(r0$gap.l2u), f3(r0$gap.l2u.jse)) else sprintf("%s [%s]", f3(r0$gap.l2u), f3(r0$gap.l2u.matched))
  pred.s <- sprintf("%s (%s)", f3(r0$pred.l2u), f3(r0$pred.l2u.mcse))
  neff.s <- if (!is.null(shac.cmp)) {
    nx <- shac.cmp$neff[shac.cmp$type == ty & shac.cmp$weighting == w & shac.cmp$estimator == "U.eq" & shac.cmp$H == 0]
    sprintf("%.1f", nx[1])
  } else "--"
  panel.b.lines <- c(panel.b.lines, sprintf("\\quad %s, %s & %s & %s & %s & %d/%d & %s & %s \\\\",
    TYPES[[ty]], w, gap.s, pred.s, r0$iv.root.list, r0$iv.loo.no.root, r0$V.full, neff.s, r0$ar.set))
}
e1.lines <- character(0)
for (ty in names(TYPES)) {
  e1.lines <- c(e1.lines, sprintf("\\multicolumn{11}{l}{\\textit{%s}} \\\\", TYPES[[ty]]))
  for (m in H.MAIN) { key <- data.frame(type = ty, H = m / 12)
    v <- character(0)
    for (w in names(WTS)) { k <- cbind(key, weighting = w)
      cu <- cell.u(ty, w, m, "MKT", "MKT"); cl <- cell.l2(ty, w, m, "MKT", "MKT")
      p1.str <- if (cu$bound) "--" else fp(cu$p1)
      u.c <- c(paste0(f3(cu$est), if (cu$bound) "$^{\\dagger}$" else stars(cu$p)), if (cu$bound) "--" else f3(cu$se))
      l.c <- c(paste0(f3(cl$est), if (cl$bound) "$^{\\dagger}$" else ""), if (cl$bound) "--" else f3(cl$se))
      v <- c(v, u.c, p1.str, l.c) }
    e1.lines <- c(e1.lines, sprintf("\\quad %s & %s \\\\", hlab(m), paste(v, collapse = " & "))) }
  if (ty == "BO") e1.lines <- c(e1.lines, "\\addlinespace")
}
write.tex("table_E1_beta_jbes.tex", c("\\begin{table}[htbp]", "\\centering", "\\footnotesize",
  "\\caption{Single-factor $\\hat\\beta_{\\mathrm{MKT}}$ ($\\alpha$ fixed at zero): $U_{v,H}$, L2, and inception ($H=0$) identification diagnostics}",
  "\\label{tab:u_l2_beta}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{l rrr rr rrr rr}", "\\toprule",
  "\\multicolumn{11}{l}{\\textit{Panel A: Single-factor $\\hat\\beta_{\\mathrm{MKT}}$ ($\\alpha$ fixed at zero) across horizons}} \\\\",
  "\\addlinespace[2pt]",
  " & \\multicolumn{5}{c}{Equally weighted funds (EW)} & \\multicolumn{5}{c}{Fund-size weighted funds (FW)} \\\\",
  "\\cmidrule(lr){2-6} \\cmidrule(lr){7-11}",
  " & \\multicolumn{3}{c}{$U_{v,H}$ (main)} & \\multicolumn{2}{c}{L2} & \\multicolumn{3}{c}{$U_{v,H}$ (main)} & \\multicolumn{2}{c}{L2} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-9} \\cmidrule(lr){10-11}",
  "Horizon & $\\hat\\beta$ & SE & $p_1$ & $\\hat\\beta$ & SE & $\\hat\\beta$ & SE & $p_1$ & $\\hat\\beta$ & SE \\\\", "\\midrule",
  e1.lines, "\\bottomrule", "\\end{tabular}%", "}",
  "\\par\\vspace{4pt}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{l r r r c r c}", "\\toprule",
  "\\multicolumn{7}{l}{\\textit{Panel B: Inception ($H=0$, $\\alpha=0$) estimator gaps, IV root fragility, $N_{\\mathrm{eff}}$, and nominal 95\\% AR sets}} \\\\",
  "\\addlinespace[2pt]",
  "Series & L2$-U_{v,0}$ & Expected core gap & IV roots & LOO no-root & $N_{\\mathrm{eff}}(U_{v,0})$ & Nominal 95\\% AR set \\\\", "\\midrule",
  panel.b.lines, "\\bottomrule", "\\end{tabular}%", "}",
  "\\par\\vspace{2pt}",
  paste0("{\\scriptsize\\textit{Notes:} Single-factor $\\hat\\beta_{\\mathrm{MKT}}$ ($\\alpha$ fixed at zero). ", NOTE.COMMON, NOTE.STARS.E1,
    "Panel~B reports inception ($H=0$) L2$-U_{v,0}$ gaps (Buyout: LOO jackknife SE in parentheses; VC: full-sample 27-vintage L2 minus 23-vintage $U_{v,0}$, followed in brackets by the matched 23-vintage difference), paired expected core MC gaps (MC SE), all roots of the inception mean pricing moment $V^{-1}\\sum_v PV_v(\\beta)=0$ on $\\Theta$ (using all 25 Buyout / 27 VC vintages), delete-one-vintage refits losing all roots, $U_{v,0}$ score-concentration count $N_{\\mathrm{eff}}$, and centered Bartlett fixed-$b$ nominal 95\\% AR confidence sets \\eqref{eq:ar_set} (all 27 VC vintages). Full 12-horizon estimates and four-horizon IV/AR diagnostics: Online Supplement Tables~\\ref{tab:u_l2_beta_all} and~\\ref{tab:iv_comparison}.}"),
  "\\end{table}"))

# ---- Table E2: (alpha, beta) model (stars suppressed per calibration audit) -------------
make.e2 <- function(H.vec, file, label, caption, record = FALSE) {
  lines <- character(0)
  for (ty in names(TYPES)) for (w in names(WTS)) {
    lines <- c(lines, sprintf("\\multicolumn{5}{l}{\\textit{%s, %s}} \\\\", TYPES[[ty]], WTS[[w]]))
    for (m in H.vec) { k <- data.frame(type = ty, weighting = w, H = m / 12)
      if (record) {
        cu.b <- est.cell("E2", k, cell.u(ty, w, m, "Alpha", "MKT"), "U.beta", marks = FALSE)
        cu.a <- est.cell("E2", k, cell.u(ty, w, m, "Alpha", "Alpha"), "U.alpha", scale = 100, marks = FALSE)
        cl.b <- est.cell("E2", k, cell.l2(ty, w, m, "Alpha", "MKT"), "L2.beta", marks = FALSE)
        cl.a <- est.cell("E2", k, cell.l2(ty, w, m, "Alpha", "Alpha"), "L2.alpha", scale = 100, marks = FALSE)
      } else {
        fmt.c <- function(c, scale = 1) c(paste0(f3(c$est * scale), if (c$bound) "$^{\\dagger}$" else ""), if (c$bound) "--" else f3(c$se * scale))
        cu.b <- fmt.c(cell.u(ty, w, m, "Alpha", "MKT"))
        cu.a <- fmt.c(cell.u(ty, w, m, "Alpha", "Alpha"), 100)
        cl.b <- fmt.c(cell.l2(ty, w, m, "Alpha", "MKT"))
        cl.a <- fmt.c(cell.l2(ty, w, m, "Alpha", "Alpha"), 100)
      }
      v <- vapply(list(cu.b, cu.a, cl.b, cl.a), function(x) sprintf("%s (%s)", x[1], x[2]), "")
      lines <- c(lines, sprintf("\\quad %s & %s \\\\", hlab(m), paste(v, collapse = " & "))) }
    if (!(ty == "VC" && w == "FW")) lines <- c(lines, "\\addlinespace")
  }
  write.tex(file, c("\\begin{table}[htbp]", "\\centering", "\\footnotesize",
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
    "\\begin{tabular}{l rr rr}", "\\toprule",
    " & \\multicolumn{2}{c}{$U_{v,H}$ (main)} & \\multicolumn{2}{c}{L2} \\\\", "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
    "Horizon & $\\hat\\beta$ (SE) & $\\hat\\alpha$ (SE) & $\\hat\\beta$ (SE) & $\\hat\\alpha$ (SE) \\\\", "\\midrule",
    lines, "\\bottomrule", "\\end{tabular}", "\\par\\vspace{2pt}",
    paste0("{\\scriptsize\\textit{Notes:} $\\hat\\alpha$ in percent per month. ", NOTE.COMMON,
      "Significance stars are omitted in this table: in the EW core Monte Carlo at $H=0$, the nominal 5\\% fixed-$b$ test of $\\alpha=0$ for $U_{v,H}$ rejects 12.1\\% (Buyout) and 10.9\\% (VC) of draws (Table~\\ref{tab:mc_alpha_beta}); FW joint calibration is not separately simulated; and at 30 years conditional rejection collapses to 0.4--0.5\\% for $\\alpha$ and 0.3--1.1\\% for $\\beta$ with 6.4--8.0\\% on the boundary. ",
      "SEs are therefore reported as approximate precision measures only. $^{\\dagger}$: estimate on a bound of the admissible domain, no SE. ",
      "In VC the estimate at $H=0$ sits on the $+2\\%$ monthly $\\alpha$ bound under both criteria because the criterion keeps falling in $\\alpha$; that comes from the 1993--1997 vintages as a group, and without them $\\hat\\alpha$ is negative at every horizon (Online Supplement Table~\\ref{tab:vc_without_1993_1997}). ",
      if (length(H.vec) < length(H.E2)) "Complete 7-horizon estimates: Online Supplement Table~\\ref{tab:u_l2_alpha_beta_all}.}" else "The horizon then imposes a minimum in $\\alpha$ without identifying it.}"),
    "\\end{table}"))
}
make.e2(H.E2, "table_E2_alpha_beta.tex", "tab:u_l2_alpha_beta_legacy", "The $(\\alpha,\\beta)$ model: $U_{v,H}$ and L2", record = TRUE)
make.e2(H.E2, "table_E2_alpha_beta_all.tex", "tab:u_l2_alpha_beta_all", "The $(\\alpha,\\beta)$ model: $U_{v,H}$ and L2 (all seven horizons)", record = FALSE)
make.e2(c(0, 60, 180, 360), "table_E2_alpha_beta_jbes.tex", "tab:u_l2_alpha_beta", "The $(\\alpha,\\beta)$ model: $U_{v,H}$ and L2", record = FALSE)

# ---- Online Appendix: VC (alpha, beta) without 1993-1997 -----------------------------------
chk.full <- 0
for (i in which(vcx$sample == "full")) { r <- vcx[i, ]
  ref <- if (r$estimator == "L2") cell.l2("VC", r$weighting, r$max.month, "Alpha", "MKT") else cell.u("VC", r$weighting, r$max.month, "Alpha", "MKT")
  refa <- if (r$estimator == "L2") cell.l2("VC", r$weighting, r$max.month, "Alpha", "Alpha") else cell.u("VC", r$weighting, r$max.month, "Alpha", "Alpha")
  chk.full <- max(chk.full, abs(r$beta - ref$est), abs(r$alpha - refa$est), as.numeric(r$at.bound != ref$bound)) }
stopifnot(chk.full < 1e-8, nrow(vcx) == 2 * 2 * 2 * length(H.E2))
message("VC diagnostic: full-sample fits = G1/G2 table data (largest difference ", signif(chk.full, 2), ")")
lines <- character(0)
for (w in names(WTS)) for (e in c("U", "L2")) {
  lines <- c(lines, sprintf("\\multicolumn{5}{l}{\\textit{%s, %s}} \\\\", if (e == "U") "$U_{v,H}$ ($q_v=1$)" else "L2", WTS[[w]]))
  for (m in H.E2) { k <- data.frame(type = "VC", weighting = w, estimator = e, H = m / 12)
    v <- character(0)
    for (smp in c("full", "without 1993-1997")) { r <- vcx[vcx$weighting == w & vcx$estimator == e & vcx$max.month == m & vcx$sample == smp, ]; stopifnot(nrow(r) == 1)
      kk <- cbind(k, sample = smp); dg <- if (isTRUE(r$at.bound)) "$^{\\dagger}$" else ""
      v <- c(v, paste0(rec("OA.vc", kk, "beta", r$beta, f3(r$beta), "vc_alpha_without_1993_1997.csv"), dg),
                paste0(rec("OA.vc", kk, "alpha", 100 * r$alpha, f3(100 * r$alpha), "vc_alpha_without_1993_1997.csv"), dg)) }
    lines <- c(lines, sprintf("\\quad %s & %s \\\\", hlab(m), paste(v, collapse = " & "))) }
  if (!(w == "FW" && e == "L2")) lines <- c(lines, "\\addlinespace")
}
nv <- unique(vcx[, c("sample", "n.vintages")]); stopifnot(nrow(nv) == 2)
write.tex("table_OA_vc_without_1993_1997.tex", c("\\begin{table}[htbp]", "\\centering", "\\footnotesize",
  "\\caption{Venture capital, $(\\alpha,\\beta)$ model: full sample and without the 1993--1997 vintages}", "\\label{tab:vc_without_1993_1997}",
  "\\begin{tabular}{l rr rr}", "\\toprule",
  sprintf(" & \\multicolumn{2}{c}{Full sample (%d vintages)} & \\multicolumn{2}{c}{Without 1993--1997 (%d vintages)} \\\\",
    nv$n.vintages[nv$sample == "full"], nv$n.vintages[nv$sample == "without 1993-1997"]),
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}", "Horizon & $\\hat\\beta$ & $\\hat\\alpha$ & $\\hat\\beta$ & $\\hat\\alpha$ \\\\", "\\midrule",
  lines, "\\bottomrule", "\\end{tabular}", "\\par\\vspace{2pt}",
  paste0("{\\scriptsize\\textit{Notes:} Point estimates only, as a description: the block of vintages was chosen after looking at the data (it is the group of vintages whose terms push the criterion down in $\\alpha$ at the $H=0$ estimate), so no inference is reported. ",
    "$\\hat\\alpha$ in percent per month. Full sample: the estimates of Table~\\ref{tab:u_l2_alpha_beta}. Vintage counts in the header are those of L2; $U_{v,H}$ uses the vintages with at least two funds. ",
    "Same estimators, domain and data as there; refits on the vintages outside 1993--1997. $^{\\dagger}$: estimate on a bound of the admissible domain.}"),
  "\\end{table}"))

# ---- Figure E1: single-factor beta across horizons (pointwise fixed-b 95% bands) --------
COL <- c(teal = "#00897B", orange = "#C26A00", purple = "#AB23CC")
theme_paper <- theme_minimal(base_size = 10, base_family = "serif") +
  theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.3), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5), panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(face = "bold", size = 9), strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.5),
        axis.title = element_text(face = "bold", size = 9), axis.text = element_text(size = 8, color = "grey20"),
        axis.ticks = element_line(color = "grey40", linewidth = 0.3), legend.position = "bottom",
        legend.title = element_blank(), legend.text = element_text(size = 8))
SER <- c("U (q = 1), fixed-b 95% band", "L2, vintage portfolios", "U + expected core gap")
fd <- do.call(rbind, lapply(names(TYPES), function(ty) {
  crit.ty <- if (ty == "BO") crit.BO else crit.VC
  do.call(rbind, lapply(names(WTS), function(w) do.call(rbind, lapply(H.ALL, function(m) {
    cu <- cell.u(ty, w, m, "MKT", "MKT"); cl <- cell.l2(ty, w, m, "MKT", "MKT"); k <- data.frame(type = ty, weighting = w, H = m / 12)
    rec("F.E1", k, "U.beta", cu$est, NA, cu$src); rec("F.E1", k, "U.se", if (cu$bound) NA else cu$se, NA, cu$src); rec("F.E1", k, "L2.beta", cl$est, NA, cl$src)
    pr <- pred[pred$type == ty & pred$weighting == w & pred$max.month == m, ]
    pv <- if (nrow(pr) == 1 && is.finite(pr$pred.l2u)) pr$pred.l2u else NA
    if (is.finite(pv)) { rec("F.E1", k, "pred.l2u", pv, NA, "iv_comparison_rows.csv"); stopifnot(abs(pr$u - cu$est) < 1e-12) }
    data.frame(type = TYPES[[ty]], weighting = paste0(w, ": ", WTS[[w]]), H = m / 12, series = SER,
      value = c(cu$est, cl$est, cu$est + pv),
      lo = c(if (cu$bound) NA else cu$est - crit.ty * cu$se, NA, NA),
      hi = c(if (cu$bound) NA else cu$est + crit.ty * cu$se, NA, NA),
      bound = c(cu$bound, cl$bound, FALSE))
  }))))
}))
fd <- fd[is.finite(fd$value), ]
fd$series <- factor(fd$series, levels = SER); fd$type <- factor(fd$type, levels = TYPES)
E1.COL <- setNames(c(COL[["purple"]], COL[["teal"]], COL[["orange"]]), SER)
E1.LTY <- setNames(c("solid", "dashed", "blank"), SER)
E1.SHP <- setNames(c(15, 16, 2), SER)
lines.df <- fd[fd$series != SER[3], ]
p <- ggplot(fd, aes(H, value, colour = series, linetype = series, shape = series)) +
  geom_ribbon(data = fd[fd$series == SER[1] & is.finite(fd$lo), ], aes(ymin = lo, ymax = hi), fill = COL[["purple"]], alpha = 0.13, colour = NA, show.legend = FALSE) +
  geom_line(data = lines.df, linewidth = 0.7) +
  geom_point(data = fd[!fd$bound, ], size = 1.7, stroke = 0.8) +
  geom_point(data = fd[fd$bound, ], shape = 4, size = 2.4, stroke = 0.9, show.legend = FALSE) +
  facet_grid(type ~ weighting, scales = "free_y") +
  scale_colour_manual(values = E1.COL, breaks = SER, limits = SER) +
  scale_linetype_manual(values = E1.LTY, breaks = SER, limits = SER) +
  scale_shape_manual(values = E1.SHP, breaks = SER, limits = SER) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(x = "Horizon H (years)", y = expression(bold("Single-factor ") * bolditalic(hat(beta))[MKT])) +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme_paper + theme(legend.key.width = unit(1.4, "cm"))
ggsave(file.path(out.dir, "figure_E1_beta_profiles.pdf"), p, width = 7, height = 5.2, device = "pdf")
message("wrote ", file.path(out.dir, "figure_E1_beta_profiles.pdf"))
write.tex("figure_E1_beta_profiles.tex", c("\\begin{figure}[htbp]", "\\centering",
  "\\includegraphics[width=0.85\\textwidth]{figure_E1_beta_profiles.pdf}",
  "\\caption{Single-factor $\\hat\\beta_{MKT}$ across horizons: $U_{v,H}$ and L2}", "\\label{fig:u_beta_profiles}",
  "\\par\\vspace{2pt}",
  paste0("{\\scriptsize\\textit{Notes:} ", NOTE.COMMON,
    sprintf("Shaded band: pointwise 95\\%% fixed-$b$ confidence band around $U_{v,H}$ ($\\pm%.3f$ SE for Buyout, $b=%.3f$; $\\pm%.3f$ SE for VC, $b=%.3f$), matching the table inference. ",
      crit.BO, bU.BO, crit.VC, bU.VC),
    "Triangles: $U_{v,H}$ plus the expected estimator gap (mean L2 minus $U_{v,H}$, combining target shift and finite-sample difference) in the data-shaped core Monte Carlo without vintage-level effects ",
    "(EW at 0, 5, 15 and 30 years; FW, with the data's fund-size dispersion, at 0 and 15 years). ",
    "In VC, full-sample L2 uses 27 vintages while $U_{v,H}$ uses the 23 vintages with $n_v\\ge2$; matched-vintage 23-vintage L2 gaps are reported alongside full-sample gaps in Table~\\ref{tab:iv_comparison}. ",
    "The headline estimate is $U_{v,H}$ at $H=0$; the profile across $H$ is a display, not a choice of horizon. Values: Tables~\\ref{tab:u_l2_beta} and~\\ref{tab:u_l2_beta_all}.}"),
  "\\end{figure}"))

cr <- do.call(rbind, lapply(cells, function(x) { for (k in c("type", "weighting", "estimator", "sample", "H")) if (!k %in% names(x)) x[[k]] <- NA; x[, c("exhibit", "type", "weighting", "estimator", "sample", "H", "column", "value", "printed", "source")] }))
write.csv(cr, file.path(out.dir, "exhibits_data.csv"), row.names = FALSE)
message("recorded ", nrow(cr), " numbers (", sum(!is.na(cr$printed)), " printed) in ", file.path(out.dir, "exhibits_data.csv"))

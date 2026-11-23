# Sample robustness table (revision plan, log entry "robustness runs"; reply to
# R2-4.1): North America Buyout and VC, EW and FW, on the main sample
# (complete histories), all funds, and the complete histories without early
# stoppers with terminal NAV above 20% of total value. Run from r_project/:
#   LANG=en_US.UTF-8 Rscript results/robustness_samples.R
# Reads the r5 result rows of the cached results (pipeline version of
# helper/result_cache.R, full sample, from the scenarios without CV): single-factor beta by L2 on
# vintage portfolios and by U_{v,H} (q_v = 1), and (alpha, beta) by L2, at 0, 15
# and 30 years. Writes results/figures_bovc_us_g1/robustness_samples.tex and
# robustness_samples_data.csv (every printed number with its source row).
source("helper/result_cache.R")
RS.SAMPLES <- c(main = "results/data_out_2026-emp-c", allfunds = "results/data_out_2026-emp-c-allfunds", navtrim20 = "results/data_out_2026-emp-c-navtrim20")
RS.DATA <- c(main = "empirical/data_prepared_2026_c", allfunds = "empirical/data_prepared_2026_c_allfunds", navtrim20 = "empirical/data_prepared_2026_c_navtrim20")
RS.LABEL <- c(main = "Complete histories (main)", allfunds = "All funds", navtrim20 = "Without early stoppers, NAV $>20\\%$")
RS.H <- c(0, 180, 360)
out.dir <- "results/figures_bovc_us_g1"

# r5 full-sample rows of one result folder and horizon, matched strictly to current input.md5
rs.rows <- function(sample, folder, m) {
  fl <- list.files(file.path(RS.SAMPLES[[sample]], folder), sprintf("_cvALL_m%d_l0_cached_res\\.csv$", m), full.names = TRUE)
  d <- do.call(rbind.fill.rs, lapply(fl, read.csv, stringsAsFactors = FALSE))
  is.fw <- grepl("_FW_", folder)
  is.vyp <- grepl("_VYP_", folder)
  in.file <- file.path(RS.DATA[[sample]], sprintf("preqin_cashflows_%s%s_North America_NAV.csv", if (is.fw) "FW" else "EW", if (is.vyp) "_VYP" else ""))
  cur.md5 <- unname(tools::md5sum(in.file))
  # r5, current input.md5, full sample, from the no-CV scenario
  d <- d[d$pipeline.version == PIPELINE.VERSION & d$input.md5 == cur.md5 & d$CV.key == "ALL" & grepl("no_cv", d$scenario.id) & d$Type %in% c("BO", "VC"), ]
  key <- paste(d$Type, d$Factor)
  if (anyDuplicated(key)) stop("more than one r5 row per type and model in ", folder, " m", m, " (", sample, ")")
  d
}
rbind.fill.rs <- function(...) { l <- list(...); cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(x) { for (k in setdiff(cols, names(x))) x[[k]] <- NA; x[, cols] })) }

cells <- list()
rec <- function(key, column, value, printed) { cells[[length(cells) + 1]] <<- data.frame(key, column = column, value = value, printed = printed, stringsAsFactors = FALSE); printed }
f3 <- function(x) if (is.na(x)) "--" else formatC(round(x, 3) + 0, format = "f", digits = 3)
est.se <- function(r, par, key, scale = 1) {
  b <- r[[par]] * scale; se <- r[[paste0("SE.", par)]] * scale
  e <- rec(key, par, b, f3(b)); s <- if (isTRUE(r$at.bound)) rec(key, paste0("SE.", par), NA, "--") else rec(key, paste0("SE.", par), se, f3(se))
  sprintf("%s%s (%s)", e, if (isTRUE(r$at.bound)) "$^\\dagger$" else "", s)
}
# funds and vintages (V <= 2010) of a sample, from its VYP file
rs.counts <- function(sample, w, type) {
  v <- read.csv(file.path(RS.DATA[[sample]], sprintf("preqin_cashflows_%s_VYP_North America_NAV.csv", w)))
  v <- v[v$type == type & v$Vintage <= 2010, ]; v <- v[!duplicated(v$Vintage), ]
  c(funds = sum(v$N.funds), vintages = nrow(v))
}
lines <- character(0)
for (type in c("BO", "VC")) for (w in c("EW", "FW")) {
  lines <- c(lines, sprintf("\\multicolumn{8}{l}{\\textit{%s, %s}} \\\\", c(BO = "Buyout", VC = "Venture capital")[[type]], c(EW = "equally weighted", FW = "fund-size weighted")[[w]]))
  for (smp in names(RS.SAMPLES)) {
    n <- rs.counts(smp, w, type)
    key0 <- data.frame(sample = smp, type = type, weighting = w)
    cnt <- rec(cbind(key0, estimator = NA, model = NA, H = NA), "funds", n[["funds"]], format(n[["funds"]]))
    vals <- character(0)
    for (est in c("L2", "U")) {
      folder <- if (est == "L2") sprintf("cache_q_factors_preqin_%s_VYP_North America", w) else sprintf("cache_q_factors_preqin_%s_North America_U_pair", w)
      for (m in RS.H) { r <- rs.rows(smp, folder, m); r <- r[r$Type == type & r$Factor == "MKT", ]; stopifnot(nrow(r) == 1)
        vals <- c(vals, est.se(r, "MKT", cbind(key0, estimator = est, model = "MKT", H = m / 12, source = basename(folder)))) }
    }
    lines <- c(lines, sprintf("\\quad %s & %s & %s \\\\", RS.LABEL[[smp]], cnt, paste(vals, collapse = " & ")))
  }
}
# Panel B: (alpha, beta), L2
linesB <- character(0)
for (type in c("BO", "VC")) for (w in c("EW", "FW")) {
  linesB <- c(linesB, sprintf("\\multicolumn{8}{l}{\\textit{%s, %s}} \\\\", c(BO = "Buyout", VC = "Venture capital")[[type]], c(EW = "equally weighted", FW = "fund-size weighted")[[w]]))
  for (smp in names(RS.SAMPLES)) {
    key0 <- data.frame(sample = smp, type = type, weighting = w)
    folder <- sprintf("cache_q_factors_preqin_%s_VYP_North America", w)
    vals <- character(0)
    for (m in RS.H) { r <- rs.rows(smp, folder, m); r <- r[r$Type == type & r$Factor == "Alpha", ]; stopifnot(nrow(r) == 1)
      k <- cbind(key0, estimator = "L2", model = "Alpha", H = m / 12, source = basename(folder))
      vals <- c(vals, est.se(r, "MKT", k), est.se(r, "Alpha", k, scale = 100)) }
    linesB <- c(linesB, sprintf("\\quad %s & %s & \\\\", RS.LABEL[[smp]], paste(vals, collapse = " & ")))
  }
}
# Panel C: (alpha, beta), U_{v,H} (q_v = 1)
linesC <- character(0)
for (type in c("BO", "VC")) for (w in c("EW", "FW")) {
  linesC <- c(linesC, sprintf("\\multicolumn{8}{l}{\\textit{%s, %s}} \\\\", c(BO = "Buyout", VC = "Venture capital")[[type]], c(EW = "equally weighted", FW = "fund-size weighted")[[w]]))
  for (smp in names(RS.SAMPLES)) {
    key0 <- data.frame(sample = smp, type = type, weighting = w)
    folder <- sprintf("cache_q_factors_preqin_%s_North America_U_pair", w)
    vals <- character(0)
    for (m in RS.H) { r <- rs.rows(smp, folder, m); r <- r[r$Type == type & r$Factor == "Alpha", ]; stopifnot(nrow(r) == 1)
      k <- cbind(key0, estimator = "U", model = "Alpha", H = m / 12, source = basename(folder))
      vals <- c(vals, est.se(r, "MKT", k), est.se(r, "Alpha", k, scale = 100)) }
    linesC <- c(linesC, sprintf("\\quad %s & %s & \\\\", RS.LABEL[[smp]], paste(vals, collapse = " & ")))
  }
}
tex <- c(sprintf("%% Auto-generated by r_project/results/robustness_samples.R at commit %s%s", system("git rev-parse --short HEAD", intern = TRUE),
    if (length(system("git status --porcelain --untracked-files=no", intern = TRUE))) " (uncommitted changes)" else ""),
  "\\begin{table}[htbp]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{0.92}",
  "\\caption{Sample robustness: complete histories, all funds, and without early stoppers with large terminal NAV}", "\\label{tab:robustness_samples}",
  "\\resizebox{\\linewidth}{!}{%", "\\begin{tabular}{l r rrr rrr}", "\\toprule",
  "\\multicolumn{8}{l}{\\textit{Panel A: single-factor $\\hat\\beta_{MKT}$ (SE)}} \\\\",
  " & & \\multicolumn{3}{c}{L2, vintage portfolios} & \\multicolumn{3}{c}{$U_{v,H}$, $q_v=1$} \\\\",
  "\\cmidrule(lr){3-5} \\cmidrule(lr){6-8}", "Sample & Funds & $H=0$ & 15 yr & 30 yr & $H=0$ & 15 yr & 30 yr \\\\", "\\midrule",
  lines, "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Panel B: $(\\alpha,\\beta)$ model, L2 on vintage portfolios: $\\hat\\beta_{MKT}$ (SE) and $\\hat\\alpha$ in monthly \\% (SE)}} \\\\",
  " & \\multicolumn{2}{c}{$H=0$} & \\multicolumn{2}{c}{15 yr} & \\multicolumn{2}{c}{30 yr} & \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}", "Sample & $\\hat\\beta$ & $\\hat\\alpha$ & $\\hat\\beta$ & $\\hat\\alpha$ & $\\hat\\beta$ & $\\hat\\alpha$ & \\\\", "\\midrule",
  linesB, "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Panel C: $(\\alpha,\\beta)$ model, $U_{v,H}$ ($q_v=1$): $\\hat\\beta_{MKT}$ (SE) and $\\hat\\alpha$ in monthly \\% (SE)}} \\\\",
  " & \\multicolumn{2}{c}{$H=0$} & \\multicolumn{2}{c}{15 yr} & \\multicolumn{2}{c}{30 yr} & \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}", "Sample & $\\hat\\beta$ & $\\hat\\alpha$ & $\\hat\\beta$ & $\\hat\\alpha$ & $\\hat\\beta$ & $\\hat\\alpha$ & \\\\", "\\midrule",
  linesC, "\\bottomrule", "\\end{tabular}%", "}", "\\par\\vspace{2pt}",
  paste0("{\\scriptsize\\textit{Notes:} North America, vintages up to 2010, the estimators of the main tables (r5 domain, Bartlett sandwich SE with $S=\\lceil1.3\\sqrt V\\rceil$, no cross-validation). ",
    "\\textit{Complete histories}: the main sample (a fund's last record lies more than 5 years after its vintage year starts, or it still reports at the data end). \\textit{All funds}: adds the truncated histories. ",
    "\\textit{Without early stoppers}: the main sample without funds that stop reporting after age 5 and before 2021-03-31 with a last reported NAV above 20\\% of total value (distributions plus that NAV). ",
    "In Panel~C, VC $U_{v,H}$ sits on the $+2\\%$ monthly $\\alpha$ cap at $H=0$ across all three samples while FW $\\hat\\beta$ shifts noticeably under the NAV-threshold trim. ",
    "Funds: funds in the sample's vintage portfolios. $^\\dagger$: estimate on the boundary of the admissible domain, no SE.}"),
  "\\end{table}")
writeLines(tex, file.path(out.dir, "robustness_samples.tex"))
write.csv(do.call(rbind.fill.rs, cells), file.path(out.dir, "robustness_samples_data.csv"), row.names = FALSE)
cat("wrote", file.path(out.dir, "robustness_samples.tex"), "with", length(cells), "recorded numbers\n")

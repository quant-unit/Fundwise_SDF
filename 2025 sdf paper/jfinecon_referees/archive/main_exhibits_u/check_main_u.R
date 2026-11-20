# Independent check of results/main_exhibits_u.R output (log entry "U-first main exhibits: built"). Run from r_project/:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/main_exhibits_u/check_main_u.R"
src <- "results/figures_bovc_us_g1"; out <- file.path(src, "main_u")
x <- read.csv(file.path(out, "exhibits_data.csv"), stringsAsFactors = FALSE, colClasses = c(printed = "character"))
cat("recorded", nrow(x), "| printed", sum(!is.na(x$printed)), "\n")
# 1) every value against its source, looked up independently
L <- rbind(read.csv(file.path(src, "empirical_BO_US_MKT_data.csv")), read.csv(file.path(src, "empirical_VC_US_MKT_data.csv"))); L <- L[L$method == "Asymptotic", ]
U <- rbind(read.csv(file.path(src, "upair_BO_US_data.csv")), read.csv(file.path(src, "upair_VC_US_data.csv"))); U <- U[U$vintage.weights == "equal", ]
P <- read.csv(file.path(src, "iv_comparison_rows.csv")); X <- read.csv(file.path(src, "vc_alpha_without_1993_1997.csv"))
look <- function(r) {
  m <- round(r$H * 12); ex <- r$exhibit; col <- r$column
  model <- if (ex == "E2") "Alpha" else "MKT"
  base <- sub("\\.(se|p)$", "", col); what <- if (grepl("\\.se$", col)) "se" else if (grepl("\\.p$", col)) "p" else "est"
  if (ex == "OA.vc") { y <- X[X$weighting == r$weighting & X$estimator == r$estimator & X$sample == r$sample & X$max.month == m, ]
    return(if (col == "beta") y$beta else 100 * y$alpha) }
  if (col == "pred.l2u") return(P$pred.l2u[P$type == r$type & P$weighting == r$weighting & P$max.month == m])
  if (col == "U.se") base <- "U.beta"
  par <- if (grepl("alpha", base)) "Alpha" else "MKT"; sc <- if (par == "Alpha") 100 else 1
  if (startsWith(base, "U")) { y <- U[U$Type == r$type & U$weighting == r$weighting & U$max.month == m & U$model == model & U$parameter == par, ]
    v <- switch(what, est = y$estimate * sc, se = if (y$at.bound) NA else y$SE.A * sc, p = y$p) }
  else { y <- L[L$Type == r$type & L$weighting == r$weighting & L$max.month == m & L$Factor == model, ]
    v <- switch(what, est = (if (par == "MKT") y$MKT else y$Coef) * sc, se = if (y$at.bound) NA else (if (par == "MKT") y$SE.MKT else y$SE.Coef) * sc, p = if (par == "MKT") y$p.MKT else y$p.Coef) }
  stopifnot(length(v) == 1); v
}
ref <- vapply(seq_len(nrow(x)), function(i) look(x[i, ]), 0)
d <- abs(ref - x$value); d[is.na(ref) & is.na(x$value)] <- 0
cat("1) values vs sources: max |diff|", max(d), "| NA mismatches", sum(is.na(ref) != is.na(x$value)), "\n")
# 2) printed strings: 3-decimal rounding of the value; stars from p
f3 <- function(v) ifelse(is.na(v), "--", formatC(round(v, 3) + 0, format = "f", digits = 3))
num <- !is.na(x$printed) & !grepl("\\.p$", x$column)
cat("2a) printed = f3(value):", all(x$printed[num] == f3(x$value[num])), "\n")
pp <- !is.na(x$printed) & grepl("\\.p$", x$column)
st <- ifelse(x$value[pp] < 0.001, "***", ifelse(x$value[pp] < 0.01, "**", ifelse(x$value[pp] < 0.05, "*", "")))
cat("2b) stars = p thresholds:", all(st == x$printed[pp]), "(", sum(pp), "markers;", sum(st != ""), "non-empty )\n")
# 3) the numbers appear in the tex files in recorded order
for (f in c("table_E1_beta.tex", "table_E1_beta_all.tex", "table_E2_alpha_beta.tex", "table_OA_vc_without_1993_1997.tex")) {
  ex <- c(table_E1_beta.tex = "E1", table_E1_beta_all.tex = "E1.all", table_E2_alpha_beta.tex = "E2", table_OA_vc_without_1993_1997.tex = "OA.vc")[[f]]
  body <- readLines(file.path(out, f)); body <- body[grepl("^\\\\quad", body)]
  toks <- unlist(regmatches(body, gregexpr("-?[0-9]+\\.[0-9]{3}|--", body)))
  want <- x$printed[x$exhibit == ex & !grepl("\\.p$", x$column)]
  cat(sprintf("3) %-36s %3d numbers in tex, %3d recorded, same sequence: %s\n", f, length(toks), length(want), identical(toks, want)))
}
# 4) cross-check against the other exhibits: iv_comparison rows and robustness table (main sample)
R <- read.csv(file.path(src, "robustness_samples_data.csv"), stringsAsFactors = FALSE)
e1 <- x[x$exhibit == "E1.all" & x$column %in% c("U.beta", "L2.beta"), ]
cmp <- merge(transform(e1, estimator = ifelse(column == "U.beta", "U", "L2")),
  R[R$sample == "main" & R$column == "MKT" & R$model == "MKT", c("type", "weighting", "estimator", "H", "value")], by = c("type", "weighting", "estimator", "H"))
cat("4a) E1 vs robustness table (main sample, MKT):", nrow(cmp), "cells, max |diff|", max(abs(cmp$value.x - cmp$value.y)), "\n")
cmp2 <- merge(e1[e1$column == "U.beta", ], P[, c("type", "weighting", "H", "u")], by = c("type", "weighting", "H"))
cat("4b) E1 U vs iv_comparison:", nrow(cmp2), "cells, max |diff|", max(abs(cmp2$value - cmp2$u)), "\n")
b <- x[x$exhibit %in% c("E1.all", "E2", "OA.vc") & x$printed %in% "--", ]
cat("5) cells printed as '--' (bound SEs):\n"); print(b[, c("exhibit", "type", "weighting", "estimator", "sample", "H", "column")], row.names = FALSE)

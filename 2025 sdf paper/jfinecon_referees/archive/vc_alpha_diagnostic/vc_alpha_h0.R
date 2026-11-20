# Diagnostic of the (alpha, beta) model at H = 0 (revision plan, Section I,
# "Diagnostics not run, optional"; log entry "VC alpha at H = 0: diagnostic").
# Why is the VC (alpha, beta) estimate at H = 0 on the alpha bound (+2%/month)
# under L2 on vintage portfolios and under U_{v,H}, while Buyout is interior?
# Run from r_project/:
#   LANG=en_US.UTF-8 Rscript "../2025 sdf paper/jfinecon_referees/archive/vc_alpha_diagnostic/vc_alpha_h0.R"
# Uses the Step 4 fast criterion (helper/objective.R; equal to the pipeline's,
# test T1) on the prepared North America files; prints aggregates only
# (estimates, criterion values, vintage years and fund counts).
# (A) reproduce the r5 estimates at H = 0 (G1 tables); (B) profile the
# criterion in alpha from -2% to +6% per month, minimizing over beta within the
# rest of the domain (beta in [-1, 4], gross floor 0.08): does it keep falling
# beyond the bound?; (C) leave-one-vintage-out refits within Theta: which
# vintages hold alpha at the bound?; (D) each vintage's share of the
# alpha-derivative of the criterion at the bound estimate; (E) refits without
# the mid-1990s block that (D) points to; (F) the same comparison at H > 0
# (2.5-30 years), full-sample fits checked against the G1 tables first.
# Part (F) also writes its fits to results/figures_bovc_us_g1/
# vc_alpha_without_1993_1997.csv (estimates only), read by
# results/main_exhibits_u.R for the Online Appendix table (log entry
# "U-first main exhibits: built").
source("simulation/mc_driver.R")
options(width = 200)
pub <- local({ p <- read.csv("empirical/data_prepared_2026/q_factors.csv"); p$Date <- as.Date(p$Date)
  p$Date <- as.Date(format(p$Date, "%Y-%m-01")); p$Date <- seq.Date(p$Date[1], by = "month", length.out = nrow(p))
  p$Date <- as.Date(format(p$Date + 32, "%Y-%m-01")) - 1; p$Alpha <- 1; p })
fc <- unique(c("RF", setdiff(colnames(pub), "Date")))
load.fp <- function(w, vyp, type) {
  d <- read.csv(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s%s_North America_NAV.csv", w, if (vyp) "_VYP" else ""), stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date); d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  d <- pad.all.units(d, "first.cf", max(d$Date))
  d <- merge(d, pub, by = "Date", all.x = TRUE); d <- d[d$type == type & d$Vintage <= 2010, ]
  d$Fund.ID <- droplevels(as.factor(as.character(d$Fund.ID)))
  prepare.fast(d, fc)
}
EST <- list(L2 = MC.ESTIMATORS$L2.vyp, U = MC.ESTIMATORS$U.eq)
# criterion as the pipeline evaluates it: for FW, U normalizes each fund pair by
# the funds' sizes (estim_model_optimized.R, use.size); the Step 4 criterion
# (mc.criterion) assumes size 1, which holds for EW and for simulated funds
crit <- function(k, fp, m = 0) {
  e <- EST[[sub(".* ", "", k)]]; sz <- if (e$criterion == "U" && grepl(" FW ", k)) fp$size else rep(1, fp$n.units)
  function(p) { u <- fast.unit.errors(p, fp, m); if (any(!is.finite(u))) return(MC.SETTINGS$wall)
    v <- if (e$criterion == "L2") sqrt(sum(u^2)) / length(u) else fast.u.criterion(u, fp, sz, e$vintage.weights)
    if (is.finite(v)) v else MC.SETTINGS$wall }
}
# (alpha, beta) fit over Theta as mc.fit() does, with that criterion
fit2 <- function(k, fp, alpha.bounds = MC.SETTINGS$alpha.bounds, m = 0) {
  MC.SETTINGS$alpha.bounds <- alpha.bounds; par <- c(MKT = 1, Alpha = 0)
  con <- mc.domain(fp, names(par)); ps <- optim.parscale(par)
  r <- fit.admissible(crit(k, fp, m), admissible.starts(par, con, ps), con, ps, smooth = TRUE)
  act <- active.constraints(r$par, con)
  data.frame(MKT = r$par[["MKT"]], Coef = r$par[["Alpha"]], value = r$value, at.bound = length(act) > 0, active = paste(act, collapse = ","))
}
cases <- expand.grid(type = c("VC", "BO"), w = c("EW", "FW"), est = c("L2", "U"), stringsAsFactors = FALSE)
fps <- list()
for (i in seq_len(nrow(cases))) { k <- paste(cases[i, ], collapse = " "); fps[[k]] <- load.fp(cases$w[i], cases$est[i] == "L2", cases$type[i]) }

cat("\n(A) r5 estimates at H = 0, (alpha, beta), domain alpha in [-2%, 2%]\n")
fitA <- list()
for (k in names(fps)) { r <- fit2(k, fps[[k]]); fitA[[k]] <- r
  cat(sprintf("  %-9s beta %7.3f  alpha %6.3f%%  criterion %.6g  on bound: %s %s\n", k, r$MKT, 100 * r$Coef, r$value, r$at.bound, r$active)) }

cat("\n(B) profile in alpha (beta minimized over [-1, 4] and the gross floor): criterion minus its value at the in-domain estimate, in % of |that value| (negative = lower)\n")
agrid <- c(-0.02, -0.01, 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06)
prof <- function(k, fp, a) {
  f <- crit(k, fp)
  rows <- min(fp$grids$k0):nrow(fp$X); mo <- data.frame(Date = rows, fp$X[rows, , drop = FALSE]); mo$RF <- mo$RF + a # gross = 1 + rf + alpha + beta MKT
  con <- admissible.constraints(mo, "MKT", MC.SETTINGS$gross.floor, MC.SETTINGS$mkt.bounds, c(-Inf, Inf)); br <- feasible.interval(con)
  if (br[1] >= br[2]) return(c(NA, NA))
  o <- optimize(function(b) f(c(MKT = b, Alpha = a)), br); c(o$objective, o$minimum)
}
cat(sprintf("  %-9s %s\n", "alpha(%)", paste(sprintf("%7.1f", 100 * agrid), collapse = "")))
for (k in names(fps)) { p <- sapply(agrid, function(a) prof(k, fps[[k]], a))
  rel <- 100 * (p[1, ] - fitA[[k]]$value) / abs(fitA[[k]]$value)
  cat(sprintf("  %-9s %s\n", k, paste(sprintf("%7.0f", rel), collapse = "")))
  cat(sprintf("  %-9s %s   (beta at the profile minimum)\n", "", paste(sprintf("%7.2f", p[2, ]), collapse = ""))) }

cat("\n(C) leave-one-vintage-out refits at H = 0 within Theta (VC); vintages whose removal brings alpha inside the bound\n")
for (k in grep("^VC", names(fps), value = TRUE)) { fp <- fps[[k]]
  vs <- sort(unique(fp$vintage)); nf <- tabulate(match(fp$vintage, vs), length(vs))
  r <- do.call(rbind, lapply(vs, function(v) { x <- fit2(k, fast.subset(fp, fp$vintage != v)); data.frame(v = v, beta = x$MKT, alpha = 100 * x$Coef, bound = x$at.bound, active = x$active) }))
  inside <- r[!r$bound, ]
  cat(sprintf("  %-9s %d vintages; on the bound after dropping %d of them; interior after dropping: %s\n", k, length(vs), sum(r$bound),
    if (nrow(inside)) paste(sprintf("%d (beta %.2f, alpha %.2f%%)", inside$v, inside$beta, inside$alpha), collapse = "; ") else "none"))
  cat(sprintf("  %-9s funds per vintage (U units; VYP: 1 unit each): %s\n", "", paste(sprintf("%d:%d", vs, nf), collapse = " ")))
}

cat("\n(D) vintages pushing alpha up at the VC estimate: each vintage's alpha-derivative of its criterion term, in % of the sum of all negative terms\n")
for (k in grep("^VC", names(fps), value = TRUE)) { e <- EST[[sub(".* ", "", k)]]; fp <- fps[[k]]; th <- c(MKT = fitA[[k]]$MKT, Alpha = fitA[[k]]$Coef)
  sz.all <- if (e$criterion == "U" && grepl(" FW ", k)) fp$size else rep(1, fp$n.units)
  vs <- sort(unique(fp$vintage))
  g <- vapply(vs, function(v) { s <- fast.subset(fp, fp$vintage == v); h <- 1e-6
    # vintage v's term: e_v^2 (L2, up to the common monotone transform) or U_v (equal vintage weights; vintages with one fund drop out)
    lv <- function(p) { u <- fast.unit.errors(p, s, 0); if (e$criterion == "L2") sum(u^2) else if (s$n.units < 2) 0 else fast.u.criterion(u, s, sz.all[fp$vintage == v]) }
    (lv(th + c(0, h)) - lv(th - c(0, h))) / (2 * h) }, 0)
  neg <- g < 0; sh <- ifelse(neg, g / sum(g[neg]), 0); o <- order(sh, decreasing = TRUE)
  cat(sprintf("  %-9s derivative of the criterion in alpha: %s (sum of pushing-up terms %.3g, of pushing-down terms %.3g); vintages pushing up: %d; largest: %s\n", k,
    if (sum(g) < 0) "negative, the criterion falls with alpha" else "positive", sum(g[neg]), sum(g[!neg]), sum(neg),
    paste(sprintf("%d %.0f%%", vs[o][1:5], 100 * sh[o][1:5]), collapse = ", ")))
}

cat("\n(E) VC without the mid-1990s block that pushes alpha up in (D) (1993-1997 jointly; also 1994-1996), refit within Theta\n")
for (k in grep("^VC", names(fps), value = TRUE)) for (blk in list(1993:1997, 1994:1996)) { fp <- fps[[k]]
  x <- fit2(k, fast.subset(fp, !(fp$vintage %in% blk)))
  cat(sprintf("  %-9s without %d-%d: beta %6.3f  alpha %6.3f%%  on bound: %s %s\n", k, min(blk), max(blk), x$MKT, 100 * x$Coef, x$at.bound, x$active)) }

cat("\n(F) VC (alpha, beta) by horizon, full sample (checked against the G1 tables) and without 1993-1997; beta / alpha in %/month, * = on a bound\n")
g1 <- read.csv("results/figures_bovc_us_g1/empirical_VC_US_Alpha_data.csv", stringsAsFactors = FALSE); g1 <- g1[g1$Factor == "Alpha" & g1$method == "Asymptotic", ]
gu <- read.csv("results/figures_bovc_us_g1/upair_VC_US_data.csv", stringsAsFactors = FALSE); gu <- gu[gu$model == "Alpha" & gu$vintage.weights == "equal", ]
HF <- c(0, 30, 60, 120, 180, 240, 360); fmt <- function(x) sprintf("%5.2f/%5.2f%s", x$MKT, 100 * x$Coef, if (x$at.bound) "*" else " ")
worst <- 0; F.rows <- list()
cat(sprintf("  %-26s %s\n", "H (years)", paste(sprintf("%13s", HF / 12), collapse = "")))
for (k in grep("^VC", names(fps), value = TRUE)) { fp <- fps[[k]]; w <- strsplit(k, " ")[[1]][2]
  full <- lapply(HF, function(m) fit2(k, fp, m = m)); red <- lapply(HF, function(m) fit2(k, fast.subset(fp, !(fp$vintage %in% 1993:1997)), m = m))
  ref <- vapply(seq_along(HF), function(i) { m <- HF[i]
    if (grepl(" L2$", k)) { y <- g1[g1$weighting == w & g1$max.month == m, ]; c(y$MKT, y$Coef) } else { y <- gu[gu$weighting == w & gu$max.month == m, ]; c(y$estimate[y$parameter == "MKT"], y$estimate[y$parameter == "Alpha"]) } }, c(0, 0))
  worst <- max(worst, abs(sapply(full, `[[`, "MKT") - ref[1, ]), abs(sapply(full, `[[`, "Coef") - ref[2, ]) * 100)
  fp.red <- fast.subset(fp, !(fp$vintage %in% 1993:1997)); both <- c(full, red)
  F.rows[[length(F.rows) + 1]] <- data.frame(type = "VC", weighting = w, estimator = sub(".* ", "", k), sample = rep(c("full", "without 1993-1997"), each = length(HF)),
    n.vintages = rep(c(length(unique(fp$vintage)), length(unique(fp.red$vintage))), each = length(HF)), max.month = rep(HF, 2), H = rep(HF / 12, 2),
    beta = sapply(both, `[[`, "MKT"), alpha = sapply(both, `[[`, "Coef"), criterion = sapply(both, `[[`, "value"), at.bound = sapply(both, `[[`, "at.bound"),
    active = sapply(both, `[[`, "active"), stringsAsFactors = FALSE)
  cat(sprintf("  %-26s %s\n", paste(k, "full"), paste(sprintf("%13s", sapply(full, fmt)), collapse = "")))
  cat(sprintf("  %-26s %s\n", paste(k, "without 1993-97"), paste(sprintf("%13s", sapply(red, fmt)), collapse = ""))) }
cat(sprintf("  full-sample fits against the G1 tables: largest difference %.1e (beta; alpha in %%)\n", worst))
stopifnot(worst < 1e-8)
write.csv(do.call(rbind, F.rows), "results/figures_bovc_us_g1/vc_alpha_without_1993_1997.csv", row.names = FALSE)
cat("  wrote results/figures_bovc_us_g1/vc_alpha_without_1993_1997.csv\n")

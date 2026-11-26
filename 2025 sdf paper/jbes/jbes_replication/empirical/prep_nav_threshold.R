# Robustness sample: the complete-history main sample without funds that stop
# reporting early with a large terminal NAV (revision plan, Section B,
# "Terminal NAV"; log entries "robustness runs: scope" and "robustness runs").
# Run from r_project/:  LANG=en_US.UTF-8 Rscript empirical/prep_nav_threshold.R [threshold]
#
# A fund is dropped if it stops reporting after age 5 (Age.Last > 5) and before
# 2021-03-31, and its terminal NAV (last reported NAV) exceeds `threshold`
# (default 0.20) of its total value, sum of distributions (DIS, per-row flows)
# plus terminal NAV. The rule is applied to every strategy type; the paper uses
# Buyout and VC. Input: the prepared main fund-level files (EW and FW, North
# America) in empirical/data_prepared_2026_c/, so the licensed raw export is not
# read again. The reduced sample is re-pooled with pool.all.vintages() (helper/
# pool_vintage.R): N_v and the vintage commitments are those of the kept funds.
# Output: empirical/data_prepared_2026_c_navtrim<100 x threshold>/ with the
# fund-level and VYP files under the main sample's names (git-ignored).
# Guard: before writing, the pools of the full sample (nothing dropped) must
# equal the main VYP files (relative difference below 1e-12, the CSV round
# trip); otherwise the script stops.
source("helper/pool_vintage.R")

NT.MAIN <- "empirical/data_prepared_2026_c"
NT.STOP.BEFORE <- as.Date("2021-03-31")
NT.MIN.AGE <- 5
nt.file <- function(folder, weighting, vyp) file.path(folder, sprintf("preqin_cashflows_%s%s_North America_NAV.csv", weighting, if (vyp) "_VYP" else ""))

nt.read <- function(weighting) {
  d <- read.csv(nt.file(NT.MAIN, weighting, FALSE), stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date); d$Last.Date <- as.Date(d$Last.Date)
  d
}

# One row per fund and type (a fund can belong to several types, e.g. BO,
# GroBO and PE, with the same records in each): terminal NAV, distributions,
# share, and the drop flag, which must then be the same in every type
nt.flags <- function(d, threshold) {
  f <- do.call(rbind, lapply(split(d, paste(d$type, d$Fund.ID), drop = TRUE), function(x) {
    x <- x[order(x$Date), ]; v <- x$NAV[!is.na(x$NAV)]
    data.frame(Fund.ID = x$Fund.ID[1], type = x$type[1], Vintage = x$Vintage[1], last = x$Last.Date[1], age = x$Age.Last[1],
               nav = if (length(v)) v[length(v)] else NA_real_, dis = sum(x$DIS, na.rm = TRUE), stringsAsFactors = FALSE)
  }))
  f$share <- f$nav / (f$dis + f$nav)
  f$drop <- f$age > NT.MIN.AGE & f$last < NT.STOP.BEFORE & !is.na(f$share) & f$share > threshold
  rownames(f) <- NULL
  per.fund <- tapply(f$drop, f$Fund.ID, function(x) length(unique(x)))
  if (any(per.fund > 1)) stop("funds flagged differently across types: ", paste(names(per.fund)[per.fund > 1], collapse = ", "))
  f
}

# Pools of a fund-level sample, in the order of the main VYP file (types as
# there, then vintage and date)
nt.pool <- function(d, weighting, type.order) {
  v <- do.call(rbind, lapply(split(d, d$type, drop = TRUE), pool.all.vintages, fund.size.weighting = weighting == "FW"))
  v <- v[order(match(v$type, type.order), v$Vintage, v$Date), ]
  rownames(v) <- NULL
  v
}

nt.rel.diff <- function(a, b) max(abs(a - b)) / max(abs(b))

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  threshold <- if (length(args)) as.numeric(args[1]) else 0.20
  stopifnot(is.finite(threshold), threshold > 0, threshold < 1)
  out <- paste0(NT.MAIN, "_navtrim", round(100 * threshold))
  dir.create(out, showWarnings = FALSE)
  counts <- list()
  for (w in c("EW", "FW")) {
    d <- nt.read(w)
    main.vyp <- read.csv(nt.file(NT.MAIN, w, TRUE), stringsAsFactors = FALSE); main.vyp$Date <- as.Date(main.vyp$Date)
    type.order <- unique(main.vyp$type)
    # guard: nothing dropped -> the main pools
    full <- nt.pool(d, w, type.order)
    same <- identical(dim(full), dim(main.vyp)) && identical(full$type, main.vyp$type) && identical(full$Vintage, main.vyp$Vintage) &&
      identical(full$Date, main.vyp$Date) && identical(full$N.funds, main.vyp$N.funds) &&
      all(vapply(c("CF", "NAV", "CON", "DIS", "Commitment"), function(k) nt.rel.diff(full[[k]], main.vyp[[k]]) < 1e-12, TRUE))
    if (!same) stop(w, ": the pools of the full sample do not reproduce ", nt.file(NT.MAIN, w, TRUE))
    cat(w, ": guard passed (full-sample pools = main VYP file; max relative difference ",
        signif(max(vapply(c("CF", "NAV", "CON", "DIS", "Commitment"), function(k) nt.rel.diff(full[[k]], main.vyp[[k]]), 0)), 2), ")\n", sep = "")
    f <- nt.flags(d, threshold)
    kept <- d[!(d$Fund.ID %in% f$Fund.ID[f$drop]), ] # flag is per fund (checked above)
    write.csv(kept, nt.file(out, w, FALSE), row.names = FALSE)
    write.csv(nt.pool(kept, w, type.order), nt.file(out, w, TRUE), row.names = FALSE)
    counts[[w]] <- f
    cat(sprintf("%s: %d of %d funds dropped (all types; vintages <= 2010: BO %d of %d, VC %d of %d)\n", w, length(unique(f$Fund.ID[f$drop])), length(unique(f$Fund.ID)),
      sum(f$drop & f$type == "BO" & f$Vintage <= 2010), sum(f$type == "BO" & f$Vintage <= 2010),
      sum(f$drop & f$type == "VC" & f$Vintage <= 2010), sum(f$type == "VC" & f$Vintage <= 2010)))
  }
  # the rule does not depend on the weighting: the same funds are dropped
  stopifnot(identical(sort(unique(counts$EW$Fund.ID[counts$EW$drop])), sort(unique(counts$FW$Fund.ID[counts$FW$drop]))))
  cat("same funds dropped under EW and FW; written to", out, "\n")
}

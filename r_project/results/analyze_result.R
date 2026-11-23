# analyze results -----

if (!exists("source.internally", envir = .GlobalEnv)) {
  source.internally <- TRUE
}

if (source.internally) {
  if (sys.nframe() == 0L) rm(list = ls())

  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  getwd()

  prefix <- "q_factors_preqin_"
  suffix <- "FW_VYP"
  data.out.folder <- "results/data_out_2026-emp-F-max-vin-2019"
  data.out.folder <- "results/data_out_2026_02_18"
}

dir.cache <- here(paste0(data.out.folder, "/cache_", prefix, suffix))
dir.cache

# Read cached results through the shared filter (helper/result_cache.R): only
# rows of the current pipeline version, from a single input file; stale or
# pre-repair rows in the folder are ignored (review fix 2026-11-16)
source(here("helper/result_cache.R"))
source(here("helper/result_reporting.R"))
source(here("helper/fixedb.R")) # fixed-b p-values (decision 2026-11-17)
df.f <- read.result.cache(dir.cache)
if (nrow(df.f) == 0) stop("No current cached results; rerun estimation before analysis: ", dir.cache)
if (nrow(df.f) > 0) df.f$CV.key <- as.character(df.f$CV.key)
rownames(df.f) <- NULL
df.f$X <- NULL
if (nrow(df.f) > 0) {
  df.f$RF <- 1
}
# df.f[is.na(df.f)] <- 0
# df.f$validation.error[is.na(df.f$validation.error)] <- Inf

# Determine public factor file from prefix
# prefix is e.g. "q_factors_preqin_" or "ff3_factors_preqin_"
if (grepl("^ff3_factors", prefix)) {
  public_factor_file <- "ff3_factors.csv"
} else {
  public_factor_file <- "q_factors.csv"
}
file_path <- here("empirical", "data_prepared_2026", public_factor_file)

df.q <- read.csv(file_path)
df.q$Date <- as.Date(df.q$Date)
df.q$Alpha <- 1 # Alpha = 1 is the constant term in the SDF

rbind.all.columns <- function(x, y) {
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))

  x[, c(as.character(y.diff))] <- NA

  y[, c(as.character(x.diff))] <- NA

  return(rbind(x, y))
}


# summarize cross-validation coefs -----
require(xtable)
q.factors <- levels(as.factor(df.f$Factor))

cv.res <- list()
for (Factor in q.factors) {
  # Select columns without duplicating MKT when Factor == "MKT"
  if (Factor == "MKT") {
    select_cols <- c("MKT", "Type", "max.month", "validation.error", "validation.error.H0")
  } else if (Factor == "ALL" || grepl("\\+", Factor)) {
    available_factors <- setdiff(intersect(c("Alpha", "SMB", "HML", "ME", "ROE", "IA", "EG"), colnames(df.f)), "MKT")
    select_cols <- c("MKT", available_factors, "Type", "max.month", "validation.error", "validation.error.H0")
  } else {
    select_cols <- c("MKT", Factor, "Type", "max.month", "validation.error", "validation.error.H0")
  }
  # Only select columns that actually exist in df.f (e.g. no validation.error if no CV runs)
  select_cols <- intersect(select_cols, colnames(df.f))

  fold.rows <- df.f[(df.f$CV.key != "ALL") & (df.f$Factor == Factor), , drop = FALSE]
  # Repeated executions of a fold are not additional independent folds.
  fold.rows <- fold.rows[order(fold.rows$datetime, decreasing = TRUE), , drop = FALSE]
  fold.rows <- fold.rows[!duplicated(paste(fold.rows$Type, fold.rows$max.month, fold.rows$CV.key, fold.rows$lambda)), , drop = FALSE]
  df <- fold.rows[, select_cols, drop = FALSE]
  if (nrow(df) == 0) next
  ddf <- df[, !is.na(df[1, ])]
  df$id <- paste0(df$Type, "_", df$max.month)
  provenance <- aggregate.result.provenance(fold.rows, df$id)
  df$Type <- NULL
  df$max.month <- NULL

  if (is.data.frame(df)) {
    m <- aggregate(. ~ id, mean, data = df)
    s <- aggregate(. ~ id, sd, data = df)
    s["validation.error"] <- NULL
    s["validation.error.H0"] <- NULL # B8: common H = 0 held-out loss (fold mean kept in m)
  }
  n <- ncol(s)
  a <- 2
  b <- 3
  if (Factor == "MKT") {
    # One-factor model: Coef = MKT, SE.Coef = SE.MKT
    m$Coef <- m$MKT
    s$Coef <- s$MKT
    names(s)[a:ncol(s)] <- paste0("SE.", names(s)[a:ncol(s)])
  } else if (Factor == "ALL" || grepl("\\+", Factor)) {
    names(m)[a:n] <- paste0(names(m)[a:n])
    names(s)[a:n] <- paste0("SE.", names(s)[a:n])
  } else if (n == b) {
    names(m)[b] <- "Coef"
    names(s)[b] <- "Coef"
    names(m)[a:n] <- paste0(names(m)[a:n])
    names(s)[a:n] <- paste0("SE.", names(s)[a:n])
  }

  df <- merge(m, s, by = "id")
  df <- merge(df, provenance, by = "id")
  df$Factor <- Factor
  cv.res[[Factor]] <- df
}

if (length(cv.res) > 0) {
  df.cv <- data.frame(Reduce(rbind.all.columns, cv.res))
  df.cv$id <- as.character(df.cv$id)
  for (i in seq_len(nrow(df.cv))) {
    df.cv[i, "Type"] <- strsplit(df.cv$id, "_")[[i]][1]
    df.cv[i, "max.month"] <- strsplit(df.cv$id, "_")[[i]][2]
  }
  cv_select_cols <- c("Type", "max.month", "MKT", "SE.MKT", "Factor", "validation.error", "validation.error.H0")
  if (any(df.cv$Factor != "ALL" & !grepl("\\+", df.cv$Factor))) {
    cv_select_cols <- c(cv_select_cols, "Coef", "SE.Coef")
  }
  if (any(df.cv$Factor == "ALL" | grepl("\\+", df.cv$Factor))) {
    available_factors <- setdiff(intersect(c("Alpha", "SMB", "HML", "ME", "ROE", "IA", "EG"), colnames(df.cv)), "MKT")
    for (f in available_factors) {
      cv_select_cols <- c(cv_select_cols, f, paste0("SE.", f))
    }
  }
  cv_select_cols <- c(cv_select_cols, RESULT.PROVENANCE.COLUMNS, "at.bound", "boundary.folds")
  cv_select_cols <- intersect(cv_select_cols, colnames(df.cv))
  
  for (col in cv_select_cols) {
      if (!col %in% colnames(df.cv)) df.cv[[col]] <- NA
  }
  df.cv <- df.cv[, cv_select_cols]
  df.cv <- df.cv[order(df.cv$Type, df.cv$Factor, as.numeric(df.cv$max.month)), ]

  # Filter horizons: use caller-defined max.months.to.keep if available, otherwise keep all
  if (!exists("max.months.to.keep", envir = .GlobalEnv)) {
    max.months.to.keep <- unique(df.cv$max.month)
  }
  df.cv <- df.cv[df.cv$max.month %in% as.character(max.months.to.keep), ]
} else {
  df.cv <- data.frame(
    Type = character(), max.month = character(), MKT = numeric(), SE.MKT = numeric(),
    Factor = character(), Coef = numeric(), SE.Coef = numeric(), validation.error = numeric()
  )
}

df.cv.rank <- df.cv
df.cv.rank$key <- paste(df.cv.rank$Type, df.cv.rank$Factor, as.numeric(df.cv.rank$max.month))
df.cv.rank <- df.cv.rank[, c("key", "validation.error")]


# summarize ALL coefs with asymptotic inference ----
df.all <- df.f[(df.f$CV.key == "ALL"), ]

# Deduplicate: keep only the most recent entry per (Type, max.month, Factor)
# This handles cases where multiple cached files contain overlapping ALL entries
if (nrow(df.all) > 0 && "datetime" %in% colnames(df.all)) {
  df.all <- df.all[order(df.all$Type, df.all$max.month, df.all$Factor, df.all$datetime, decreasing = TRUE), ]
  df.all <- df.all[!duplicated(paste(df.all$Type, df.all$max.month, df.all$Factor)), ]
}
for (i in seq_len(nrow(df.all))) {
  factor <- as.character(df.all$Factor[i])
  # Only map Coef and SE.Coef if 'factor' is a specific factor name (e.g., "Alpha", "SMB") and not "ALL"
  # For "ALL" models, we rely on the specific factor columns directly.
  if (factor != "ALL" && factor %in% colnames(df.all)) {
    df.all[i, "Coef"] <- df.all[i, factor]
    df.all[i, "SE.Coef"] <- df.all[i, paste0("SE.", factor)]
    if (paste0("SE.", factor, ".indep") %in% colnames(df.all)) {
      df.all[i, "SE.Coef.indep"] <- df.all[i, paste0("SE.", factor, ".indep")]
    }
    if (paste0("SE.", factor, ".jk") %in% colnames(df.all)) {
      df.all[i, "SE.Coef.jk"] <- df.all[i, paste0("SE.", factor, ".jk")]
    }
    # effective number of units of the second coefficient (influence-based, r5)
    if (paste0("neff.", factor) %in% colnames(df.all)) {
      df.all[i, "neff.Coef"] <- df.all[i, paste0("neff.", factor)]
    }
  }
}

# Select columns defensively
all_select_cols <- c(
  "Type", "max.month", "MKT", "SE.MKT", "SE.MKT.indep", "SE.MKT.jk", "Factor"
)

# For single factor runs we need Coef and SE.Coef. For ALL runs we keep the individual factor columns.
if (any(df.all$Factor != "ALL")) {
  all_select_cols <- c(all_select_cols, "Coef", "SE.Coef", "SE.Coef.indep", "SE.Coef.jk")
}
if (any(df.all$Factor == "ALL")) {
  available_factors <- setdiff(intersect(c("Alpha", "SMB", "HML", "ME", "ROE", "IA", "EG"), colnames(df.all)), c("MKT"))
  for (f in available_factors) {
    all_select_cols <- c(all_select_cols, f, paste0("SE.", f), paste0("SE.", f, ".indep"), paste0("SE.", f, ".jk"))
  }
}

# keep the boundary flag (SEs not valid on a bound) and the provenance columns
# for the table headers (review fix: these were dropped here)
all_select_cols <- c(
  all_select_cols, "at.bound", RESULT.PROVENANCE.COLUMNS,
  # inference (decision 2026-11-17): bandwidth S, vintages V, b = S / V,
  # effective number of units, jackknife diagnostics, Wald p-values
  "hac.S", "hac.V", "hac.b", "neff.units", "neff.MKT", "neff.Coef", "jk.units", "jk.selfcheck",
  "jk.n.bound", "jk.n.switch", "Wald.stat.MKT_1", "Wald.q", "Wald.p.value.MKT_1", "Wald.p.value.MKT_1.chisq",
  # admissible domain (assumption A7, r5): active constraints and the lowest
  # monthly gross benchmark return at the estimate; inference method
  "active.constraints", "min.G", "min.G.date", "gross.floor", "inference", "fit.method"
)
all_select_cols <- intersect(all_select_cols, colnames(df.all))
if (nrow(df.all) > 0) {
  df.all <- df.all[, all_select_cols]
  df.all <- df.all[order(df.all$Type, df.all$Factor, as.numeric(df.all$max.month)), ]
} else {
  df.all <- df.all[, all_select_cols, drop = FALSE]
}

# print 2 LaTeX tables -----
spec <- paste0(strsplit(suffix, "_")[[1]], collapse = "-")
types2print <- c("PE", "BO", "GroBO", "VC", "RE", "PD", "INF", "NATRES")

# Define header for repeated pages
add.to.row <- list(pos = list(0), command = "\\hline ")

# Helper to suppress repeated values for cleaner tables and rename cols
# Helper to suppress repeated values for cleaner tables and rename cols
prepare_and_print_tables <- function(df, types, table_type, suffix, spec) {
  # Determine weighting scheme for caption
  weighting_text <- "Vintage-Year Portfolios"
  if (grepl("FW", suffix)) {
    weighting_text <- "Fund-Size Weighted Vintage-Year Portfolios"
  } else if (grepl("EW", suffix)) {
    weighting_text <- "Equal-Weighted Vintage-Year Portfolios"
  }

  for (type in types) {
    df.sub <- df[df$Type == type, ]
    if (nrow(df.sub) == 0) next

    # Sort by Factor then MCM (numeric)
    # Ensure max.month is treated as numeric for sorting
    df.sub <- df.sub[order(df.sub$Factor, as.numeric(df.sub$max.month)), ]

    # Remove Type column as it's now in the caption
    df.sub$Type <- NULL

    # Rename columns based on table type
    if (table_type == "AI") {
      on.bound <- boundary.mask(df.sub)
      for (col in intersect(c("MKT", "Coef"), names(df.sub))) {
        value <- df.sub[[col]]
        df.sub[[col]] <- ifelse(is.na(value), NA_character_, paste0(formatC(value, format = "f", digits = 3),
          ifelse(on.bound, "BOUNDARYMARKER", "")))
      }
      for (col in intersect(c("SE.MKT", "SE.MKT.indep", "SE.Coef", "SE.Coef.indep"), names(df.sub))) {
        df.sub[[col]][on.bound] <- NA_real_
      }
      # Asymptotic Inference: select available columns
      has_indep <- "SE.MKT.indep" %in% colnames(df.sub)
      if (has_indep) {
          ai_cols_req <- c("max.month", "MKT", "SE.MKT", "SE.MKT.indep", "Factor", "Coef", "SE.Coef", "SE.Coef.indep")
      } else {
          ai_cols_req <- c("max.month", "MKT", "SE.MKT", "Factor", "Coef", "SE.Coef")
      }
      for (col in ai_cols_req) {
          if (!col %in% colnames(df.sub)) df.sub[[col]] <- NA
      }
      df.sub <- df.sub[, ai_cols_req]

      cap <- paste0(
        "Asymptotic results for ", type, " Funds. ",
        "This table reports asymptotic parameter estimates for linear two-factor models using ", weighting_text, ". ",
        "Standard errors are estimated by a SHAC estimator (Equation \\ref{eq:asy_se}) with a Bartlett kernel over vintage years and truncation $S=\\lceil 1.3\\sqrt{V}\\rceil$ (Lazarus, Lewis, Stock and Watson, 2018), $V$ the number of vintages; tests use Kiefer--Vogelsang fixed-$b$ critical values with $b=S/V$. ",
        "MCM (Maximum Compounding Month) denotes the maximum cash flow compounding horizon. ",
        "$^{\\dagger}$ denotes an estimate on a parameter bound; its asymptotic standard errors are not valid and are not reported."
      )
      lab <- paste0("tab:ai_", suffix, "_", type)

      # Build header based on available columns
      if (has_indep) {
        header_row <- paste0(
          "\\hline\n",
          " & Market& Market & Market & Second & Second & Second & Second \\\\ \n",
          "\\textbf{MCM} & \\textbf{Loading} & \\textbf{SE} & \\textbf{SE Ind.} & \\textbf{Name} & \\textbf{Loading} & \\textbf{SE} & \\textbf{SE Ind.} \\\\ \n",
          "\\hline\n"
        )
      } else {
        header_row <- paste0(
          "\\hline\n",
          " & Market& Market & Second & Second & Second \\\\ \n",
          "\\textbf{MCM} & \\textbf{Loading} & \\textbf{SE} & \\textbf{Name} & \\textbf{Loading} & \\textbf{SE} \\\\ \n",
          "\\hline\n"
        )
      }
    } else {
      # Cross Validation
      # Desired columns: MCM, Market Loading, Market SE, Second Name, Second Loading, Second SE, CV-Error

      # Original cols: MKT, SE.MKT, Factor, Coef, SE.Coef, validation.error, max.month
      cv_table_cols <- c("max.month", "MKT", "SE.MKT", "Factor", "Coef", "SE.Coef", "validation.error")
      for (col in cv_table_cols) {
          if (!col %in% colnames(df.sub)) df.sub[[col]] <- NA
      }
      df.sub <- df.sub[, cv_table_cols]

      cap <- paste0(
        "Cross-validation results for ", type, " Funds. ",
        "The table reports the validation error and average parameter estimates from $hv$-block cross-validation. ",
        "The estimation uses ", weighting_text, ". ",
        "MCM (Maximum Compounding Month) denotes the maximum cash flow compounding horizon."
      )
      lab <- paste0("tab:cv_", suffix, "_", type)

      # Custom Header for CV
      # MCM | Market Loading | Market SE | Second Name | Second Loading | Second SE | CV-Error
      header_row <- paste0(
        "\\hline\n",
        " & Market & Market & Second & Second & Second &  \\\\ \n",
        "\\textbf{MCM} & \\textbf{Loading} & \\textbf{SE} & \\textbf{Name} & \\textbf{Loading} & \\textbf{SE} & \\textbf{CV-Error} \\\\ \n",
        "\\hline\n"
      )
    }

    # No longer suppressing repeated MCM values
    # df.sub$MCM <- as.character(df.sub$MCM) # Not needed if we keep original max.month

    # Define add.to.row list
    add.to.row <- list(pos = list(0), command = header_row)

    print(
      xtable::xtable(df.sub,
        caption = cap,
        label = lab, digits = 3
      ),
      include.rownames = FALSE,
      include.colnames = FALSE, # We provide our own headers
      floating = TRUE,
      # tabular.environment = "longtable", # we want to use normal tables!
      add.to.row = add.to.row,
      hline.after = c(nrow(df.sub)), # Only hline at bottom, top handled by header
      NA.string = "--",
      sanitize.text.function = function(x) gsub("BOUNDARYMARKER", "$^{\\dagger}$", xtable::sanitize(x), fixed = TRUE),
      file = outfile,
      append = TRUE
    )
  }
}

outfile <- paste0(dir.cache, "/empirical_tables.tex")
# Create cache directory if it doesn't exist (e.g., FF3 runs analyzed before estimation)
if (!dir.exists(dir.cache)) dir.create(dir.cache, recursive = TRUE)
# Clear file if it exists
writeLines(result.table.provenance(df.f), outfile)

prepare_and_print_tables(df.all, types2print, "AI", suffix, spec)
# cat("\\clearpage\n", file = outfile, append = TRUE)
prepare_and_print_tables(df.cv, types2print, "CV", suffix, spec)


# abs summary ----
sum.abs <- function(df, inference = "") {
  if ("SE.MKT.indep" %in% colnames(df)) {
    cols <- c("MKT", "SE.MKT", "SE.MKT.indep", "Coef", "SE.Coef", "SE.Coef.indep")
  } else {
    cols <- c("MKT", "SE.MKT", "Coef", "SE.Coef")
  }

  if (nrow(df) == 0) {
    # Return empty row with correct structure
    df.out <- data.frame(matrix(NA, nrow = 1, ncol = length(cols)))
    colnames(df.out) <- cols
    df.out$Weighting <- suffix
    df.out$Inference <- inference
    df.out <- df.out[, c("Weighting", "Inference", cols)]
    return(df.out)
  }

  df.out <- data.frame(apply(df[, intersect(cols, names(df)), drop = FALSE], 2, function(x) {
    sum(abs(as.numeric(x))) / nrow(df)
  }))
  colnames(df.out) <- "sum.abs"
  df.out <- data.frame(t(df.out))

  # Only select columns that exist in df.out
  cols <- intersect(cols, colnames(df.out))
  df.out <- df.out[, cols]
  df.out$Weighting <- suffix
  df.out$Inference <- inference
  df.out <- df.out[, c("Weighting", "Inference", cols)]
  return(df.out)
}

df.cv.abs <- sum.abs(df.cv, "cross-validation")
df.all.abs <- sum.abs(df.all, "asymptotic")

df.abs <- rbind.all.columns(df.all.abs, df.cv.abs)

print(xtable::xtable(df.abs,
  caption = "Sum of absolute values.",
  label = "tab:ai_sum_abs"
), include.rownames = FALSE)
## write csvs ----
write.csv(df.cv.abs, paste0(dir.cache, "/0_cross_validation_sumabs.csv"))
write.csv(df.all.abs, paste0(dir.cache, "/0_asymptotic_inference_sumabs.csv"))

# Write Cross-Validation
if (nrow(df.cv) > 0) {
  df.cv$t.MKT <- df.cv$MKT / df.cv$SE.MKT
  if ("Coef" %in% colnames(df.cv)) {
    df.cv$t.Coef <- df.cv$Coef / df.cv$SE.Coef
    df.cv$sig <- ((abs(df.cv$t.MKT) > 1.96) & (abs(df.cv$t.Coef) > 1.96))
    sum(df.cv$sig, na.rm = TRUE)
  }
  nrow(df.cv)

  df.cv <- df.cv[order(df.cv$Type, df.cv$max.month, df.cv$Factor), ]
}
write.csv(df.cv, paste0(dir.cache, "/0_cross_validation_summary.csv"))

# Write asymptotic inference
if (nrow(df.all) > 0) {
  df.all$t.MKT <- df.all$MKT / df.all$SE.MKT
  # two-sided p-values of t (H0: coefficient = 0) from the Kiefer-Vogelsang
  # fixed-b distribution at b = hac.b; normal p-values if b is not recorded
  # (results from before the decision of 17 November 2026)
  p.t <- function(t) if ("hac.b" %in% colnames(df.all)) fixedb.pvalue(t^2, df.all$hac.b, 1) else 2 * pnorm(-abs(t))
  # Wald test of beta = 1 (others 0): fixed-b p-value from the stored statistic
  # (r5 rows store Wald.stat.MKT_1, Wald.q and hac.b; older rows the p-value)
  if (all(c("Wald.stat.MKT_1", "Wald.q", "hac.b") %in% colnames(df.all))) {
    df.all$Wald.p.value.MKT_1 <- fixedb.pvalue(df.all$Wald.stat.MKT_1, df.all$hac.b, df.all$Wald.q)
  }
  df.all$p.MKT <- p.t(df.all$t.MKT)
  if ("Coef" %in% colnames(df.all)) {
    df.all$t.Coef <- df.all$Coef / df.all$SE.Coef
    df.all$p.Coef <- p.t(df.all$t.Coef)
    df.all$sig <- (df.all$p.MKT < 0.05) & (df.all$p.Coef < 0.05)
    sum(df.all$sig, na.rm = TRUE)
  }
  nrow(df.all)
  # Independent SE t-stats (only if .indep columns exist)
  if ("SE.MKT.indep" %in% colnames(df.all)) {
    df.all$t.MKT.indep <- df.all$MKT / df.all$SE.MKT.indep
    if ("Coef" %in% colnames(df.all) && "SE.Coef.indep" %in% colnames(df.all)) {
      df.all$t.Coef.indep <- df.all$Coef / df.all$SE.Coef.indep
      df.all$sig.indep <- ((abs(df.all$t.MKT.indep) > 1.96) & (abs(df.all$t.Coef.indep) > 1.96))
      sum(df.all$sig.indep)
    }
  }
  nrow(df.all)
  # View(df.all[df.all$Factor == "MKT", ])

  df.all <- df.all[order(df.all$Type, df.all$max.month, df.all$Factor), ]
}
write.csv(df.all, paste0(dir.cache, "/0_asymptotic_inference_summary.csv"))


# summarize both best ----
df.cv.best <- list()
for (Type in unique(df.cv$Type)) {
  for (max.month in unique(df.cv$max.month)) {
    df.ss <- df.cv[(df.cv$Type == Type) & (df.cv$max.month == max.month), ]
    df.ss <- df.ss[!is.na(df.ss$validation.error), ] # remove NA validation errors
    if (nrow(df.ss) == 0) next
    df.ss <- df.ss[df.ss$validation.error == min(df.ss$validation.error), ]
    if (nrow(df.ss) == 0) next
    # we just want "one best" model
    if (nrow(df.ss) > 1) {
      # df.ss <- df.ss[df.ss$MKT > 0, ] # hope that exactly one model remains
      if ("t.Coef" %in% colnames(df.ss) && !all(is.na(df.ss$t.Coef))) {
          df.ss <- df.ss[df.ss$t.Coef == max(df.ss$t.Coef, na.rm=TRUE), ] # choose most significant model
      } else {
          df.ss <- df.ss[1, , drop = FALSE]
      }
    }

    if (nrow(df.ss) > 0 && df.ss$Factor[1] == "Alpha") {
      df.ss <- df.cv[(df.cv$Type == Type) & (df.cv$max.month == max.month), ]
      df.ss <- df.ss[!is.na(df.ss$validation.error), ]
      if (nrow(df.ss) >= 2) {
        df.ss <- df.ss[df.ss$validation.error <= sort(df.ss$validation.error)[2], ]
      }
    }
    df.ss <- df.ss[, c("Type", "max.month", "Factor")]
    df.cv.best[[paste(Type, max.month)]] <- df.ss
  }
}
if (length(df.cv.best) > 0) {
  df.cv.best <- data.frame(do.call(rbind, df.cv.best), row.names = NULL)
  df.cv.best$best.key <- paste(df.cv.best$Type, df.cv.best$max.month, df.cv.best$Factor)

  df.best <- df.f[(df.f$CV.key == "ALL"), ]
  df.best$best.key <- paste(df.best$Type, df.best$max.month, df.best$Factor)
  df.best <- df.best[df.best$best.key %in% df.cv.best$best.key, ]
  df.best$Tables <- NA
  df.best$ID <- paste(df.best$weighting, df.best$max.month, sep = " - ")

  # Determine output columns based on available factor columns
  if ("Alpha" %in% colnames(df.best)) {
    # Alpha is a monthly rate; annualize with ^12 (not ^4)
    df.best$Alpha.p.a. <- (1 + df.best$Alpha)^12 - 1
  }
  if (grepl("^ff3_factors", prefix)) {
    # FF3 factors: MKT, SMB, HML (+ Alpha if present)
    cols <- c("ID", "Type", "MKT")
    if ("Alpha" %in% colnames(df.best)) cols <- c(cols, "Alpha.p.a.")
    cols <- c(cols, intersect(c("SMB", "HML"), colnames(df.best)), "Tables")
  } else if ("ME" %in% colnames(df.best)) {
    # q-factor model
    cols <- c("ID", "Type", "MKT", "Alpha.p.a.", "ME", "IA", "ROE", "EG", "Tables")
  } else {
    # MSCI factors
    df.best$Alpha.p.a. <- NA
    cols <- c("ID", "Type", "MKT", "SMB", "HML", "HDY", "QLT", "MOM", "LOV", "Tables")
  }
  # Ensure all cols exist in df.best (add NA for missing columns)
  for (col in cols) {
    if (!col %in% colnames(df.best)) df.best[[col]] <- NA
  }

  if (nrow(df.best) > 0) {
    print(xtable::xtable(df.best[, cols],
      caption = "SDF models with best out-of-sample performance as determined by $hv$-block cross-validation. ID consists of cash flow weighting - maximum quarter",
      label = "tab:result_summary"
    ), include.rownames = FALSE)
  }

  write.csv(df.best, paste0(dir.cache, "/0_best_models_summary.csv"))
} else {
  # No CV results — create empty df.best for downstream compatibility
  df.best <- data.frame()
  message("No cross-validation results available for best model selection.")
}


# plot log returns ----
plot.log.return <- function(type, df.f) {
  df.f <- df.f[(df.f$Type == type), ]
  df.f <- df.f[(df.f$CV.key == "ALL"), ]
  # df.f <- df.f[df.f$Factor != "Alpha", ]

  # Dynamically determine factor columns present in both df.f and df.q
  all_possible_factors <- c("RF", "MKT", "ME", "IA", "ROE", "EG", "SMB", "HML", "Alpha")
  cols <- intersect(all_possible_factors, intersect(colnames(df.f), colnames(df.q)))

  if (nrow(df.f) == 0 || length(cols) == 0) {
    message("No data for plot.log.return: ", type)
    return(invisible(NULL))
  }

  # Alpha is already a monthly rate matching df.q's monthly frequency (no quarterly-to-monthly conversion)
  df.f[is.na(df.f)] <- 0
  # Use drop = FALSE to ensure matrix result even with single row
  df <- data.frame((as.matrix(df.q[, cols])) %*% t(as.matrix(
    df.f[, cols, drop = FALSE]
  )))
  colnames(df) <- paste(df.f$Factor, df.f$max.month, df.f$weighting)
  cols <- colnames(df)

  df$Date <- df.q$Date
  df <- df[df$Date > as.Date("1989-12-31"), ]

  # plot
  plot(
    x = df$Date, y = log(cumprod(1 + df[, 1])), main = type, type = "l",
    xlab = "Date", ylab = "Cum. Log Return", ylim = c(0, 15)
  )
  color <- 0
  i <- 1
  for (col in cols) {
    color <- color + 1
    lines(df$Date, log(cumprod(1 + df[, col])), col = color)
  }
  legend("topleft", bty = "n", legend = cols, col = 1:length(cols), lty = 1)
}

if (FALSE) {
  df.best1 <- df.best2 <- NULL
  df.best1 <- read.csv("data_out/cache_q_factors_EW_VYP_SL/0_best_models_summary.csv")
  df.best2 <- read.csv("data_out/cache_q_factors_FW_VYP_SL/0_best_models_summary.csv")
  # df.best1 <- read.csv("data_out/cache_q_factors_EW/0_best_models_summary.csv")
  # df.best2 <- read.csv("data_out/cache_q_factors_FW/0_best_models_summary.csv")
  df.best <- rbind(df.best1, df.best2)
}


if (nrow(df.best) > 0) {
  plot.log.return("VC", df.best)

  # setEPS()
  # postscript(paste0(dir.cache, "/", "0_SDF_realizations.eps"),
  #           width = 6, height = 7.5, family = "Helvetica", pointsize = 5)
  # par(mfrow=c(3,2), cex=1.2, mar = c(4.1, 4.1, 3.1, 1.1))

  for (type in c("PE", "VC", "PD", "RE", "NATRES", "INF")) {
    plot.log.return(type, df.best)
    abline(h = c(0), col = "grey", lty = 3)
  }
}

# par(mfrow=c(1,1), cex=1)
# dev.off()

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
  data.out.folder <- "data_out_2026-emp-B"
}

list.cache <- list()
dir.cache <- paste0(data.out.folder, "/cache_", prefix, suffix)

for (file in list.files(dir.cache)) {
  if (substr(file, 1, 1) == 0) next
  if (!grepl("\\.csv$", file)) next
  df.f <- read.csv(paste0(dir.cache, "/", file))
  df.f$CV.key <- as.character(df.f$CV.key)
  list.cache[[file]] <- df.f
}
df.f <- data.frame(do.call(rbind, list.cache))
rownames(df.f) <- NULL
df.f$X <- NULL
df.f$RF <- 1
# df.f[is.na(df.f)] <- 0
# df.f$validation.error[is.na(df.f$validation.error)] <- Inf

current_dir <- getwd()
parent_dir <- dirname(current_dir)
file_path <- file.path(parent_dir, "empirical", "data_prepared", "q_factors.csv")

df.q <- read.csv(file_path)
df.q$Date <- as.Date(df.q$Date)
df.q$Alpha <- 1

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
  df <- df.f[
    (df.f$CV.key != "ALL") & (df.f$Factor == Factor),
    c("MKT", Factor, "Type", "max.month", "validation.error")
  ]
  ddf <- df[, !is.na(df[1, ])]
  df$id <- paste0(df$Type, "_", df$max.month)
  df$Type <- NULL
  df$max.month <- NULL

  if (is.data.frame(df)) {
    m <- aggregate(. ~ id, mean, data = df)
    s <- aggregate(. ~ id, sd, data = df)
    s["validation.error"] <- NULL
  }
  n <- ncol(s)
  a <- 2
  b <- 3
  if (n == b) {
    names(m)[b] <- "Coef"
    names(s)[b] <- "Coef"
  }

  names(m)[a:n] <- paste0(names(m)[a:n])
  names(s)[a:n] <- paste0("SE.", names(s)[a:n])

  df <- merge(m, s, by = "id")
  df$Factor <- Factor
  cv.res[[Factor]] <- df
}

df.cv <- data.frame(Reduce(rbind.all.columns, cv.res))
for (i in 1:nrow(df.cv)) {
  df.cv[i, "Type"] <- strsplit(df.cv$id, "_")[[i]][1]
  df.cv[i, "max.month"] <- strsplit(df.cv$id, "_")[[i]][2]
}
df.cv <- df.cv[, c("Type", "max.month", "MKT", "SE.MKT", "Factor", "Coef", "SE.Coef", "validation.error")]
df.cv <- df.cv[order(df.cv$Type, df.cv$Factor, as.numeric(df.cv$max.month)), ]

max.months.to.keep <- c("1", "30", "60", "120", "180")
df.cv <- df.cv[df.cv$max.month %in% max.months.to.keep, ]

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
for (i in 1:nrow(df.all)) {
  factor <- as.character(df.all$Factor[i])
  df.all[i, "Coef"] <- df.all[i, factor]
  df.all[i, "SE.Coef"] <- df.all[i, paste0("SE.", factor)]
  df.all[i, "SE.Coef.indep"] <- df.all[i, paste0("SE.", factor, ".indep")]
}
df.all <- df.all[, c(
  "Type", "max.month", "MKT", "SE.MKT", "SE.MKT.indep",
  "Factor", "Coef", "SE.Coef", "SE.Coef.indep"
)]

df.all <- df.all[order(df.all$Type, df.all$Factor, as.numeric(df.all$max.month)), ]

# print 2 LaTeX tables -----
spec <- paste0(strsplit(suffix, "_")[[1]], collapse = "-")
types2print <- c("PE")

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

    # Sort by Factor then MCM (numeric)
    # Ensure max.month is treated as numeric for sorting
    df.sub <- df.sub[order(df.sub$Factor, as.numeric(df.sub$max.month)), ]

    # Remove Type column as it's now in the caption
    df.sub$Type <- NULL

    # Rename columns based on table type
    if (table_type == "AI") {
      # Asymptotic Inference
      # Desired columns: MCM, Market Loading, Market SE, Market SE Ind., Second Name, Second Loading, Second SE, Second SE Ind.

      # Rename for internal consistency before printing (though headers will be custom)
      # We need correct column order for the data rows

      # Map original names to temp names to ensure correct ordering
      # Original cols: MKT, SE.MKT, SE.MKT.indep, Factor, Coef, SE.Coef, SE.Coef.indep, max.month

      df.sub <- df.sub[, c("max.month", "MKT", "SE.MKT", "SE.MKT.indep", "Factor", "Coef", "SE.Coef", "SE.Coef.indep")]

      cap <- paste0(
        "Asymptotic results for ", type, " Funds. ",
        "This table reports asymptotic parameter estimates for linear two-factor models using ", weighting_text, ". ",
        "Standard errors are estimated by a SHAC estimator (Equation \\ref{eq:asy_se}) with a Bartlett kernel bandwidth of $D=12$ vintage years to correct for dependence between overlapping funds. ",
        "MCM (Maximum Compounding Month) denotes the maximum cash flow compounding horizon."
      )
      lab <- paste0("tab:ai_", suffix, "_", type)

      # Custom Header for AI
      # MCM | Market Loading | Market SE | Market SE Ind. | Second Name | Second Loading | Second SE | Second SE Ind.
      header_row <- paste0(
        "\\hline\n",
        " & Market& Market & Market & Second & Second & Second & Second \\\\ \n",
        "\\textbf{MCM} & \\textbf{Loading} & \\textbf{SE} & \\textbf{SE Ind.} & \\textbf{Name} & \\textbf{Loading} & \\textbf{SE} & \\textbf{SE Ind.} \\\\ \n",
        "\\hline\n"
      )
    } else {
      # Cross Validation
      # Desired columns: MCM, Market Loading, Market SE, Second Name, Second Loading, Second SE, CV-Error

      # Original cols: MKT, SE.MKT, Factor, Coef, SE.Coef, validation.error, max.month
      df.sub <- df.sub[, c("max.month", "MKT", "SE.MKT", "Factor", "Coef", "SE.Coef", "validation.error")]

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
      file = outfile,
      append = TRUE
    )
  }
}

outfile <- paste0(dir.cache, "/empirical_tables.tex")
# Clear file if it exists
cat("", file = outfile)

prepare_and_print_tables(df.all, types2print, "AI", suffix, spec)
cat("\\clearpage\n", file = outfile, append = TRUE)
prepare_and_print_tables(df.cv, types2print, "CV", suffix, spec)


# abs summary ----
sum.abs <- function(df, inference = "") {
  df.out <- data.frame(apply(df, 2, function(x) {
    sum(abs(as.numeric(x))) / nrow(df)
  }))
  colnames(df.out) <- "sum.abs"
  df.out <- data.frame(t(df.out))

  if ("SE.MKT.indep" %in% colnames(df)) {
    cols <- c("MKT", "SE.MKT", "SE.MKT.indep", "Coef", "SE.Coef", "SE.Coef.indep")
  } else {
    cols <- c("MKT", "SE.MKT", "Coef", "SE.Coef")
  }
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

# Write Corss-Validation
df.cv$t.MKT <- df.cv$MKT / df.cv$SE.MKT
df.cv$t.Coef <- df.cv$Coef / df.cv$SE.Coef
df.cv$sig <- ((abs(df.cv$t.MKT) > 1.96) & (abs(df.cv$t.Coef) > 1.96))
sum(df.cv$sig)
nrow(df.cv)

df.cv <- df.cv[order(df.cv$Type, df.cv$max.month, df.cv$Factor), ]
write.csv(df.cv, paste0(dir.cache, "/0_cross_validation_summary.csv"))

# Write asymptotic inference
df.all$t.MKT <- df.all$MKT / df.all$SE.MKT
df.all$t.MKT.indep <- df.all$MKT / df.all$SE.MKT.indep
df.all$t.Coef <- df.all$MKT / df.all$SE.Coef
df.all$t.Coef.indep <- df.all$MKT / df.all$SE.Coef.indep
df.all$sig <- ((abs(df.all$t.MKT) > 1.96) & (abs(df.all$t.Coef) > 1.96))
sum(df.all$sig)
nrow(df.all)
df.all$sig.indep <- ((abs(df.all$t.MKT.indep) > 1.96) & (abs(df.all$t.Coef.indep) > 1.96))
sum(df.all$sig.indep)
nrow(df.all)
# View(df.all[df.all$Factor == "MKT", ])

df.all <- df.all[order(df.all$Type, df.all$max.month, df.all$Factor), ]
write.csv(df.all, paste0(dir.cache, "/0_asymptotic_inference_summary.csv"))


# summarize both best ----
df.cv.best <- list()
for (Type in unique(df.cv$Type)) {
  for (max.month in unique(df.cv$max.month)) {
    df.ss <- df.cv[(df.cv$Type == Type) & (df.cv$max.month == max.month), ]
    df.ss <- df.ss[df.ss$validation.error == min(df.ss$validation.error), ]
    # we just want "one best" model
    if (nrow(df.ss) > 1) {
      # df.ss <- df.ss[df.ss$MKT > 0, ] # hope that exactly one model remains
      df.ss <- df.ss[df.ss$t.Coef == max(df.ss$t.Coef), ] # choose most significant model
    }

    if (df.ss$Factor == "Alpha") {
      df.ss <- df.cv[(df.cv$Type == Type) & (df.cv$max.month == max.month), ]
      df.ss <- df.ss[df.ss$validation.error <= sort(df.ss$validation.error)[2], ]
    }
    df.ss <- df.ss[, c("Type", "max.month", "Factor")]
    df.cv.best[[paste(Type, max.month)]] <- df.ss
  }
}
df.cv.best <- data.frame(do.call(rbind, df.cv.best), row.names = NULL)
df.cv.best$best.key <- paste(df.cv.best$Type, df.cv.best$max.month, df.cv.best$Factor)


df.best <- df.f[(df.f$CV.key == "ALL"), ]
df.best$best.key <- paste(df.best$Type, df.best$max.month, df.best$Factor)
df.best <- df.best[df.best$best.key %in% df.cv.best$best.key, ]
df.best$Tables <- NA
df.best$ID <- paste(df.best$weighting, df.best$max.month, sep = " - ")

# Alpha just included in q-factor models, not MSCI models (so far)
if ("Alpha" %in% colnames(df.best)) {
  df.best$Alpha.p.a. <- (1 + df.best$Alpha)^4 - 1
  cols <- c("ID", "Type", "MKT", "Alpha.p.a.", "ME", "IA", "ROE", "EG", "Tables") # q factors
} else {
  # assuming MSCI factors
  df.best$Alpha.p.a. <- NA
  cols <- c("ID", "Type", "MKT", "SMB", "HML", "HDY", "QLT", "MOM", "LOV", "Tables") # msci factors
}

print(xtable::xtable(df.best[, cols],
  caption = "SDF models with best out-of-sample performance as determined by $hv$-block cross-validation. ID consists of cash flow weighting - maximum quarter",
  label = "tab:result_summary"
), include.rownames = FALSE)


write.csv(df.best, paste0(dir.cache, "/0_best_models_summary.csv"))


# plot log returns ----
plot.log.return <- function(type, df.f) {
  df.f <- df.f[(df.f$Type == type), ]
  df.f <- df.f[(df.f$CV.key == "ALL"), ]
  # df.f <- df.f[df.f$Factor != "Alpha", ]

  cols <- c("RF", "MKT", "ME", "IA", "ROE", "EG", "Alpha")
  cols <- cols[cols %in% colnames(df.f)]

  if ("Alpha" %in% cols) {
    df.f["Alpha"] <- (1 + df.f["Alpha"])^(1 / 3) - 1 # quarterly to monthly
  }

  df.f[is.na(df.f)] <- 0
  df <- data.frame((as.matrix(df.q[, cols])) %*% t(as.matrix(
    df.f[, cols]
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


plot.log.return("VC", df.best)


# setEPS()
# postscript(paste0(dir.cache, "/", "0_SDF_realizations.eps"),
#           width = 6, height = 7.5, family = "Helvetica", pointsize = 5)
# par(mfrow=c(3,2), cex=1.2, mar = c(4.1, 4.1, 3.1, 1.1))

for (type in c("PE", "VC", "PD", "RE", "NATRES", "INF")) {
  plot.log.return(type, df.best)
  abline(h = c(0), col = "grey", lty = 3)
}

# par(mfrow=c(1,1), cex=1)
# dev.off()

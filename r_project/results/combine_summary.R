# combine summary ----
if (sys.nframe() == 0L) rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


data.out.folder <- "data_out_2026-emp"
prefix <- "q_factors_"
suffixes <- c("EW", "EW_VYP", "FW", "FW_VYP")
suffixes <- c("preqin_EW_VYP", "preqin_FW_VYP")

files <- c("asymptotic_inference", "cross_validation")

list.cache <- list()
for (suffix in suffixes) {
  for (file in files) {
    path <- paste0(data.out.folder, "/cache_", prefix, suffix, "/0_", file, "_summary.csv")
    df <- read.csv(path)
    df$X <- NULL
    # Handle different column names for time horizon:
    # - New format: max.month (from updated analyze_result.R)
    # - Legacy format: max.quarter or Max.Quarter
    if ("max.month" %in% colnames(df)) {
      df$MQ <- df$max.month
      df$max.month <- NULL
    } else if ("Max.Quarter" %in% colnames(df)) {
      df$MQ <- df$Max.Quarter
      df$Max.Quarter <- NULL
    } else if ("max.quarter" %in% colnames(df)) {
      df$MQ <- df$max.quarter
      df$max.quarter <- NULL
    }
    df$PF <- ifelse(grepl("VYP", suffix), "vy", "sf")
    df$WE <- ifelse(grepl("EW", suffix), "eq", "fs")
    df$IN <- ifelse(grepl("asymptotic_inference", file), "ai", "cv")

    list.cache[[paste0(suffix, file)]] <- df
  }
}

rbind.all.columns <- function(x, y) {
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))

  x[, c(as.character(y.diff))] <- NA

  y[, c(as.character(x.diff))] <- NA

  return(rbind(x, y))
}
df <- data.frame(Reduce(rbind.all.columns, list.cache))
colnames(df)
col.order <- c(
  "IN", "PF", "MQ", "WE", "Type",
  "MKT", "SE.MKT", "SE.MKT.indep",
  "Factor", "Coef", "SE.Coef", "SE.Coef.indep",
  "validation.error"
)
df <- df[order(df$IN, df$PF, df$MQ, df$WE, df$Type, df$Factor), col.order]

fts <- c("PE", "VC")
for (fund.type in fts) {
  df1 <- df[(df$Type == fund.type) & (df$Factor == "MKT"), ]
  df1[, c("Factor", "Coef", "SE.Coef", "SE.Coef.indep", "validation.error")] <- NULL

  print(xtable::xtable(df1,
    caption = paste0("MKT-factor models for fund type ", fund.type, "."),
    label = paste0("tab:result_", fund.type)
  ), include.rownames = FALSE)
}

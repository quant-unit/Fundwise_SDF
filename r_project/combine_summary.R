# combine summary ----
if(sys.nframe() == 0L) rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


list.cache <- list()
prefix <- "q_factors_"
suffixes <- c("EW", "EW_VYP", "FW", "FW_VYP")
files <- c("asymptotic_inference", "cross_validation")
for(suffix in suffixes) {
  for(file in files) {
    path <- paste0("data_out/cache_", prefix, suffix, "/0_", file, "_summary.csv")
    df <- read.csv(path)
    df$X <- NULL
    if("Max.Quarter" %in% colnames(df)) {
      df$max.quarter <- df$Max.Quarter
      df$Max.Quarter <- NULL
    }
    df$MQ <- df$max.quarter
    df$max.quarter <- NULL
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
col.order <- c("IN", "PF", "MQ", "WE", "Type", 
               "MKT", "SE.MKT", "SE.MKT.indep", 
               "Factor", "Coef", "SE.Coef", "SE.Coef.indep", 
               "validation.error")
df <- df[order(df$IN, df$PF, df$MQ, df$WE, df$Type, df$Factor), col.order]

fts <- c("PE", "VC")
for (fund.type in fts) {
  df1 <- df[(df$Type == fund.type) & (df$Factor == "MKT"), ]
  df1[, c("Factor", "Coef", "SE.Coef", "SE.Coef.indep", "validation.error")] <- NULL
  
  print(xtable::xtable(df1, 
                       caption = paste0("MKT-factor models for fund type ", fund.type, "."),
                       label = paste0("tab:result_", fund.type)), include.rownames=FALSE)
}


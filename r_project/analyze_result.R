# analyze results -----
if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

list.cache <- list()
dir.cache <- "data_out/cache_q_factors_EW_VYP"
for(file in list.files(dir.cache)) {
  if(substr(file,1,1) == 0) next
  df.f <- read.csv(paste0(dir.cache, "/", file))
  df.f$CV.key <- as.character(df.f$CV.key)
  list.cache[[file]] <- df.f
}
df.f <- data.frame(do.call(rbind, list.cache))
rownames(df.f) <- NULL
df.f$X <- NULL
df.f$RF <- 1
#df.f[is.na(df.f)] <- 0
#df.f$validation.error[is.na(df.f$validation.error)] <- Inf

df.q <- read.csv("data_prepared/q_factors.csv")
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
q.factors <- levels(df.f$Factor)

cv.res <- list()
for(Factor in q.factors) {
  df <- df.f[(df.f$CV.key != "ALL") & (df.f$Factor == Factor), 
             c("MKT", Factor, "Type", "max.quarter", "validation.error")]
  ddf <- df[, !is.na(df[1, ])]
  df$id <- paste0(df$Type, "_", df$max.quarter)
  df$Type <- NULL
  df$max.quarter <- NULL
  
  if( is.data.frame(df) ) {
    m <-   aggregate(. ~ id, mean, data = df)
    s <-   aggregate(. ~ id, sd, data = df)
    s["validation.error"] <- NULL
  }
  n <- ncol(s)
  a <- 2
  b <- 3
  if(n == b) {
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
df.cv$validation.error <- as.character(round(df.cv$validation.error))
for(i in 1:nrow(df.cv)) {
  df.cv[i, "Type"] <- strsplit(df.cv$id, "_")[[i]][1]
  df.cv[i, "Max.Quarter"] <- strsplit(df.cv$id, "_")[[i]][2]
}
df.cv <- df.cv[, c("Type", "Max.Quarter", "MKT", "SE.MKT", "Factor", "Coef", "SE.Coef", "validation.error")]
df.order <- data.frame(Type = c("PE", "VC", "PD", "RE", "NATRES", "INF"), Order = seq(1,6))
df.cv <- merge(df.cv, df.order, by ="Type")
df.cv <- df.cv[order(df.cv$Order, df.cv$Max.Quarter), ]
df.cv$Order <- NULL
df.cv.40 <- df.cv[df.cv$Max.Quarter == 40, ]
df.cv.60 <- df.cv[df.cv$Max.Quarter == 60, ]
df.cv.40$Max.Quarter <- df.cv.60$Max.Quarter <- NULL

write.csv(df.cv, paste0(dir.cache,"/0_cross_validation_summary.csv"))




# summarize all coefs ----
df.all <- df.f[(df.f$CV.key == "ALL"), ]
for(factor in q.factors) {
  for(i in 1:nrow(df.all)) {
    if(!is.na(df.all[i, factor])) {
      df.all[i, "Coef"] <- df.all[i, factor]
      df.all[i, "SE.Coef"] <- df.all[i, paste0("SE.", factor)]
      df.all[i, "SE.Coef.indep"] <- df.all[i, paste0("SE.", factor, ".indep")] 
    }
  }
}
df.all <- df.all[, c("Type", "max.quarter", "MKT", "SE.MKT", "SE.MKT.indep", 
                     "Factor", "Coef", "SE.Coef", "SE.Coef.indep")]

df.all40 <- df.all[df.all$max.quarter == 40, ]
df.all60 <- df.all[df.all$max.quarter == 60, ]
df.all40$max.quarter <- df.all60$max.quarter <- NULL




print(xtable::xtable(df.all40, 
                     caption = "Asymptotic inference with XXX weighting, max quarter 40, and $D=12$.",
                     label = "tab:ai_40"), include.rownames=FALSE)
print(xtable::xtable(df.cv.40, 
                     caption = "$hv$-block cross-validation with XXX weighting and max quarter 40.",
                     label = "tab:cv_40"), include.rownames=FALSE)
print(xtable::xtable(df.all60, 
                     caption = "Asymptotic inference with XXX weighting, max quarter 60, and $D=12$.",
                     label = "tab:ai_60"), include.rownames=FALSE)
print(xtable::xtable(df.cv.60, 
                     caption = "$hv$-block cross-validation with XXX weighting and max quarter 60.",
                     label = "tab:cv_60"), include.rownames=FALSE)


write.csv(df.all, paste0(dir.cache,"/0_asymptotic_inference_summary.csv"))


# summarize both best ----
df.cv$validation.error <- as.numeric(df.cv$validation.error)
df.cv.best <- list()
for(Type in df.cv$Type) {
  for(Max.Quarter in df.cv$Max.Quarter) {
    df.ss <- df.cv[(df.cv$Type == Type) & (df.cv$Max.Quarter == Max.Quarter), ]
    df.ss <- df.ss[df.ss$validation.error == min(df.ss$validation.error), ]
    df.ss <- df.ss[, c("Type", "Max.Quarter", "Factor")]
    df.cv.best[[paste(Type, Max.Quarter)]] <- df.ss
  }
}
df.cv.best <- data.frame(do.call(rbind, df.cv.best), row.names = NULL)
df.cv.best$best.key <- paste(df.cv.best$Type, df.cv.best$Max.Quarter, df.cv.best$Factor)


df.best <- df.f[(df.f$CV.key == "ALL"), ]
df.best$best.key <- paste(df.best$Type, df.best$max.quarter, df.best$Factor)
df.best <- df.best[df.best$best.key %in% df.cv.best$best.key, ]
df.best$Tables <- NA
df.best$Alpha.p.a. <- (1 + df.best$Alpha)^4 - 1
df.best$ID <- paste(df.best$weighting, df.best$max.quarter, sep = " - ")
cols <- c("ID", "Type", "MKT", "Alpha.p.a.", "ME", "IA", "ROE", "EG", "Tables")

print(xtable::xtable(df.best[, cols], 
                     caption = "SDF models with best out-of-sample performance as determined by $hv$-block cross-validation. ID consists of cash flow weighting - maximum quarter",
                     label = "tab:result_summary"), include.rownames=FALSE)

write.csv(df.best, paste0(dir.cache,"/0_best_models_summary.csv"))


# plot log returns ----
plot.log.return <- function(type, df.f) {
  df.f <- df.f[(df.f$Type == type), ]
  df.f <- df.f[(df.f$CV.key == "ALL"), ]
  
  cols <- c("RF", "MKT", "ME", "IA", "ROE", "EG", "Alpha")
  cols <- cols[cols %in% colnames(df.f)]
  
  if("Alpha" %in% cols) {
    df.f["Alpha"] <- (1 + df.f["Alpha"])^(1/3) - 1 # quarterly to monthly
  }
  
  df.f[is.na(df.f)] <- 0
  df <- data.frame((as.matrix(df.q[, cols])) %*% t(as.matrix(
    df.f[, cols]
  )))
  colnames(df) <- paste(df.f$Factor, df.f$max.quarter, df.f$weighting)
  cols <- colnames(df)
  
  df$Date <- df.q$Date
  df <- df[df$Date > as.Date("1989-12-31"), ]
  
  # plot
  plot(x = df$Date, y = log(cumprod(1 + df[, 1])), main = type, type = "l",
       xlab = "Date", ylab = "Cum. Log Return", ylim = c(0, 15))
  color <- 0 ; i <- 1
  for(col in cols) {
    color <- color + 1
    lines(df$Date, log(cumprod(1 + df[, col])), col = color)
  }
  legend("topleft", bty="n", legend = cols, col = 1:length(cols), lty = 1)
  
}

df.best1 <- df.best2 <- NULL
df.best1 <- read.csv("data_out/cache_q_factors_EW_VYP/0_best_models_summary.csv")
df.best2 <- read.csv("data_out/cache_q_factors_FW_VYP/0_best_models_summary.csv")
df.best <- rbind(df.best1, df.best2)

plot.log.return("VC", df.best)


setEPS()
postscript(paste0(dir.cache, "/", "0_SDF_realizations.eps"), 
           width = 6, height = 4, family = "Helvetica", pointsize = 5)
par(mfrow=c(2,2), cex=1.2, mar = c(4.1, 4.1, 3.1, 1.1))

for(type in c("PE", "VC", "PD", "RE")) {
  plot.log.return(type, df.best)
  abline(h=c(0), col = "grey", lty = 3)
}

par(mfrow=c(1,1), cex=1)
dev.off() 

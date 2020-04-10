# analyze results -----
if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

list.cache <- list()
dir.cache <- "data_out/cache_q_factors"
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


rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}


# summarize coefficient estimates -----
require(xtable)
q.factors <- c("MKT", "ME", "IA", "ROE", "EG")

cv.res <- list()
for(Factor in q.factors) {
  df <- df.f[(df.f$CV.key != "ALL") & (df.f$Factor == Factor), 
             c(q.factors, "Type", "max.quarter", "validation.error")]
  df <- df[, !is.na(df[1, ])]
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
df.cv <- df.cv[, c(8, 9, 2, 4:7, 3)]
df.cv <- df.cv[order(df.cv$Type, df.cv$Max.Quarter), ]
df.cv.40 <- df.cv[df.cv$Max.Quarter == 40, ]
df.cv.60 <- df.cv[df.cv$Max.Quarter == 60, ]
df.cv.40$Max.Quarter <- df.cv.60$Max.Quarter <- NULL

write.csv(df.cv, paste0(dir.cache,"/0_cross_validation_summary.csv"))

df.all <- df.f[(df.f$CV.key == "ALL"), ]
for(factor in q.factors[2:5]) {
  for(i in 1:nrow(df.all)) {
    if(!is.na(df.all[i, factor])) {
      df.all[i, "Coef"] <- df.all[i, factor]
      df.all[i, "SE.Coef"] <- df.all[i, paste0("SE.", factor)] 
    }
  }
}
df.all <- df.all[, c("Type", "max.quarter", "MKT", "SE.MKT", "Factor", "Coef", "SE.Coef", "Wald.p.value.MKT_1")]

df.all40 <- df.all[df.all$max.quarter == 40, ]
df.all60 <- df.all[df.all$max.quarter == 60, ]
df.all40$max.quarter <- df.all60$max.quarter <- NULL





print(xtable::xtable(df.all40, caption = "Asymptotic inference with max quarter 40. The Wald test null hypothesis is a MKT coefficient of one, and zeros for all other coefficients."), include.rownames=FALSE)
print(xtable::xtable(df.cv.40, caption = "$hv$-block cross-validation with max quarter 40."), include.rownames=FALSE)
print(xtable::xtable(df.all60, caption = "Asymptotic inference with max quarter 60. The Wald test null hypothesis is a MKT coefficient of one, and zeros for all other coefficients."), include.rownames=FALSE)
print(xtable::xtable(df.cv.60, caption = "$hv$-block cross-validation with max quarter 60."), include.rownames=FALSE)


write.csv(df.all, paste0(dir.cache,"/0_asymptotic_inference_summary.csv"))


# plot log returns ----
plot.log.return <- function(type, df.f) {
  df.f <- df.f[(df.f$Type == type), ]
  df.f <- df.f[(df.f$CV.key == "ALL"), ]
  
  cols <- c("RF", "MKT", "ME", "IA", "ROE", "EG")
  df.f[is.na(df.f)] <- 0
  df <- data.frame((as.matrix(df.q[, cols])) %*% t(as.matrix(
    df.f[, cols]
  )))
  colnames(df) <- df.f$Factor
  
  df$Date <- df.q$Date
  df <- df[df$Date > as.Date("1989-12-31"), ]
  
  # plot
  plot(x = df$Date, y = log(cumprod(1 + df[, 1])), main = type, type = "l",
       xlab = "Date", ylab = "Cum. Log Return", ylim = c(-5, 15))
  color <- 0 ; i <- 1
  for(col in cols[2:6]) {
    df.ss <- df[, colnames(df) == col]
    color <- color + 1
    for(i in 1:ncol(df.ss)) {
      lines(df$Date, log(cumprod(1 + df.ss[, i])), col = color)
    }
  }
  legend("topleft", bty="n", legend = cols[2:6], col = 1:5, lty = 1)
  
}

setEPS()
postscript(paste0(dir.cache, "/", "0_SDF_realizations.eps"), 
           width = 5.5, height = 4, family = "Helvetica", pointsize = 5)
par(mfrow=c(2,2), cex=1.0, mar = c(4.1, 4.1, 3.1, 1.1))

for(type in c("PE", "VC", "PD", "RE")) {
  plot.log.return(type, df.f)
  abline(h=c(0,5), col = "grey", lty = 3)
}

par(mfrow=c(1,1), cex=1)
dev.off() 

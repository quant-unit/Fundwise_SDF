# analyze balanced ensembles ----
rm(list=ls()) # remove workspace objects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

if (TRUE) {
  path <- "data/cache_msci_market_factors"
  l <- list()
  for(file in list.files(path)) {
    df <- read.csv(paste0(path, "/",file))
    if("xtime" %in% colnames(df)) {
      names(df)[names(df) == 'xtime'] <- 'xtimes'
    }
    l[[file]] <- df
  }
  df <- data.frame(do.call(rbind, l))
  df[is.na(df)] <- 0
  rm(path, l)
} else {
  df <- read.csv("data/Benchmark_SDF_MSCI_coefs.csv")
  
}
df <- df[!(df$Factor %in% c("ESG", "LOV", "MOM")), ]
df <- df[df$lambda == 0, ]
nrow(df) / length(levels(df$Type))
5 * 2 * 2 * 4
table(df$max.month)

factors <- c("MKT", "HML", "SMB", "HDY", "QLT") # , "MOM", "ESG", "LOV)
mean.sd.4 <- function(x) {paste0(round(mean(x[x!=0])/4, 2), " (", round(sd(x[x!=0])/4, 2), ")")}
mean.sd <- function(x) {paste0(round(mean(x[x!=0]), 2), " (", round(sd(x[x!=0]), 2), ")")}
df.ms <- aggregate(. ~ Type, df, mean.sd.4)[, c("Type", factors[-1])]
df.mkt <- aggregate(. ~ Type, df, mean.sd)[, c("Type", "MKT")]
df.ms <- merge(df.mkt, df.ms, by = "Type")
df.ms$Type <- as.character(df.ms$Type)
df.ms <- rbind(df.ms, c("MKT" ,"1", "0", "0", "0", "0"))
df <- aggregate(. ~ Type, df, mean)[, c("Type", factors)]
#View(df.ms)

print(xtable::xtable(df.ms, 
                     caption = "Average coefficients of MSCI factor models with standard deviation in parentheses.",
                     label = "tab:average_coefs"), include.rownames=FALSE)
#print(xtable::xtable(df), include.rownames=FALSE)
df1 <- df
rownames(df1) <- df1$Type ; df1$Type <- NULL
# write.csv(round(df1, 2), "coefs.csv")

# cumulative returns ----
df.r <- read.csv("data/msci_market_factor_returns.csv")
df.r$Date <- as.Date(df.r$Date)
rownames(df.r) <- df.r$Date
df.r <- df.r[df.r$Date > as.Date("1996-01-01"), ]
df.r <- df.r[, grep("World" ,colnames(df.r))]
colnames(df.r) <- sub("_World", "", colnames(df.r))

df$RF <- 1
rownames(df) <- df$Type
df$Type <- NULL
col <- colnames(df)

df0 <- data.frame((as.matrix(df.r[, col])  %*%  t(as.matrix(df))))
df0 <- df0[, !(colnames(df0) %in% c("PE", "FOF", "SEC"))]
df0 <- df0[, names(sort(apply(df0, 2, function(x) prod(1+x)), decreasing = TRUE))]


# plot it ----
do.eps <- FALSE
if (do.eps) {
  if(TRUE) {
    setEPS()
    postscript(paste0("Cumulative_Returns.eps"), width = 5.5, height = 3, family = "Helvetica", pointsize = 11)
  } else {
    png(filename = "Cumulative_Returns.png",
         width = 480*5, height = 360*5, units = "px", pointsize = 16, res = 400)
    par(lwd=2)
  }
  par(mar=c(4.2,4.2,1,2))
}

plot(as.Date(rownames(df0)), cumprod(1+df0[, 1])-1, type="l", 
     col = "black", ylab = "Cumulative Return", xlab = "Date", ylim = c(0, 20))
i <- 1
for(col in colnames(df0)[-1]) {
  i <- i + 1
  lines(as.Date(rownames(df0)), cumprod(1+df0[, col])-1, col=i)
}
lines(as.Date(rownames(df0)), cumprod(1+df.r$MKT + df.r$RF)-1, lty=3, lwd=3, col = 1)
legend("topleft", bty="n", legend = colnames(df0), col= 1:ncol(df0), lty=1, cex=0.9)
legend("top", bty="n", legend = "MKT", col = 1, lty = 3, lwd = 3, cex = 0.9)

if (do.eps) {
  par(mfrow=c(1,1), cex=1, lwd=2)
  dev.off() 
}
# annulaized returns ----
annualize <- function(x) {(1 + x)^12 - 1}

x.col <- c("MKT", "HML", "SMB", "HDY", "QLT")
df.xr <- data.frame((as.matrix(df.r[, x.col])  %*%  t(as.matrix(df[, x.col]))))
df.xr$MKT <- df.r$MKT

col <- c("MKT", "HML", "SMB", "HDY", "QLT", "RF")
df0 <- data.frame((as.matrix(df.r[, col])  %*%  t(as.matrix(df[, col]))))
df0$MKT <- df.r$MKT + df.r$RF

df.an <- data.frame(
  'Mean' = annualize(colMeans(df0)),
  'Stdv' = sqrt(12) * apply(df0, 2, sd),
  'Mean.ex' = annualize(colMeans(df.xr)),
  'Stdv.ex'= sqrt(12) * apply(df.xr, 2, sd)
)
df.an$Sharpe <- df.an$Mean.ex / df.an$Stdv.ex
print(xtable::xtable(df.an, digits =3, caption = "Annualized returns of five-factor models.",
                     label = "tab:ann_returns"), include.rownames=TRUE)


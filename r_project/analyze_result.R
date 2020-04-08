# analyze results -----
if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

list.cache <- list()
dir.cache <- "data_out/cache_q_factors"
for(file in list.files(dir.cache)) {
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
q.factors <- c("MKT", "ME", "IA", "ROE", "EG")

cv.res <- list()
for(Factor in q.factors) {
  df <- df.f[(df.f$CV.key != "ALL") & (df.f$Factor == Factor), c(q.factors, "Type", "max.quarter")]
  df <- df[, !is.na(df[1, ])]
  df$id <- paste0(df$Type, "_", df$max.quarter)
  df$Type <- NULL
  df$max.quarter <- NULL
  
  if( is.data.frame(df) ) {
    m <-   aggregate(. ~ id, mean, data = df)
    s <-   aggregate(. ~ id, sd, data = df)
  }
  n <- ncol(df)
  a <- 2
  b <- 3
  if(n == b) {
    names(m)[b] <- "Second"
    names(s)[b] <- "Second"
  }
  
  names(m)[a:n] <- paste0(names(m)[a:n], ".mean")
  names(s)[a:n] <- paste0(names(s)[a:n], ".sd")
  
  df <- merge(m, s, by = "id")
  df$Factor <- Factor
  cv.res[[Factor]] <- df
}

df.cv <- data.frame(Reduce(rbind.all.columns, cv.res))
for(i in 1:nrow(df.cv)) {
  df.cv[i, "Type"] <- strsplit(df.cv$id, "_")[[i]][1]
  df.cv[i, "Max.Quarter"] <- strsplit(df.cv$id, "_")[[i]][2]
}
df.cv <- df.cv[, c(7, 8, 2:6)]
df.cv <- df.cv[order(df.cv$Type, df.cv$Max.Quarter), ]
df.cv.40 <- df.cv[df.cv$Max.Quarter == 40, ]
df.cv.60 <- df.cv[df.cv$Max.Quarter == 60, ]
df.cv.40$Max.Quarter <- df.cv.60$Max.Quarter <- NULL
print(xtable::xtable(df.cv.40, caption = "Cross validation with max quarter is 40."), include.rownames=FALSE)
print(xtable::xtable(df.cv.60, caption = "Cross validation with max quarter is 60."), include.rownames=FALSE)


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
print(xtable::xtable(df.all40, caption = "Asymptotic inference with max quarter 40."), include.rownames=FALSE)
print(xtable::xtable(df.all60, caption = "Asymptotic inference with max quarter 60."), include.rownames=FALSE)



# cross validation ----

df.cv <- data.frame(
  aggregate(validation.error ~ Type + lambda, FUN=mean, data=df.f[df.f$CV.key != "ALL", ])
)

df.cv <- df.cv[order(df.cv$validation.error), ]
df.cv <- df.cv[!duplicated(df.cv$Type), ]

df <- df.f[((df.f$CV.key == "ALL") & (df.f$lambda == 2)), ]

aggregate(validation.error ~ Type, FUN=min, data=df.cv[is.finite(df.cv$validation.error), ])

# calc returns ----
type <- "BO"
factor <- "ME"
cols <- c("RF", "MKT", "ME", "IA", "ROE", "EG")
df <- data.frame((as.matrix(df.q[, cols])) %*% t(as.matrix(df.f[(df.f$Type == type) & (df.f[, factor] != 0), cols])))
colnames(df) <- df.f$Factor[df.f$Type == type]

df$Date <- df.q$Date
df <- df[df$Date > as.Date("1989-12-31"), ]

i <- 1
plot(x = df$Date, y = cumprod(1 + df[, i]), xlab = "Date", ylab = "Cum. Return", ylim = c(0, 5000))
max.c <- (ncol(df)-1)
for( i in 2:max.c) {
  lines(df$Date, cumprod(1 + df[, i]), col = "grey")
}
#plot(x = df$Date, y = df[, i], type="l", xlab = "Date", ylab = "Return")


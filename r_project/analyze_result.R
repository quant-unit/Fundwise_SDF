# analyze results -----
if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

df.f <- read.csv("data_out/result_q_factors_4060_2015.csv")
df.f$RF <- 1
df.f[is.na(df.f)] <- 0

df.q <- read.csv("data_prepared/q_factors.csv")
df.q$Date <- as.Date(df.q$Date)

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

# summarize coefficient estimates -----

df.f$Q <- apply(df.f[, c("ME", "IA", "ROE", "EG")], 1, sum)
summi <- function(x) {
  r <- 2
  m <- round(mean(x), r)
  s <- round(sd(x), r)
  t <- round(mean(x) / sd(x), r)
  c(mean = m, t = t)
}
df.summary <- aggregate(. ~ Type + Factor, data=df.f[, c("MKT", "Q", "Type", "Factor")], FUN = summi)
df.summary <- do.call(data.frame, df.summary)
print(xtable::xtable(df.summary), include.rownames=FALSE)

      
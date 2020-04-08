# corona model ----
inf.dyn <- function() {
  population <- 85 * 1000 * 1000
  infections.now <- 77776
  horizon <- 52 # in weeks
  R0_start <- 2.5
  trend <- - (2.5-1) / 52 / sample(seq(0.7, 1.3, 1), 1)
  stdev <- 0.15
  base <- rnorm(horizon, trend, stdev) - (1.01)^(seq(1, horizon * sample(c(1,2), 1)/52, length.out= horizon)) + 1
  
  measure.strength <- trend * seq(2,6)
  measure.weeks <- sample(seq(2,8), 1)
  measure <- sample(measure.strength, measure.weeks, replace = TRUE)
  measure <- c(measure, rep(0, horizon - measure.weeks))
  
  change <- exp(base + measure)
  
  # assumption 1: virus has one week to infect other people
  # assumption 2: all R0 infections in first week
  R0 <- R0_start
  y <- c(infections.now)
  for(x in 1:length(change)) {
    infected.lastweek <- tail(y, 1)
    R0 <- R0 * change[x]
    new <- (population - sum(y)) / population * infected.lastweek * R0
    y <- c(y, new)
  }
  
  #y <- log(y)
  return(y)
}
inf.dyn.runner <- function(o) {
  y <- Inf
  while(sum(y) > 1000 * 1000 * 1000) {
    y <- inf.dyn()
  }
  return(y)
}

# simulate for 1 year horizon ----
set.seed(100)

scenarios <- 1000 * 100
system.time(
  ts.infections <- lapply(1:scenarios, inf.dyn.runner)
)
total.infections <- sapply(ts.infections, sum)  / 1000 / 1000

df <- data.frame(do.call(rbind, ts.infections))
f <- function(x) {
  quantile(x, c(0.5, 0.75, 0.9, 0.99))
}
df.quantile <- apply(df, 2, f)

# analyze weekl infetions over 1 year -----
png(file = "corona_pdf.png", bg = "transparent", width = 480 * 1.25, height = 480)

density.cap <- min(30 * 1000 * 1000, max(unlist(ts.infections)))
plot(ts.infections[[1]], type = "b", 
     ylim = c(0, density.cap), xlim = c(1, 53),
     main = "Weekly infections over 1 year (in Germany, starting 2020-04-01)",
     xlab = "Weeks", ylab = "People infected per week")
for(x in 2:min(2000, scenarios)) {
  lines(ts.infections[[x]], col = "grey")
}
lines(df.quantile[1, ], col = "blue", lwd = 3)
lines(df.quantile[2, ], col = "orange", lwd = 3)
lines(df.quantile[3, ], col = "red", lwd = 3)
lines(df.quantile[4, ], col = "darkred", lwd = 3)
legend("topright", bty="n", cex = 1.4, legend = paste0(c(50, 75, 90, 99), "% Quantile"), lty=1, lwd=3, col = c("blue", "orange", "red", "darkred"))

dev.off()

# analye total infetions after 1 year ----
png(file = "corona_cdf.png", bg = "transparent", width = 480 * 1.25, height = 480)

plot.ecdf(total.infections, xlim = c(0, 90), lwd = 1,
          main = "CDF of Total Infections in 1 Year (in Germany, ending 2021-04-01)",
          ylab = "Probability of less than X infections", 
          xlab = "X in million infected people")

qs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
q.leg <- round(quantile(total.infections, qs),2)
legend("right", bty="n", cex = 1.4, legend = c(paste(names(q.leg), "Quantile:", q.leg)))


quantile(total.infections, seq(0,1,0.01))
summary(total.infections)

dev.off()

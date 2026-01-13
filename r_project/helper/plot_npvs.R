# plot net present value setting ####
## generic cash flow stream ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(latex2exp)

# generate data
set.seed(0)
df.cf <- data.frame(t = seq(0,8), CF = c(-3, -2, -4, -1, 5, 2, 5, 0, 1))
npvs <- data.frame(t = seq(0,8), NPV = sum(df.cf$CF) / 2 + rnorm(9, 0, 0.3))

# plot chart as EPS
do.eps <- FALSE
if(do.eps) {
  setEPS()
  postscript("chart/npvs2.eps", 
             width = 5.5, height = 3, 
             family = "Helvetica", pointsize = 11)
}

par( mar = c(4.2, 4.2, 1, 1) )

plot(df.cf,
     type = "h", xlim = c(-1, 9), ylim = c(-6, 6), ylab = "(Discounted) Cash Flow", xlab = "Time")
abline(h=0, col="grey")
points(npvs, pch = 22, cex = 1, lwd = 1, col= "blue")
points(df.cf, cex = 1, pch= 19, lwd = 1)
legend("bottomright", bty = "n", 
       legend = c(latex2exp::TeX('net cash flow: $\\CF_t$'), 
                  latex2exp::TeX('pricing error: $\\epsilon_{\\tau}$'),
                  latex2exp::TeX('average pricing error: $\\bar{\\epsilon}$')), 
       pch = c(19, 22, NA), lty = c(NA, NA, 1), col = c("black", "blue", "blue"))
avg.npv <- mean(npvs$NPV)
lines(c(0, 8), c(avg.npv, avg.npv), col = "blue")

if(do.eps){ 
  dev.off() 
}

## two deal example ----
# two projects with different investment date 
# only one investment and divestment cash flow per project

# two‐project cash‐flow chart + pricing errors

# define project cash flows
df_proj <- data.frame(
  project = c("A","A","B","B"),
  type    = c("Invest","Divest","Invest","Divest"),
  t       = c(0, 4, 2, 6),
  CF      = c(-3, 5, -4, 6),
  stringsAsFactors = FALSE
)

# generate pricing‐error data (as before)
set.seed(0)
df_all <- data.frame(
  t  = seq(0, 6),
  CF = c(-3, 0, -4, 0, 5, 0, 6)
)
npvs <- data.frame(
  t   = df_all$t,
  NPV = sum(df_all$CF) / 2 + rnorm(nrow(df_all), 0, 0.3)
)


if(do.eps) {
  setEPS()
  postscript("chart/npvs3.eps", 
             width = 5.5, height = 3, 
             family = "Helvetica", pointsize = 10)
}

# plot setup
par(mar = c(4.2, 4.2, 1, 1))
plot(NULL,
     xlim = c(-1, 7),
     ylim = c(min(c(df_proj$CF, npvs$NPV)) - 1,
              max(c(df_proj$CF, npvs$NPV)) + 1),
     xlab = "Time",
     ylab = "(Discounted) Cash Flow")
abline(h = 0, col = "grey")

# color map
cols <- c(A = "red", B = "darkgreen")

# draw project A & B cash flows
for(i in seq_len(nrow(df_proj))) {
  segments(x0 = df_proj$t[i], y0 = 0,
           x1 = df_proj$t[i], y1 = df_proj$CF[i],
           col = cols[df_proj$project[i]], lwd = 2)
  points(df_proj$t[i], df_proj$CF[i],
         pch = 19, col = cols[df_proj$project[i]])
  text(x    = df_proj$t[i],
       y    = df_proj$CF[i] + sign(df_proj$CF[i]) * 0.8,
       labels = paste0(df_proj$project[i], " ", df_proj$type[i]),
       col    = cols[df_proj$project[i]])
}

# add pricing‐error points
points(npvs$t, npvs$NPV,
       pch = 22, cex = 1, lwd = 1, col = "blue")

# average pricing‐error line
avg.npv <- mean(npvs$NPV)
lines(c(0, 6), c(avg.npv, avg.npv), col = "blue")

# legend
legend("bottomright", bty = "n",
       legend = c("Project A cash flows", 
                  "Project B cash flows",
                  latex2exp::TeX('Pricing error per date: $\\epsilon_{\\tau}$'),
                  latex2exp::TeX('Average pricing error: $\\bar{\\epsilon}$')),
       pch    = c(19, 19, 22, NA),
       lty    = c(NA, NA, NA, 1),
       col    = c(cols["A"], cols["B"], "blue", "blue"))

if(do.eps){ 
  dev.off() 
}


# plot net present value setting ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(latex2exp)
set.seed(0)
df.cf <- data.frame(t = seq(0,8), CF = c(-3, -2, -4, -1, 5, 2, 5, 0, 1))
npvs <- data.frame(t = seq(0,8), NPV = sum(df.cf$CF) / 2 + rnorm(9, 0, 0.3))

do.eps <- TRUE
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

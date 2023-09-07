# Prologue ----
# https://www.r-bloggers.com/basic-waterfall-graphs-in-r/
rm(list=ls()) # remove workspace objects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
library(dplyr)
library(ggplot2)
# SET OPTIONS
do.eps <- FALSE
include.irr <- TRUE
sdf.type <- "best"
period <- 5

portfolio <- "Seafarers.full"
portfolio <- "Casagrande.L2.full"
portfolio <- "Certior.L2.full"
portfolio <- "VER_Inflexion 2010 Buyout Fund"
portfolio <- "VER_Inveni Life Science Fund I Ky (VC EU)"
portfolio <- "Unigestion.EC3"
portfolio <- "veritas.full" 

name <- switch(portfolio, 
               "Unigestion.EC3" = "40_EU_BO",
               "Seafarers.full" = "99_Funds",
               "Casagrande.L2.full" = "29_Funds",
               "Certior.L2.full" = "20_PD_Funds",
               "VER_Inveni Life Science Fund I Ky (VC EU)" = "VC_EU",
               "VER_Inflexion 2010 Buyout Fund" = "BO_EU",
               "veritas.full" = "100_Pofo")
period <- switch(portfolio, 
               "Unigestion.EC3" = 15,
               "Seafarers.full" = 15,
               "Casagrande.L2.full" = 15,
               "Certior.L2.full" = 5,
               "VER_Inveni Life Science Fund I Ky (VC EU)" = 10,
               "VER_Inflexion 2010 Buyout Fund" = 10,
               "veritas.full" = 15)

# Get data ------
# as of 2019-12-31 / Region: Europe
df.average.return <- data.frame(Factor = c("RF", "MKT.RF", "SMB", "HML", "QLT.MKT", "HDY.MKT", "COST"))
df.average.return$average15 <- c(0.01411958,0.062593222,0.05868715,-0.006967423,0.044793267,0.016339359,-0.02)
df.average.return$average10 <- c(0.006017912,0.085588926,0.075695219,-0.004235658,0.05535489,0.029023313,-0.02)
df.average.return$average5 <- c(0.010913105,0.070146568,0.057489858,-0.020186355,0.043768153,0.021167862,-0.02)
df.average.return$average50 <- c(0.047772926,0.068512997,0.029043128,0.016431317,0.018717367,0.015804019,-0.02)

if (portfolio == "Unigestion.EC3") {
  df <- data.frame(RF = c(1, 1),
                   'MKT.RF' = c(1.087301001, 1.06509775), 
                   SMB = c(0.460365778, 0.165908979), 
                   HML = c(0.11626654, 0.089034275), 
                   'QLT.MKT' = c(0.060883247, 0.064590687), 
                   'HDY.MKT' = c(0.180889146, 0.197375191),
                   COST = c(1, 1))
  irr <- list()
  irr[["best"]][["mPME"]] <- 0.117888414
  irr[["balanced"]][["mPME"]] <- 0.091960254
  irr[["realized"]] <- 0.088225853

}

if (portfolio == "Seafarers.full") {
  df <- data.frame(RF = c(1, 1),
                   'MKT.RF' = c(0.987319177, 0.815990412), 
                   SMB = c(0.213771369, 0.108750749), 
                   HML = c(0.129148449, 0.069690094), 
                   'QLT.MKT' = c(0.184503545, 0.062332091), 
                   'HDY.MKT' = c(0.390942624, 0.159740729),
                   COST = c(1, 1))
  irr <- list()
  irr[["best"]][["mPME"]] <- 0.078915373
  irr[["balanced"]][["mPME"]] <- 0.052880812
  irr[["realized"]] <- 0.066413344
}

if (portfolio == "Casagrande.L2.full") {
  df <- data.frame(RF = c(1, 1),
                   'MKT.RF' = c(0.789611175, 0.950554896), 
                   SMB = c(0.06134811, 0.135644578), 
                   HML = c(0.126306046, 0.071922732), 
                   'QLT.MKT' = c(0.026360618, 0.060392933), 
                   'HDY.MKT' = c(0.31687394, 0.169181574),
                   COST = c(1, 1))
  irr <- list()
  irr[["best"]][["mPME"]] <- 0.060925164
  irr[["balanced"]][["mPME"]] <- 0.075593808
  irr[["realized"]] <- 0.031143312
}

if (portfolio == "Certior.L2.full") {
  df <- data.frame(RF = c(1, 1),
                   'MKT.RF' = c(0.990469925, 0.779636548), 
                   SMB = c(0.313196402, 0.12307946), 
                   HML = c(0.369853338, 0.09270626), 
                   'QLT.MKT' = c(0.104427164, 0.042912642), 
                   'HDY.MKT' = c(0.30330726, 0.163042788),
                   COST = c(1, 1))
  irr <- list()
  irr[["best"]][["mPME"]] <- 0.058573841
  irr[["balanced"]][["mPME"]] <- 0.049532156
  irr[["realized"]] <- 0.084608483
}


if (portfolio == "VER_Inflexion 2010 Buyout Fund") {
  df <- data.frame(RF = c(1, 1),
                   'MKT.RF' = c(1.307901641, 1.090380151), 
                   SMB = c(1.226881555, 0.169702255), 
                   HML = c(0, 0.088709059), 
                   'QLT.MKT' = c(0.053016844, 0.066510643), 
                   'HDY.MKT' = c(1.055844559, 0.200415904),
                   COST = c(1, 1))
  
  irr <- list()
  irr[["best"]][["mPME"]] <- 0.30326749144699
  irr[["balanced"]][["mPME"]] <- 0.111718719345244
  irr[["realized"]] <- 0.279108054379121
}



if (portfolio == "VER_Inveni Life Science Fund I Ky (VC EU)") {
  df <- data.frame(RF = c(1, 1),
                         'MKT-RF' = c(0.711138396, 1.058079251), 
                         SMB = c(-0.187615279, -0.052997201), 
                         HML = c(-0.02140567, -0.071975159), 
                         'QLT-MKT' = c(0, 0.131547401), 
                         'HDY-MKT' = c(-0.027233445, -0.035453423),
                         COST = c(1, 1))
  irr <- list()
  irr[["best"]][["mPME"]] <- 0.030692279
  irr[["balanced"]][["mPME"]] <- 0.072784701
  irr[["realized"]] <- -0.079196334837057
}



if (portfolio == "veritas.full") {
  df <- data.frame(RF = c(1, 1),
                   'MKT-RF' = c(1.170892538, 1.032887102), 
                   SMB = c(0.211346749, 0.138343884), 
                   HML = c(0.100534286, 0.094225983), 
                   'QLT-MKT' = c(0.256032758, 0.059147616), 
                   'HDY-MKT' = c(0.242805277, 0.181054425),
                   COST = c(1, 1))
  irr <- list()
  irr[["best"]][["mPME"]] <- 0.139695039
  irr[["balanced"]][["mPME"]] <- 0.109698499
  irr[["realized"]] <- 0.125192054809775
}



df[3, ] <- df.average.return[, paste0("average", period)]
rownames(df) <- c("best-ensemble coef", "balanced-ensemble coef", "average factor return")
df["best", ] <- df["best-ensemble coef", ] * df["average factor return", ]
df["balanced", ] <- df["balanced-ensemble coef", ] * df["average factor return", ]
df

# Prepare data for waterfall ----
df1 <- data.frame(t(df[sdf.type, ]))
colnames(df1) <- "Value"
df1["Expected Return", ] <- sum(df1$Value)
df1$Category <- rownames(df1)
df1

data1 <- df1  %>%
  mutate(Category = factor(Category, levels = df1$Category),
         ymin = cumsum(Value),
         ymax = lag(cumsum(Value), default = 0),
         xmin = c(head(Category, -1), NA),
         xmax = c(tail(Category, -1), NA),
         Impact = ifelse(Category %in% c("Expected Return"),"Time-weighted Return",
                         ifelse(Value > 0, "Increase", "Decrease")
         ))
data1[nrow(data1), c("ymin", "xmin", "xmax")] <- c(0, 8, 9)

if (include.irr) {
  data1$Category <- as.character(data1$Category)
  data1[(nrow(data1)+1), ] <- c(0, "Timing Gap", irr[[sdf.type]], 
                                data1$ymax[data1$Category == "Expected Return"], 9, 10, "Timing Gap")
  
  data1[(nrow(data1)+1), ] <- c(0, "IRR mPME", 0, irr[[sdf.type]], 9, 10, "Money-weighted Return")
  
  illiquidity.premium  <- 0.01
  data1[(nrow(data1)+1), ] <- c(0, "Illiquidity", irr[[sdf.type]], 
                                irr[[sdf.type]] + illiquidity.premium, 10, 11, "Replication Gap")
  private.cost <- - 0.01
  data1[(nrow(data1)+1), ] <- c(0, "Private Cost", irr[[sdf.type]] + illiquidity.premium, 
                                irr[[sdf.type]] + illiquidity.premium + private.cost, 11, 12, "Replication Gap")
  
  residual <- irr[[sdf.type]] + illiquidity.premium + private.cost - irr$realized
  data1[(nrow(data1)+1), ] <- c(0, "Residual", irr$realized, irr$realized + residual, 12, 13, "Replication Gap")
  
  data1[(nrow(data1)+1), ] <- c(0, "IRR realized", 0, irr$realized, NA, NA, "Money-weighted Return")
  data1$Category[data1$Category == "COST"] <- "Public Cost"
  data1$Category <- factor(data1$Category, data1$Category)
  data1$xmin <- 1:nrow(data1)
  data1$xmax <- data1$xmin + 1
  for (col in c("Value", "ymin", "ymax")) data1[, col] <- as.numeric(as.character(data1[, col]))
  for (col in c("xmin", "xmax")) data1[, col] <- as.integer(as.character(data1[, col]))
}

str(data1)
data1

# Plot waterfall ----
title <- ifelse(include.irr, portfolio, name)

g <- ggplot(data1) +
  theme_bw()+
  theme(legend.position = "right", 
        panel.grid = element_blank(), 
        plot.title = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "Return (p.a.)", x = "Return Components", 
       title = paste0(title, " (", sdf.type, ", ", period, " year average)"))

w <- 0.8  #use to set width of bars

g <- g +
  geom_rect(aes(xmin = as.integer(Category) - w/2,
                xmax = as.integer(Category) + w/2, ymin = ymin, ymax = ymax,
                fill = Impact), colour = NA) +
  scale_x_discrete(limits = data1$Category) +
  geom_hline(yintercept = 0, colour = "grey", linetype=3) +
  scale_fill_manual(values = (c("Increase" = "green", "Decrease" = "red", 
                                "Timing Gap" = "blue",
                                "Replication Gap" = "navy",
                                "Time-weighted Return" = "gray60",
                                "Money-weighted Return" = "gray40")))
g

if (do.eps) {
  setEPS()
  postscript(paste0("Waterfall_", name,".eps"), width = 5, height = 3, family = "Helvetica", pointsize = 7)
  plot(g)
  par(mfrow=c(1,1), cex=1)
  dev.off() 
  dev.off() 
}


# Plot coefs ----
library(latex2exp)
df2 <- data.frame(t(df[c("best-ensemble coef", "balanced-ensemble coef"), ]))
df2$Factors <- rownames(df2)
df2$fill.col <- "abc"
df2 <- df2[!(rownames(df2) %in% c("ALPHA", "RF", "COST")), ]
df2$Factors <- c("MKT-RF", "SMB", "HML", "QLT-MKT", "HDY-MKT")
df2

do.eps <- TRUE
if (do.eps) {
  setEPS()
  postscript(paste0("Coefs_", name,".eps"), width = 5.5, height = 2.5, family = "Helvetica", pointsize = 11)
}

par(mar = c(5,6,2,4), xpd=TRUE)
df.bar <- barplot(df2$best.ensemble.coef, xlim = c(0,1.2), names.arg = df2$Factors, horiz = TRUE, las=2, 
                  xlab= "Coefficient Estimate")
points(y = df.bar, x = df2$balanced.ensemble.coef, col = "red", pch= 16, cex=1.5)
legend("bottomright", bty="n", legend = c("Best Ensemble: M**", "Balanced Ensemble: M*"), 
       fill=c("grey", NA), border = c("black", NA), col = c(NA, "red"), pch = c(NA, 19),
       inset=c(0,1), xpd=TRUE, horiz=TRUE)

if (do.eps) {
  par(mfrow=c(1,1), cex=1)
  dev.off() 
}


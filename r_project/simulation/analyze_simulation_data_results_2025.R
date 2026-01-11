## analyze simulation ----
if(sys.nframe() == 0L) rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

bivar <- function(files, digi=3) {
  l <- list()
  for (file in files) {
    name <- paste0("data_results_sim/data_out_2025/", file, ".csv")
    l[[file]] <- read.csv(name)
  }
  df0 <- data.frame(do.call(rbind, l))
  
  print(nrow(df0))
  print(paste("max.month:", df0$max.month[1]))
  
  res <- list()
  for (factor in levels(as.factor(df0$Factor))) {
    coefs <- "MKT"
    if(factor != "MKT") coefs <- c(coefs, factor)
    print(coefs)
    
    df <- df0[df0$Factor == factor, ]
    print(paste("# of models:", nrow(df)))
    for(coef in coefs) {
      plot(hist(df[, coef]), main = paste(factor, coef))
      m <- mean(df[, coef])
      s <- sd(df[, coef])
      x <- ifelse(coef == "Alpha", 100, 1)
      res[[factor]][[paste0(coef)]] <- list(Mean = m*x, SD = s*x)
      y <- paste0(factor, "-", coef, ":  ", round(m*x, digi), " (", round(s*x, digi), ")      ")
      print(y)
    }
  }
  invisible(res)
}

l <- list()
## cross-sectional unit -----
# vintage year portfolios
run <- "cache_q_factors_20250808_222540_simulated_cashflows_EW_VYP/"
l$sl_1$one <-   bivar(paste0(run, "2025-08-12_010504_cached_res")) # 1
l$sl_60$one <-  bivar(paste0(run, "2025-08-12_010958_cached_res")) # 60
l$sl_120$one <- bivar(paste0(run, "2025-08-12_011508_cached_res")) # 120
l$sl_150$one <- bivar(paste0(run, "2025-08-12_012028_cached_res")) # 150
l$sl_180$one <- bivar(paste0(run, "2025-08-12_012558_cached_res")) # 180
l$sl_210$one <- bivar(paste0(run, "2025-08-12_013127_cached_res")) # 210
l$sl_240$one <- bivar(paste0(run, "2025-08-12_013656_cached_res")) # 240
l$sl_300$one <- bivar(paste0(run, "2025-08-12_014225_cached_res")) # 300
l$sl_360$one <- bivar(paste0(run, "2025-08-12_014752_cached_res")) # 360

# individual funds (new in 2025: also 1000 iterations)
run <- "cache_q_factors_20250808_222540_simulated_cashflows_EW/"
l$sl_1$funds <- bivar(
  c(
    paste0(run, "2025-08-13_000858_cached_res"),
    paste0(run, "2025-08-13_011913_cached_res"),
    paste0(run, "2025-08-13_021810_cached_res"),
    paste0(run, "2025-08-13_031732_cached_res"),
    paste0(run, "2025-08-13_041901_cached_res"),
    paste0(run, "2025-08-13_051818_cached_res"),
    paste0(run, "2025-08-13_061953_cached_res"),
    paste0(run, "2025-08-13_072647_cached_res"),
    paste0(run, "2025-08-13_083327_cached_res"),
    paste0(run, "2025-08-13_093948_cached_res")
                        )
  ) # 1
l$sl_60$funds <- bivar(
  c(
    paste0(run, "2025-08-13_001431_cached_res"),
    paste0(run, "2025-08-13_012439_cached_res"),
    paste0(run, "2025-08-13_022345_cached_res"),
    paste0(run, "2025-08-13_032312_cached_res"),
    paste0(run, "2025-08-13_042436_cached_res"),
    paste0(run, "2025-08-13_052401_cached_res"),
    paste0(run, "2025-08-13_062548_cached_res"),
    paste0(run, "2025-08-13_073255_cached_res"),
    paste0(run, "2025-08-13_083945_cached_res"),
    paste0(run, "2025-08-13_094617_cached_res")
  )
) # 60
l$sl_120$funds <- bivar(
  c(
    paste0(run, "2025-08-13_002013_cached_res"),
    paste0(run, "2025-08-13_013022_cached_res"),
    paste0(run, "2025-08-13_022926_cached_res"),
    paste0(run, "2025-08-13_032905_cached_res"),
    paste0(run, "2025-08-13_043011_cached_res"),
    paste0(run, "2025-08-13_052953_cached_res"),
    paste0(run, "2025-08-13_063210_cached_res"),
    paste0(run, "2025-08-13_073915_cached_res"),
    paste0(run, "2025-08-13_084604_cached_res"),
    paste0(run, "2025-08-13_095310_cached_res")
  )
) # 120
l$sl_150$funds <- bivar(
  c(
    paste0(run, "2025-08-13_002637_cached_res"),
    paste0(run, "2025-08-13_013633_cached_res"),
    paste0(run, "2025-08-13_023540_cached_res"),
    paste0(run, "2025-08-13_033535_cached_res"),
    paste0(run, "2025-08-13_043621_cached_res"),
    paste0(run, "2025-08-13_053619_cached_res"),
    paste0(run, "2025-08-13_063914_cached_res"),
    paste0(run, "2025-08-13_074617_cached_res"),
    paste0(run, "2025-08-13_085254_cached_res"),
    paste0(run, "2025-08-13_100035_cached_res")
  )
) # 150
l$sl_180$funds <- bivar(
  c(
    paste0(run, "2025-08-13_003314_cached_res"),
    paste0(run, "2025-08-13_014253_cached_res"),
    paste0(run, "2025-08-13_024200_cached_res"),
    paste0(run, "2025-08-13_034214_cached_res"),
    paste0(run, "2025-08-13_044239_cached_res"),
    paste0(run, "2025-08-13_054258_cached_res"),
    paste0(run, "2025-08-13_064634_cached_res"),
    paste0(run, "2025-08-13_075332_cached_res"),
    paste0(run, "2025-08-13_090005_cached_res"),
    paste0(run, "2025-08-13_100819_cached_res")
  )
) # 180
l$sl_210$funds <- bivar(
  c(
    paste0(run, "2025-08-13_003953_cached_res"),
    paste0(run, "2025-08-13_014906_cached_res"),
    paste0(run, "2025-08-13_024819_cached_res"),
    paste0(run, "2025-08-13_034853_cached_res"),
    paste0(run, "2025-08-13_044857_cached_res"),
    paste0(run, "2025-08-13_054937_cached_res"),
    paste0(run, "2025-08-13_065353_cached_res"),
    paste0(run, "2025-08-13_080047_cached_res"),
    paste0(run, "2025-08-13_090715_cached_res"),
    paste0(run, "2025-08-13_101603_cached_res")
  )
) # 210
l$sl_240$funds <- bivar(
  c(
    paste0(run, "2025-08-13_004632_cached_res"),
    paste0(run, "2025-08-13_015520_cached_res"),
    paste0(run, "2025-08-13_025439_cached_res"),
    paste0(run, "2025-08-13_035531_cached_res"),
    paste0(run, "2025-08-13_045515_cached_res"),
    paste0(run, "2025-08-13_055616_cached_res"),
    paste0(run, "2025-08-13_070118_cached_res"),
    paste0(run, "2025-08-13_080802_cached_res"),
    paste0(run, "2025-08-13_091425_cached_res"),
    paste0(run, "2025-08-13_102348_cached_res")
  )
) # 240
l$sl_300$funds <- bivar(
  c(
    paste0(run, "2025-08-13_005315_cached_res"),
    paste0(run, "2025-08-13_020133_cached_res"),
    paste0(run, "2025-08-13_030058_cached_res"),
    paste0(run, "2025-08-13_040210_cached_res"),
    paste0(run, "2025-08-13_050133_cached_res"),
    paste0(run, "2025-08-13_060254_cached_res"),
    paste0(run, "2025-08-13_070854_cached_res"),
    paste0(run, "2025-08-13_081518_cached_res"),
    paste0(run, "2025-08-13_092135_cached_res"),
    paste0(run, "2025-08-13_103133_cached_res")
  )
) # 300
l$sl_360$funds <- bivar(
  c(
    paste0(run, "2025-08-13_010004_cached_res"),
    paste0(run, "2025-08-13_020752_cached_res"),
    paste0(run, "2025-08-13_030719_cached_res"),
    paste0(run, "2025-08-13_040849_cached_res"),
    paste0(run, "2025-08-13_050751_cached_res"),
    paste0(run, "2025-08-13_060933_cached_res"),
    paste0(run, "2025-08-13_071614_cached_res"),
    paste0(run, "2025-08-13_082233_cached_res"),
    paste0(run, "2025-08-13_092844_cached_res"),
    paste0(run, "2025-08-13_103918_cached_res")
  )
) # 360

# Aggregate for plotting
months <- c(1, 60, 120, 150, 180, 210, 240, 300, 360)
df1 <- data.frame(
  months = months,
  VYP_mean = sapply(months, function(x) (l[[paste0("sl_", x)]]$one$MKT$MKT$Mean)),
  VYP_sd = sapply(months, function(x) (l[[paste0("sl_", x)]]$one$MKT$MKT$SD)),
  Funds_mean = sapply(months, function(x) (l[[paste0("sl_", x)]]$funds$MKT$MKT$Mean)),
  Funds_sd = sapply(months, function(x) (l[[paste0("sl_", x)]]$funds$MKT$MKT$SD))
)
df1


# Plot
do.eps <- FALSE
if(do.eps) {
  setEPS()
  postscript("chart/Simulation_funds_vs_vyps.eps", 
             width = 5.5, height = 3, 
             family = "Helvetica", pointsize = 10)
}

par( mar = c(4.2, 4.2, 1, 1) )

plot(df1$months, df1$VYP_mean, ylim = c(0, 2), type = "b",
     xlab = "Max Months", ylab = "Mean and Stdv of Estimates")
lines(df1$months, df1$Funds_mean, type = "b", col = "blue")
lines(df1$months, df1$VYP_sd, type = "b", lty = 3)
lines(df1$months, df1$Funds_sd, type = "b", lty = 3, col = "blue")
abline(h=c(0, 1), col="grey")
legend("topright", bty = "n", legend = c("Mean Funds", "Mean VYP", "Stdv Funds", "Stdv VYP"), 
       col = c("blue", "black"), lty = c(1,1,3,3), cex = 0.8)

if(do.eps){ 
  dev.off() 
}

## Shorter fund lifetime max(5 + 5) ----
run <- "cache_q_factors_20250808_223358_simulated_cashflows_EW_VYP/"
l$sl_1$short <-   bivar(paste0(run, "2025-08-13_104842_cached_res")) # 1
l$sl_60$short <-  bivar(paste0(run, "2025-08-13_105657_cached_res")) # 60
l$sl_120$short <- bivar(paste0(run, "2025-08-13_110543_cached_res")) # 120
l$sl_150$short <- bivar(paste0(run, "2025-08-13_111426_cached_res")) # 150
l$sl_180$short <- bivar(paste0(run, "2025-08-13_112309_cached_res")) # 180
l$sl_210$short <- bivar(paste0(run, "2025-08-13_113152_cached_res")) # 210
l$sl_240$short <- bivar(paste0(run, "2025-08-13_114034_cached_res")) # 240
l$sl_300$short <- bivar(paste0(run, "2025-08-13_114917_cached_res")) # 300
l$sl_360$short <- bivar(paste0(run, "2025-08-13_115806_cached_res")) # 360

## High Beta, Negative Alpha ------
# SL (Beta 2.5, Alpha -0.0025) # -3% p.a. Alpha
run <- "cache_q_factors_20250808_224228_simulated_cashflows_EW_VYP/"

l$sl_1$high <-   bivar(paste0(run, "2025-08-13_141005_cached_res")) # 1
l$sl_60$high <-  bivar(paste0(run, "2025-08-13_142305_cached_res")) # 60
l$sl_120$high <- bivar(paste0(run, "2025-08-13_143608_cached_res")) # 120
l$sl_150$high <- bivar(paste0(run, "2025-08-13_144833_cached_res")) # 150
l$sl_180$high <- bivar(paste0(run, "2025-08-13_150111_cached_res")) # 180
l$sl_210$high <- bivar(paste0(run, "2025-08-13_151331_cached_res")) # 210
l$sl_240$high <- bivar(paste0(run, "2025-08-13_152547_cached_res")) # 240
l$sl_300$high <- bivar(paste0(run, "2025-08-13_153800_cached_res")) # 300
l$sl_360$high <- bivar(paste0(run, "2025-08-13_155012_cached_res")) # 360


## Exp aff models ----
# EA, true DGP (Beta 1.0, Alpha 0)
run <- "cache_q_factors_20250808_225102_simulated_cashflows_EW_VYP/"

l$ea_1$one <-   bivar(paste0(run, "2025-08-13_160317_cached_res")) # 1
l$ea_60$one <-  bivar(paste0(run, "2025-08-13_161507_cached_res")) # 60
l$ea_120$one <- bivar(paste0(run, "2025-08-13_162747_cached_res")) # 120
l$ea_150$one <- bivar(paste0(run, "2025-08-13_164108_cached_res")) # 150
l$ea_180$one <- bivar(paste0(run, "2025-08-13_165455_cached_res")) # 180
l$ea_210$one <- bivar(paste0(run, "2025-08-13_170843_cached_res")) # 210
l$ea_240$one <- bivar(paste0(run, "2025-08-13_172231_cached_res")) # 240
l$ea_300$one <- bivar(paste0(run, "2025-08-13_173617_cached_res")) # 300
l$ea_360$one <- bivar(paste0(run, "2025-08-13_175003_cached_res")) # 360

# EA (Beta 2.5, Alpha -0.0025) # -3% p.a. Alpha
run <- "cache_q_factors_20250808_225920_simulated_cashflows_EW_VYP/"

l$ea_1$high <-   bivar(paste0(run, "2025-08-13_180629_cached_res")) # 1
l$ea_60$high <-  bivar(paste0(run, "2025-08-13_183031_cached_res")) # 60
l$ea_120$high <- bivar(paste0(run, "2025-08-13_185404_cached_res")) # 120
l$ea_150$high <- bivar(paste0(run, "2025-08-13_191738_cached_res")) # 150
l$ea_180$high <- bivar(paste0(run, "2025-08-13_194040_cached_res")) # 180
l$ea_210$high <- bivar(paste0(run, "2025-08-13_200339_cached_res")) # 210
l$ea_240$high <- bivar(paste0(run, "2025-08-13_202648_cached_res")) # 240
l$ea_300$high <- bivar(paste0(run, "2025-08-13_204954_cached_res")) # 300
l$ea_360$high <- bivar(paste0(run, "2025-08-13_211258_cached_res")) # 360

## Sim Lin & Exp aff summary ------

make.df <- function(m = c("sl", "ea")) {
  months <- c(1, 1, 60, 60, 120, 120, 180, 180, 240, 240, 300, 300, 360, 360)
  df <- data.frame(MaxMonth = months)
  mean.sd <- unlist(lapply(1:(nrow(df)/2), function(x) return(c("mean", "stdv"))))
  df$MaxMonth <- paste(df$MaxMonth, "-", mean.sd)
  
  # beta 1
  df$beta1true <- c(l[[paste0(m, "_1")]]$one$MKT$MKT$Mean, l[[paste0(m, "_1")]]$one$MKT$MKT$SD,
                    l[[paste0(m, "_60")]]$one$MKT$MKT$Mean, l[[paste0(m, "_60")]]$one$MKT$MKT$SD,
                    l[[paste0(m, "_120")]]$one$MKT$MKT$Mean, l[[paste0(m, "_120")]]$one$MKT$MKT$SD,
                    l[[paste0(m, "_180")]]$one$MKT$MKT$Mean, l[[paste0(m, "_180")]]$one$MKT$MKT$SD,
                    l[[paste0(m, "_240")]]$one$MKT$MKT$Mean, l[[paste0(m, "_240")]]$one$MKT$MKT$SD,
                    l[[paste0(m, "_300")]]$one$MKT$MKT$Mean, l[[paste0(m, "_300")]]$one$MKT$MKT$SD,
                    l[[paste0(m, "_360")]]$one$MKT$MKT$Mean, l[[paste0(m, "_360")]]$one$MKT$MKT$SD)
  
  df$alpha1false <- c(l[[paste0(m, "_1")]]$one$Alpha$Alpha$Mean, l[[paste0(m, "_1")]]$one$Alpha$Alpha$SD,
                      l[[paste0(m, "_60")]]$one$Alpha$Alpha$Mean, l[[paste0(m, "_60")]]$one$Alpha$Alpha$SD,
                      l[[paste0(m, "_120")]]$one$Alpha$Alpha$Mean, l[[paste0(m, "_120")]]$one$Alpha$Alpha$SD,
                      l[[paste0(m, "_180")]]$one$Alpha$Alpha$Mean, l[[paste0(m, "_180")]]$one$Alpha$Alpha$SD,
                      l[[paste0(m, "_240")]]$one$Alpha$Alpha$Mean, l[[paste0(m, "_240")]]$one$Alpha$Alpha$SD,
                      l[[paste0(m, "_300")]]$one$Alpha$Alpha$Mean, l[[paste0(m, "_300")]]$one$Alpha$Alpha$SD,
                      l[[paste0(m, "_360")]]$one$Alpha$Alpha$Mean, l[[paste0(m, "_360")]]$one$Alpha$Alpha$SD)
  
  df$beta1false <- c(l[[paste0(m, "_1")]]$one$Alpha$MKT$Mean, l[[paste0(m, "_1")]]$one$Alpha$MKT$SD,
                     l[[paste0(m, "_60")]]$one$Alpha$MKT$Mean, l[[paste0(m, "_60")]]$one$Alpha$MKT$SD,
                     l[[paste0(m, "_120")]]$one$Alpha$MKT$Mean, l[[paste0(m, "_120")]]$one$Alpha$MKT$SD,
                     l[[paste0(m, "_180")]]$one$Alpha$MKT$Mean, l[[paste0(m, "_180")]]$one$Alpha$MKT$SD,
                     l[[paste0(m, "_240")]]$one$Alpha$MKT$Mean, l[[paste0(m, "_240")]]$one$Alpha$MKT$SD,
                     l[[paste0(m, "_300")]]$one$Alpha$MKT$Mean, l[[paste0(m, "_300")]]$one$Alpha$MKT$SD,
                     l[[paste0(m, "_360")]]$one$Alpha$MKT$Mean, l[[paste0(m, "_360")]]$one$Alpha$MKT$SD)
  
  # high beta
  df$beta25false <- c(l[[paste0(m, "_1")]]$high$MKT$MKT$Mean, l[[paste0(m, "_1")]]$high$MKT$MKT$SD,
                      l[[paste0(m, "_60")]]$high$MKT$MKT$Mean, l[[paste0(m, "_60")]]$high$MKT$MKT$SD,
                      l[[paste0(m, "_120")]]$high$MKT$MKT$Mean, l[[paste0(m, "_120")]]$high$MKT$MKT$SD,
                      l[[paste0(m, "_180")]]$high$MKT$MKT$Mean, l[[paste0(m, "_180")]]$high$MKT$MKT$SD,
                      l[[paste0(m, "_240")]]$high$MKT$MKT$Mean, l[[paste0(m, "_240")]]$high$MKT$MKT$SD,
                      l[[paste0(m, "_300")]]$high$MKT$MKT$Mean, l[[paste0(m, "_300")]]$high$MKT$MKT$SD,
                      l[[paste0(m, "_360")]]$high$MKT$MKT$Mean, l[[paste0(m, "_360")]]$high$MKT$MKT$SD)
  
  df$alpha25true <- c(l[[paste0(m, "_1")]]$high$Alpha$Alpha$Mean, l[[paste0(m, "_1")]]$high$Alpha$Alpha$SD,
                      l[[paste0(m, "_60")]]$high$Alpha$Alpha$Mean, l[[paste0(m, "_60")]]$high$Alpha$Alpha$SD,
                      l[[paste0(m, "_120")]]$high$Alpha$Alpha$Mean, l[[paste0(m, "_120")]]$high$Alpha$Alpha$SD,
                      l[[paste0(m, "_180")]]$high$Alpha$Alpha$Mean, l[[paste0(m, "_180")]]$high$Alpha$Alpha$SD,
                      l[[paste0(m, "_240")]]$high$Alpha$Alpha$Mean, l[[paste0(m, "_240")]]$high$Alpha$Alpha$SD,
                      l[[paste0(m, "_300")]]$high$Alpha$Alpha$Mean, l[[paste0(m, "_300")]]$high$Alpha$Alpha$SD,
                      l[[paste0(m, "_360")]]$high$Alpha$Alpha$Mean, l[[paste0(m, "_360")]]$high$Alpha$Alpha$SD)
  
  df$beta25true <- c(l[[paste0(m, "_1")]]$high$Alpha$MKT$Mean, l[[paste0(m, "_1")]]$high$Alpha$MKT$SD,
                     l[[paste0(m, "_60")]]$high$Alpha$MKT$Mean, l[[paste0(m, "_60")]]$high$Alpha$MKT$SD,
                     l[[paste0(m, "_120")]]$high$Alpha$MKT$Mean, l[[paste0(m, "_120")]]$high$Alpha$MKT$SD,
                     l[[paste0(m, "_180")]]$high$Alpha$MKT$Mean, l[[paste0(m, "_180")]]$high$Alpha$MKT$SD,
                     l[[paste0(m, "_240")]]$high$Alpha$MKT$Mean, l[[paste0(m, "_240")]]$high$Alpha$MKT$SD,
                     l[[paste0(m, "_300")]]$high$Alpha$MKT$Mean, l[[paste0(m, "_300")]]$high$Alpha$MKT$SD,
                     l[[paste0(m, "_360")]]$high$Alpha$MKT$Mean, l[[paste0(m, "_360")]]$high$Alpha$MKT$SD)
  
  print(xtable::xtable(df, digits = 3, caption = "Simulation study for max month, exp aff.", label = "tab:simulation_study_max_month"), include.rownames=FALSE)
  
  df$months <- months
  df$mean.sd <- mean.sd
  invisible(df)
}

df.sl <- make.df("sl")
df.ea <- make.df("ea")

# plot
do.eps <- FALSE
if(do.eps) {
  setEPS()
  postscript("chart/Simulation_expaff_vs_simlin.eps", 
             width = 5.5, height = 3, 
             family = "Helvetica", pointsize = 10)
}

par( mar = c(4.2, 4.2, 1, 1) )

plot(df.sl[df.sl$mean.sd == "mean", c("months", "beta1true")], ylim = c(0,2), type = "b",
     xlab = "Max Months", ylab = "Mean and Stdv of Estimates")
abline(h = c(0,1), col = "grey")

lines(df.sl[df.sl$mean.sd == "stdv", c("months", "beta1true")], type = "b", lty = 3)
lines(df.ea[df.ea$mean.sd == "mean", c("months", "beta1true")], type = "b", col = "red")
lines(df.ea[df.ea$mean.sd == "stdv", c("months", "beta1true")], type = "b", lty = 3, col = "red")

legend("topright", bty = "n", 
       legend = c("Mean Exp.Aff. SDF", "Mean Linear SDF", "Stdv Exp.Aff. SDF", "Stdv Linear SDF"), 
       col = c("red", "black"), lty = c(1,1,3,3), cex = 0.8)

if(do.eps){ 
  dev.off() 
}

## Double half models (500 iter) ----
h <- list()

## Same, same with 1000 iters (base case)
run <- "cache_q_factors_20250808_222540_simulated_cashflows_EW_VYP/"
h$base <- bivar(paste0(run, "2025-08-12_012558_cached_res")) # 180

## Double vintages, half funds per vintage (10)
run <- "cache_q_factors_20250808_232503_simulated_cashflows_EW_VYP/"
h$double.vin.half.funds <- bivar(paste0(run, "2025-08-13_235242_cached_res")) # 180

## Double vintages, same funds per vintage (20)
run <- "cache_q_factors_20250813_231430_simulated_cashflows_EW_VYP/"
h$double.vin.same.funds <- bivar(paste0(run, "2025-08-14_034617_cached_res")) # 180

## Double funds per vintage (40), same vintages
run <- "cache_q_factors_20250808_231619_simulated_cashflows_EW_VYP/"
h$same.vin.double.funds <- bivar(paste0(run, "2025-08-13_220154_cached_res")) # 180

## Half vintages, same funds per vintage (20)
run <- "cache_q_factors_20250808_233803_simulated_cashflows_EW_VYP/"
h$vin8695.same.funds <- bivar(paste0(run, "2025-08-14_115400_cached_res")) # 180, 1986-1995
run <- "cache_q_factors_20250808_234201_simulated_cashflows_EW_VYP/"
h$vin9605.same.funds <- bivar(paste0(run, "2025-08-14_122430_cached_res")) # 180, 1996-2005

## Double Half summary (500 iter) ----

df.time <- data.frame(
  Base = c(1986, 2005, 20, h$base$MKT$MKT$Mean, h$base$MKT$MKT$SD),
  Big.n_v = c(1986, 2005, 40, h$same.vin.double.funds$MKT$MKT$Mean, h$same.vin.double.funds$MKT$MKT$SD),
  Big.V = c(1967, 2005, 10, h$double.vin.half.funds$MKT$MKT$Mean, h$double.vin.half.funds$MKT$MKT$SD),
  Big.V = c(1967, 2005, 20, h$double.vin.same.funds$MKT$MKT$Mean, h$double.vin.same.funds$MKT$MKT$SD),
  Small.V = c(1986, 1995 , 20, h$vin8695.same.funds$MKT$MKT$Mean, h$vin8695.same.funds$MKT$MKT$SD),
  Small.V = c(1996, 2005, 20, h$vin9605.same.funds$MKT$MKT$Mean, h$vin9605.same.funds$MKT$MKT$SD)
  )
rownames(df.time) <- c("Start vintage", "End vintage", "#Funds per vintage", "Mean MKT estimate", "SD MKT estimate")

print(xtable::xtable(df.time, digits = 3, caption = "Simulation study for varying  number of vintages and number of funds per vintage with true beta = 1, maximum month 180, and 500 simulation iterations.", 
                     label = "tab:simulation_study_size"), include.rownames=TRUE)



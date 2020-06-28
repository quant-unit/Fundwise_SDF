## analyze simulation ----
if(sys.nframe() == 0L) rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

bivar <- function(files, digi=3) {
  l <- list()
  for (file in files) {
    name <- paste0("data_out/cache_simulation_1000/", file, ".csv")
    l[[file]] <- read.csv(name)
  }
  df0 <- data.frame(do.call(rbind, l))
  
  print(nrow(df0))
  print(paste("max.month:", df0$max.month[1]))
  
  res <- list()
  for (factor in levels(df0$Factor)) {
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
l$sl_1$one <- bivar("2020-06-21_211429_cached_res") # 1
l$sl_60$one <- bivar("2020-06-21_214048_cached_res") # 60
l$sl_120$one <- bivar("2020-06-21_220834_cached_res") # 120
l$sl_180$one <- bivar("2020-06-21_223600_cached_res") # 180
l$sl_240$one <- bivar("2020-06-21_230244_cached_res") # 240
l$sl_300$one <- bivar("2020-06-22_050528_cached_res") # 300
l$sl_360$one <- bivar("2020-06-22_053047_cached_res") # 360

# TODO: delete - duplciates
bivar("2020-06-26_135008_cached_res")
bivar("2020-06-26_141646_cached_res")
bivar("2020-06-26_144253_cached_res")

# individual funds (200 iterations)
l$sl_1$funds <- bivar(c("2020-06-19_145053_cached_res", "2020-06-19_202353_cached_res", "2020-06-20_023931_cached_res", "2020-06-20_170146_cached_res")) # 1
l$sl_60$funds <- bivar(c("2020-06-19_153628_cached_res", "2020-06-19_211021_cached_res", "2020-06-20_034103_cached_res", "2020-06-20_180831_cached_res")) # 60
l$sl_120$funds <- bivar(c("2020-06-19_162009_cached_res", "2020-06-19_215546_cached_res", "2020-06-20_045252_cached_res", "2020-06-20_190809_cached_res")) # 120
l$sl_180$funds <- bivar(c("2020-06-19_170444_cached_res", "2020-06-19_223843_cached_res", "2020-06-20_055215_cached_res", "2020-06-20_200707_cached_res")) # 180
l$sl_240$funds <- bivar(c("2020-06-19_182159_cached_res", "2020-06-19_232145_cached_res", "2020-06-20_065156_cached_res", "2020-06-20_210719_cached_res")) # 240

## Shorter fund lifetime max(5 + 5) ----
l$sl_1$short <- bivar("2020-06-22_143545_cached_res") # 1
l$sl_60$short <- bivar("2020-06-22_150254_cached_res") # 60
l$sl_120$short <- bivar("2020-06-22_153039_cached_res") # 120 min for new fundliftime
l$sl_180$short <- bivar("2020-06-22_155824_cached_res") # 180
l$sl_240$short <- bivar("2020-06-22_162550_cached_res") # 240
l$sl_300$short <- bivar("2020-06-22_165212_cached_res") # 300
l$sl_360$short <- bivar("2020-06-22_171854_cached_res") # 360

## High Beta, Negative Alpha ------
# SL (Beta 2.5, Alpha -0.0025) # -3% p.a. Alpha
l$sl_1$high <- bivar("2020-06-24_174114_cached_res") # 1
l$sl_60$high <- bivar("2020-06-24_175744_cached_res") # 60
l$sl_120$high <- bivar("2020-06-24_181407_cached_res") # 120
l$sl_180$high <- bivar("2020-06-24_183037_cached_res") # 180
l$sl_240$high <- bivar("2020-06-27_022131_cached_res") # 240
l$sl_300$high <- bivar("2020-06-27_024503_cached_res") # 300
l$sl_360$high <- bivar("2020-06-27_030915_cached_res") # 360


## Exp aff models ----
# EA, true DGP (Beta 1.0, Alpha 0)
l$ea_1$one <- bivar("2020-06-22_230535_cached_res") # 1
l$ea_60$one <- bivar("2020-06-22_235514_cached_res") # 60
l$ea_120$one <- bivar("2020-06-23_004538_cached_res") # 120
l$ea_180$one <- bivar("2020-06-23_013620_cached_res") # 180
l$ea_240$one <- bivar("2020-06-23_022506_cached_res") # 240
l$ea_300$one <- bivar("2020-06-23_031302_cached_res") # 300
l$ea_360$one <- bivar("2020-06-23_035954_cached_res") # 360

# EA (Beta 2.5, Alpha -0.0025) # -3% p.a. Alpha
l$ea_1$high <- bivar("2020-06-25_102832_cached_res") # 1
l$ea_60$high <- bivar("2020-06-25_110238_cached_res") # 60
l$ea_120$high <- bivar("2020-06-25_113642_cached_res") # 120
l$ea_180$high <- bivar("2020-06-25_121521_cached_res") # 180
l$ea_240$high <- bivar("2020-06-25_153206_cached_res") # 240
l$ea_300$high <- bivar("2020-06-25_161423_cached_res") # 300
l$ea_360$high <- bivar("2020-06-25_165925_cached_res") # 360

## Sim Lin & Exp aff summary ------
df <- data.frame(MaxMonth = c(1, 1, 60, 60, 120, 120, 180, 180, 240, 240, 300, 300, 360, 360))
df$MaxMonth <- paste(df$MaxMonth, "-", unlist(lapply(1:(nrow(df)/2), function(x) return(c("mean", "stdv")))))
m <- "sl"
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

## Double half models (500 iter) ----
h <- list()
## Same, same with 500 iters (base case)
h$base <- bivar("2020-06-23_214836_cached_res")

## Double vintages, half funds per vintage (10)
h$double.vin.half.funds <- bivar("2020-06-23_134833_cached_res") # 180

## Double vintages, same funds per vintage (20)
h$double.vin.same.funds <- bivar("2020-06-23_154337_cached_res") # 180

## Double funds per vintage (40), same vintages
h$same.vin.double.funds <- bivar("2020-06-23_191734_cached_res")

## Half vintages, same funds per vintage (20)
h$vin8695.same.funds <- bivar("2020-06-23_201744_cached_res") # 180, 1986-1995
h$vin9605.same.funds <- bivar("2020-06-23_211006_cached_res") # 180, 1996-2005

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



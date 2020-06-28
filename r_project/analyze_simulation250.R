## analyze simulation ----
if(sys.nframe() == 0L) rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

bivar <- function(file, digi=3) {
  name <- paste0("data_out/cache_q_factors/", file, ".csv")
  df0 <- read.csv(name)
  print(paste("max.month:", df0$max.month[1]))
  
  res <- list()
  for (factor in levels(df0$Factor)) {
    coefs <- "MKT"
    if(factor != "MKT") coefs <- c(coefs, factor)
    print(coefs)
    
    df <- df0[df0$Factor == factor, ]
    for(coef in coefs) {
      plot(hist(df[, coef]), main = paste(factor, coef))
      m <- mean(df[, coef])
      s <- sd(df[, coef])
      res[[factor]][[paste0(coef)]] <- list(Mean = m, SD = s)
      y <- paste0(factor, "-", coef, ":  ", round(m, digi), " (", round(s, digi), ")      ")
      print(y)
    }
  }
  invisible(res)
}

## monthly Beta 1 -----
bivar("2020-06-15_190609_cached_res")
bivar("2020-06-15_191243_cached_res")
bivar("2020-06-15_191933_cached_res")
bivar("2020-06-15_192609_cached_res")
bivar("2020-06-15_193214_cached_res")
bivar("2020-06-15_193836_cached_res")
bivar("2020-06-15_194421_cached_res")
bivar("2020-06-15_195032_cached_res")
bivar("2020-06-15_195626_cached_res")
bivar("2020-06-15_200219_cached_res")
bivar("2020-06-15_200812_cached_res")


## cross-sectional unit -----
df.vyp <- read.csv("data_out/cache_q_factors/2020-06-11_131155_cached_res.csv") # MQ 80
df.ifp.200 <- read.csv("data_out/cache_q_factors/2020-06-10_172157_cached_res.csv") # MQ 200
df.ifp <- read.csv("data_out/cache_q_factors/2020-06-10_142319_cached_res.csv") # MQ 80

bivar(df.vyp)
bivar(df.ifp.200)
bivar(df.ifp)

t.test(df.vyp$MKT, df.ifp$MKT)

## Max quarter ----
df.1 <- read.csv("data_out/cache_q_factors/2020-06-11_124900_cached_res.csv")
df.2 <- read.csv("data_out/cache_q_factors/2020-06-11_125253_cached_res.csv")
df.5 <- read.csv("data_out/cache_q_factors/2020-06-11_125744_cached_res.csv")
df.20 <- read.csv("data_out/cache_q_factors/2020-06-11_130220_cached_res.csv")
df.40 <- read.csv("data_out/cache_q_factors/2020-06-11_130655_cached_res.csv")
df.80 <- read.csv("data_out/cache_q_factors/2020-06-11_131155_cached_res.csv")
df.120 <- read.csv("data_out/cache_q_factors/2020-06-11_131625_cached_res.csv")
df.160 <- read.csv("data_out/cache_q_factors/2020-06-11_132049_cached_res.csv")
df.200 <- read.csv("data_out/cache_q_factors/2020-06-11_132532_cached_res.csv")
df.240 <- read.csv("data_out/cache_q_factors/2020-06-11_133004_cached_res.csv")

bivar(df.1)
bivar(df.2)
bivar(df.5)
bivar(df.20)
bivar(df.40)
bivar(df.80)
bivar(df.120)
bivar(df.160)
bivar(df.200)
bivar(df.240)

t.test(df.200$MKT, df.1$MKT)
t.test(df.200$MKT, df.80$MKT)

## Sim lin models ----
## Beta 0.5
df.05.1 <- read.csv("data_out/cache_q_factors/2020-06-11_201421_cached_res.csv") # MQ 1
df.05.1 <- df.05.1[df.05.1$Factor == "MKT", ]
df.05.80 <- read.csv("data_out/cache_q_factors/2020-06-11_202039_cached_res.csv") # MQ 80
df.05.80 <- df.05.80[df.05.80$Factor == "MKT", ]
df.05.200 <- read.csv("data_out/cache_q_factors/2020-06-11_202712_cached_res.csv") # MQ 200
df.05.200 <- df.05.200[df.05.200$Factor == "MKT", ]

bivar(df.05.1)
bivar(df.05.80)
bivar(df.05.200)
## Beta 2.0
df.2.1 <- read.csv("data_out/cache_q_factors/2020-06-11_204915_cached_res.csv") # MQ 1
df.2.1 <- df.2.1[df.2.1$Factor == "MKT", ]
df.2.80 <- read.csv("data_out/cache_q_factors/2020-06-11_205312_cached_res.csv") # MQ 80
df.2.80 <- df.2.80[df.2.80$Factor == "MKT", ]
df.2.200 <- read.csv("data_out/cache_q_factors/2020-06-11_205711_cached_res.csv") # MQ 200
df.2.200 <- df.2.200[df.2.200$Factor == "MKT", ]

bivar(df.2.1)
bivar(df.2.80)
bivar(df.2.200)
## Beta 0.5 & Alpha +0.001
df.050001.1 <- read.csv("data_out/cache_q_factors/2020-06-11_212218_cached_res.csv") # MQ 1
df.050001.1 <- df.050001.1[df.050001.1$Factor == "Alpha", ]
df.050001.80 <- read.csv("data_out/cache_q_factors/2020-06-11_213244_cached_res.csv") # MQ 80
df.050001.80 <- df.050001.80[df.050001.80$Factor == "Alpha", ]
df.050001.200 <- read.csv("data_out/cache_q_factors/2020-06-11_214307_cached_res.csv") # MQ 200
df.050001.200 <- df.050001.200[df.050001.200$Factor == "Alpha", ]

bivar(df.050001.1)
bivar(df.050001.80)
bivar(df.050001.200)
## Beta 2.0 & Alpha -0.001
df.20001.1 <- read.csv("data_out/cache_q_factors/2020-06-11_225902_cached_res.csv") # MQ 1
df.20001.1 <- df.20001.1[df.20001.1$Factor == "Alpha", ]
df.20001.80 <- read.csv("data_out/cache_q_factors/2020-06-11_230521_cached_res.csv") # MQ 80
df.20001.80 <- df.20001.80[df.20001.80$Factor == "Alpha", ]
df.20001.200 <- read.csv("data_out/cache_q_factors/2020-06-11_231128_cached_res.csv") # MQ 200
df.20001.200 <- df.20001.200[df.20001.200$Factor == "Alpha", ]

bivar(df.20001.1)
bivar(df.20001.80)
bivar(df.20001.200)
## Sim lin summary ----
df <- data.frame(MaxQ = c("mean", "stdv", "mean", "stdv", "mean", "stdv"))
df$MaxQ <- paste(c(1, 1, 80, 80, 200, 200), df$MaxQ, sep = " - ")
# No Alpha
df$MKT0.5 <- c(mean(df.05.1$MKT), sd(df.05.1$MKT), mean(df.05.80$MKT), sd(df.05.80$MKT), mean(df.05.200$MKT), sd(df.05.200$MKT))
df$MKT2 <- c(mean(df.2.1$MKT), sd(df.2.1$MKT), mean(df.2.80$MKT), sd(df.2.80$MKT), mean(df.2.200$MKT), sd(df.2.200$MKT))
# Beta 0.5, Alpha = 0.001
df$MKT05_alpha <- c(mean(df.050001.1$MKT), sd(df.050001.1$MKT), mean(df.050001.80$MKT), sd(df.050001.80$MKT), mean(df.050001.200$MKT), sd(df.050001.200$MKT))
df$Alpha0.001 <- c(mean(df.050001.1$Alpha), sd(df.050001.1$Alpha), mean(df.050001.80$Alpha), sd(df.050001.80$Alpha), mean(df.050001.200$Alpha), sd(df.050001.200$Alpha))
df$Alpha0.001 <- df$Alpha0.001  * 100
# Beta 2.0, Alpha = -0.001
df$MKT_alpha <- c(mean(df.20001.1$MKT), sd(df.20001.1$MKT), mean(df.20001.80$MKT), sd(df.20001.80$MKT), mean(df.20001.200$MKT), sd(df.20001.200$MKT))
df$Alpha_0.001 <- c(mean(df.20001.1$Alpha), sd(df.20001.1$Alpha), mean(df.20001.80$Alpha), sd(df.20001.80$Alpha), mean(df.20001.200$Alpha), sd(df.20001.200$Alpha))
df$Alpha_0.001 <- df$Alpha_0.001 * 100
print(xtable::xtable(df, digits = 3, caption = "Simulation study (sim.lin) for maximum quarter.", label = "tab:simulation_study_max_q"), include.rownames=FALSE)

## Exp aff models ----
# EA, true DGP (Beta 1.0)
df.exp.aff.1.80 <- read.csv("data_out/cache_q_factors/2020-06-12_173512_cached_res.csv") # MQ 80
df.exp.aff.1.200 <- read.csv("data_out/cache_q_factors/2020-06-12_175951_cached_res.csv") # MQ 200
bivar(df.exp.aff.1.80)
bivar(df.exp.aff.1.200)

# EA, true DGP (Beta 2.0)
df.exp.aff.2.1 <- read.csv("data_out/cache_q_factors/2020-06-12_184158_cached_res.csv") # MQ 1
df.exp.aff.2.80 <- read.csv("data_out/cache_q_factors/2020-06-12_185821_cached_res.csv") # MQ 80
df.exp.aff.2.200 <- read.csv("data_out/cache_q_factors/2020-06-12_191358_cached_res.csv") # MQ 200
ea.2.1 <- bivar(df.exp.aff.2.1)
ea.2.80 <- bivar(df.exp.aff.2.80)
ea.2.200 <- bivar(df.exp.aff.2.200)

# EA (Beta 0.5)
df.exp.aff.05.1 <- read.csv("data_out/cache_q_factors/2020-06-12_204714_cached_res.csv") # MQ 1
df.exp.aff.05.80 <- read.csv("data_out/cache_q_factors/2020-06-12_211638_cached_res.csv") # MQ 80
df.exp.aff.05.200 <- read.csv("data_out/cache_q_factors/2020-06-12_214558_cached_res.csv") # MQ 200
ea.05.1 <- bivar(df.exp.aff.05.1)
ea.05.80 <- bivar(df.exp.aff.05.80)
ea.05.200 <- bivar(df.exp.aff.05.200)

# EA (Beta 2.0, Alpha -0.001)
df.exp.aff.2na.1 <- read.csv("data_out/cache_q_factors/2020-06-13_003008_cached_res.csv") # MQ 1
df.exp.aff.2na.80 <- read.csv("data_out/cache_q_factors/2020-06-13_003658_cached_res.csv") # MQ 80
df.exp.aff.2na.200 <- read.csv("data_out/cache_q_factors/2020-06-13_004339_cached_res.csv") # MQ 200
ea.2na.1 <- bivar(df.exp.aff.2na.1)
ea.2na.80 <- bivar(df.exp.aff.2na.80)
ea.2na.200 <- bivar(df.exp.aff.2na.200)

# EA (Beta 0.5, Alpha 0.001)
df.exp.aff.05pa.1 <- read.csv("data_out/cache_q_factors/2020-06-13_102626_cached_res.csv") # MQ 1
df.exp.aff.05pa.80 <- read.csv("data_out/cache_q_factors/2020-06-13_103723_cached_res.csv") # MQ 80
df.exp.aff.05pa.200 <- read.csv("data_out/cache_q_factors/2020-06-13_104838_cached_res.csv") # MQ 200
ea.05pa.1 <- bivar(df.exp.aff.05pa.1)
ea.05pa.80 <- bivar(df.exp.aff.05pa.80)
ea.05pa.200 <- bivar(df.exp.aff.05pa.200)


## Exp aff summary ------
df <- data.frame(MaxQ = c("mean", "stdv", "mean", "stdv", "mean", "stdv"))
df$MaxQ <- paste(c(1, 1, 80, 80, 200, 200), df$MaxQ, sep = " - ")
# No Alpha
df$MKT0.5 <- c(ea.05.1$MKT$MKT$Mean, ea.05.1$MKT$MKT$SD, 
               ea.05.80$MKT$MKT$Mean, ea.05.80$MKT$MKT$SD, 
               ea.05.200$MKT$MKT$Mean, ea.05.200$MKT$MKT$SD)
df$MKT2 <- c(ea.2.1$MKT$MKT$Mean, ea.2.1$MKT$MKT$SD,
             ea.2.80$MKT$MKT$Mean, ea.2.80$MKT$MKT$SD,
             ea.2.200$MKT$MKT$Mean, ea.2.200$MKT$MKT$SD)
# Beta 0.5, Alpha = 0.001
df$MKT05_alpha <- c(ea.05pa.1$Alpha$MKT$Mean, ea.05pa.1$Alpha$MKT$SD, 
                    ea.05pa.80$Alpha$MKT$Mean, ea.05pa.80$Alpha$MKT$SD, 
                    ea.05pa.200$Alpha$MKT$Mean, ea.05pa.200$Alpha$MKT$SD)
df$Alpha0.001 <- c(ea.05pa.1$Alpha$Alpha$Mean, ea.05pa.1$Alpha$Alpha$SD, 
                   ea.05pa.80$Alpha$Alpha$Mean, ea.05pa.80$Alpha$Alpha$SD, 
                   ea.05pa.200$Alpha$Alpha$Mean, ea.05pa.200$Alpha$Alpha$SD)
df$Alpha0.001 <- df$Alpha0.001  * 100
# Beta 2.0, Alpha = -0.001
df$MKT_alpha <- c(ea.2na.1$Alpha$MKT$Mean, ea.2na.1$Alpha$MKT$SD, 
                  ea.2na.80$Alpha$MKT$Mean, ea.2na.80$Alpha$MKT$SD, 
                  ea.2na.200$Alpha$MKT$Mean, ea.2na.200$Alpha$MKT$SD)
df$Alpha_0.001 <- c(ea.2na.1$Alpha$Alpha$Mean, ea.2na.1$Alpha$Alpha$SD, 
                    ea.2na.80$Alpha$Alpha$Mean, ea.2na.80$Alpha$Alpha$SD, 
                    ea.2na.200$Alpha$Alpha$Mean, ea.2na.200$Alpha$Alpha$SD)
df$Alpha_0.001 <- df$Alpha_0.001 * 100

print(xtable::xtable(df, digits = 3, caption = "Simulation study (exp.aff) for maximum quarter.", label = "tab:simulation_study_max_q"), include.rownames=FALSE)

## Double half models ----
## Double vintages, half funds per vintage (10)
df.double.half <- read.csv("data_out/cache_q_factors/2020-06-11_173144_cached_res.csv") # MQ 80

bivar(df.vyp)
bivar(df.double.half)
t.test(df.vyp$MKT, df.double.half$MKT)

## Double vintages, same funds per vintage (20)
df.double.same <- read.csv("data_out/cache_q_factors/2020-06-11_175653_cached_res.csv") # MQ 80

bivar(df.vyp)
bivar(df.double.same)
t.test(df.vyp$MKT, df.double.same$MKT)
## Double funds per vintage (40), same vintages
df.same.double <- read.csv("data_out/cache_q_factors/2020-06-11_190603_cached_res.csv") # MQ 80

bivar(df.vyp)
bivar(df.same.double)
t.test(df.vyp$MKT, df.same.double$MKT)
## Half vintages, same funds per vintage (20)
df.half.same.8695 <- read.csv("data_out/cache_q_factors/2020-06-11_194248_cached_res.csv") # MQ 80
df.half.same.9605 <- read.csv("data_out/cache_q_factors/2020-06-11_195131_cached_res.csv") # MQ 80

bivar(df.vyp)
bivar(df.half.same.8695)
bivar(df.half.same.9605)
t.test(df.half.same.8695$MKT, df.half.same.9605$MKT)
## Double Half summary ----

df.time <- data.frame(
  Base = c(1986, 2005, 20, mean(df.vyp$MKT), sd(df.vyp$MKT)),
  Big.n_v = c(1986, 2005, 40, mean(df.same.double$MKT), sd(df.same.double$MKT)),
  Big.V = c(1967, 2005, 10, mean(df.double.half$MKT), sd(df.double.half$MKT)),
  Big.V = c(1967, 2005, 20, mean(df.double.same$MKT), sd(df.double.same$MKT)),
  Small.V = c(1986, 1995 , 20, mean(df.half.same.8695$MKT), sd(df.half.same.8695$MKT)),
  Small.V = c(1996, 2005, 20, mean(df.half.same.9605$MKT), sd(df.half.same.9605$MKT))
  )
rownames(df.time) <- c("Start vintage", "End vintage", "#Funds per vintage", "Mean MKT estimate", "SD MKT estimate")

print(xtable::xtable(df.time, digits = 3, caption = "Simulation study for data size by different vintage years and funds per vintage.", 
                     label = "tab:simulation_study_size"), include.rownames=TRUE)


## shifted lognormal & large stdv 108% p.a. ----
# 20 Funds (stdv 0.312) linear
df.1 <- read.csv("data_out/cache_q_factors/2020-06-13_132853_cached_res.csv")
df.80 <- read.csv("data_out/cache_q_factors/2020-06-13_133323_cached_res.csv")
df.200 <- read.csv("data_out/cache_q_factors/2020-06-13_133800_cached_res.csv")
bivar(df.200)

# 20 Funds (stdv 0.312) exp.aff
df.1.ea <- read.csv("data_out/cache_q_factors/2020-06-13_144952_cached_res.csv")
df.80.ea <- read.csv("data_out/cache_q_factors/2020-06-13_150027_cached_res.csv")
df.200.ea <- read.csv("data_out/cache_q_factors/2020-06-13_151025_cached_res.csv")
bivar(df.200.ea)

# 50 Funds (stdv 0.312) linear
df.1 <- read.csv("data_out/cache_q_factors/2020-06-13_153411_cached_res.csv")
df.80 <- read.csv("data_out/cache_q_factors/2020-06-13_153832_cached_res.csv")
df.200 <- read.csv("data_out/cache_q_factors/2020-06-13_154258_cached_res.csv")
bivar(df.200)

# 50 Funds (stdv 0.312) lognormal error - 1
df.1 <- read.csv("data_out/cache_q_factors/2020-06-13_175154_cached_res.csv")
df.80 <- read.csv("data_out/cache_q_factors/2020-06-13_175650_cached_res.csv")
df.200 <- read.csv("data_out/cache_q_factors/2020-06-13_180228_cached_res.csv")
bivar(df.200)

# 50 Funds (stdv 0.312) shifted lognormal error - 0.75
df.80 <- read.csv("data_out/cache_q_factors/2020-06-14_134402_cached_res.csv")
df.200 <- read.csv("data_out/cache_q_factors/2020-06-14_134829_cached_res.csv")
bivar(df.200)

# 20 Funds (stdv 0.2) shifted lognormal error - 0.75
df.1 <- read.csv("data_out/cache_q_factors/2020-06-14_142237_cached_res.csv")
df.80 <- read.csv("data_out/cache_q_factors/2020-06-14_142822_cached_res.csv")
df.200 <- read.csv("data_out/cache_q_factors/2020-06-14_143315_cached_res.csv")
bivar(df.200)


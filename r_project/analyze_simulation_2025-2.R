## analyze simulation ----
if(sys.nframe() == 0L) rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

bivar <- function(files, digi=3) {
  l <- list()
  for (file in files) {
    name <- paste0("data_out_2025-2/", file, ".csv")
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
l$sl_1$one <-   bivar(paste0(run, "2025-08-17_191933_cached_res")) # 1
l$sl_60$one <-  bivar(paste0(run, "2025-08-17_192506_cached_res")) # 60
l$sl_120$one <- bivar(paste0(run, "2025-08-17_193052_cached_res")) # 120
l$sl_150$one <- bivar(paste0(run, "2025-08-17_193646_cached_res")) # 150
l$sl_180$one <- bivar(paste0(run, "2025-08-17_194249_cached_res")) # 180
l$sl_210$one <- bivar(paste0(run, "2025-08-17_194859_cached_res")) # 210
l$sl_240$one <- bivar(paste0(run, "2025-08-17_195514_cached_res")) # 240
l$sl_300$one <- bivar(paste0(run, "2025-08-17_200140_cached_res")) # 300
l$sl_360$one <- bivar(paste0(run, "2025-08-17_200816_cached_res")) # 360

# Single duration date instead of date averaging
l$sl_duration$one <- bivar(paste0(run, "2025-08-18_220836_cached_res")) # duration
l$sl_duration$one <- bivar(paste0(run, "2025-08-18_222238_cached_res")) # duration (to assure max.month does not matter)

# TODO: Single date
l$sl_single_date1$one <- bivar(paste0(run, "2025-08-18_224944_cached_res")) # single date
l$sl_single_date10$one <- bivar(paste0(run, "2025-08-18_225958_cached_res")) # single date
l$sl_single_date20$one <- bivar(paste0(run, "2025-08-18_231021_cached_res")) # single date
l$sl_single_date30$one <- bivar(paste0(run, "2025-08-18_232052_cached_res")) # single date
l$sl_single_date40$one <- bivar(paste0(run, "2025-08-18_233133_cached_res")) # single date
l$sl_single_date50$one <- bivar(paste0(run, "2025-08-18_234225_cached_res")) # single date
l$sl_single_date60$one <- bivar(paste0(run, "2025-08-18_235333_cached_res")) # single date
l$sl_single_date70$one <- bivar(paste0(run, "2025-08-19_000441_cached_res")) # single date
l$sl_single_date80$one <- bivar(paste0(run, "2025-08-19_001603_cached_res")) # single date
l$sl_single_date90$one <- bivar(paste0(run, "2025-08-19_002742_cached_res")) # single date
l$sl_single_date100$one <- bivar(paste0(run, "2025-08-19_003929_cached_res")) # single date
l$sl_single_date110$one <- bivar(paste0(run, "2025-08-19_005119_cached_res")) # single date
l$sl_single_date120$one <- bivar(paste0(run, "2025-08-19_010323_cached_res")) # single date
l$sl_single_date150$one <- bivar(paste0(run, "2025-08-19_011608_cached_res")) # single date
l$sl_single_date180$one <- bivar(paste0(run, "2025-08-19_012913_cached_res")) # single date
l$sl_single_date210$one <- bivar(paste0(run, "2025-08-19_014238_cached_res")) # single date
l$sl_single_date240$one <- bivar(paste0(run, "2025-08-19_015643_cached_res")) # single date
l$sl_single_date300$one <- bivar(paste0(run, "2025-08-19_021213_cached_res")) # single date
l$sl_single_date360$one <- bivar(paste0(run, "2025-08-19_022921_cached_res")) # single date

# individual funds (new in 2025: also 1000 iterations)
run <- "cache_q_factors_20250808_222540_simulated_cashflows_EW/"
l$sl_1$funds <- bivar(
  c(
    paste0(run, "2025-08-17_201954_cached_res"),
    paste0(run, "2025-08-17_214212_cached_res"),
    paste0(run, "2025-08-17_230517_cached_res"),
    paste0(run, "2025-08-18_002926_cached_res"),
    paste0(run, "2025-08-18_015524_cached_res"),
    paste0(run, "2025-08-18_032000_cached_res"),
    paste0(run, "2025-08-18_045107_cached_res"),
    paste0(run, "2025-08-18_062123_cached_res"),
    paste0(run, "2025-08-18_075258_cached_res"),
    paste0(run, "2025-08-18_092457_cached_res")
                        )
  ) # 1
l$sl_60$funds <- bivar(
  c(
    paste0(run, "2025-08-17_202631_cached_res"),
    paste0(run, "2025-08-17_214902_cached_res"),
    paste0(run, "2025-08-17_231209_cached_res"),
    paste0(run, "2025-08-18_003621_cached_res"),
    paste0(run, "2025-08-18_020219_cached_res"),
    paste0(run, "2025-08-18_032728_cached_res"),
    paste0(run, "2025-08-18_045818_cached_res"),
    paste0(run, "2025-08-18_062852_cached_res"),
    paste0(run, "2025-08-18_080032_cached_res"),
    paste0(run, "2025-08-18_093233_cached_res")
  )
) # 60
l$sl_120$funds <- bivar(
  c(
    paste0(run, "2025-08-17_203321_cached_res"),
    paste0(run, "2025-08-17_215615_cached_res"),
    paste0(run, "2025-08-17_231916_cached_res"),
    paste0(run, "2025-08-18_004333_cached_res"),
    paste0(run, "2025-08-18_020916_cached_res"),
    paste0(run, "2025-08-18_033506_cached_res"),
    paste0(run, "2025-08-18_050605_cached_res"),
    paste0(run, "2025-08-18_063634_cached_res"),
    paste0(run, "2025-08-18_080808_cached_res"),
    paste0(run, "2025-08-18_094035_cached_res")
  )
) # 120
l$sl_150$funds <- bivar(
  c(
    paste0(run, "2025-08-17_204110_cached_res"),
    paste0(run, "2025-08-17_220413_cached_res"),
    paste0(run, "2025-08-17_232722_cached_res"),
    paste0(run, "2025-08-18_005147_cached_res"),
    paste0(run, "2025-08-18_021703_cached_res"),
    paste0(run, "2025-08-18_034346_cached_res"),
    paste0(run, "2025-08-18_051500_cached_res"),
    paste0(run, "2025-08-18_064528_cached_res"),
    paste0(run, "2025-08-18_081635_cached_res"),
    paste0(run, "2025-08-18_094939_cached_res")
  )
) # 150
l$sl_180$funds <- bivar(
  c(
    paste0(run, "2025-08-17_204950_cached_res"),
    paste0(run, "2025-08-17_221259_cached_res"),
    paste0(run, "2025-08-17_233621_cached_res"),
    paste0(run, "2025-08-18_010107_cached_res"),
    paste0(run, "2025-08-18_022604_cached_res"),
    paste0(run, "2025-08-18_035336_cached_res"),
    paste0(run, "2025-08-18_052452_cached_res"),
    paste0(run, "2025-08-18_065519_cached_res"),
    paste0(run, "2025-08-18_082617_cached_res"),
    paste0(run, "2025-08-18_095953_cached_res")
  )
) # 180
l$sl_210$funds <- bivar(
  c(
    paste0(run, "2025-08-17_205902_cached_res"),
    paste0(run, "2025-08-17_222224_cached_res"),
    paste0(run, "2025-08-17_234540_cached_res"),
    paste0(run, "2025-08-18_011057_cached_res"),
    paste0(run, "2025-08-18_023540_cached_res"),
    paste0(run, "2025-08-18_040401_cached_res"),
    paste0(run, "2025-08-18_053512_cached_res"),
    paste0(run, "2025-08-18_070539_cached_res"),
    paste0(run, "2025-08-18_083653_cached_res"),
    paste0(run, "2025-08-18_101041_cached_res")
  )
) # 210
l$sl_240$funds <- bivar(
  c(
    paste0(run, "2025-08-17_210843_cached_res"),
    paste0(run, "2025-08-17_223217_cached_res"),
    paste0(run, "2025-08-17_235530_cached_res"),
    paste0(run, "2025-08-18_012112_cached_res"),
    paste0(run, "2025-08-18_024543_cached_res"),
    paste0(run, "2025-08-18_041459_cached_res"),
    paste0(run, "2025-08-18_054600_cached_res"),
    paste0(run, "2025-08-18_071627_cached_res"),
    paste0(run, "2025-08-18_084754_cached_res"),
    paste0(run, "2025-08-18_102154_cached_res")
  )
) # 240
l$sl_300$funds <- bivar(
  c(
    paste0(run, "2025-08-17_211905_cached_res"),
    paste0(run, "2025-08-17_224251_cached_res"),
    paste0(run, "2025-08-18_000611_cached_res"),
    paste0(run, "2025-08-18_013207_cached_res"),
    paste0(run, "2025-08-18_025621_cached_res"),
    paste0(run, "2025-08-18_042639_cached_res"),
    paste0(run, "2025-08-18_055716_cached_res"),
    paste0(run, "2025-08-18_072803_cached_res"),
    paste0(run, "2025-08-18_085940_cached_res"),
    paste0(run, "2025-08-18_103356_cached_res")
  )
) # 300
l$sl_360$funds <- bivar(
  c(
    paste0(run, "2025-08-17_213002_cached_res"),
    paste0(run, "2025-08-17_225403_cached_res"),
    paste0(run, "2025-08-18_001729_cached_res"),
    paste0(run, "2025-08-18_014339_cached_res"),
    paste0(run, "2025-08-18_030734_cached_res"),
    paste0(run, "2025-08-18_043852_cached_res"),
    paste0(run, "2025-08-18_060900_cached_res"),
    paste0(run, "2025-08-18_074008_cached_res"),
    paste0(run, "2025-08-18_091208_cached_res"),
    paste0(run, "2025-08-18_104635_cached_res")
  )
) # 360

# Aggregate for plotting
months <- c(1, 60, 120, 150, 180, 210, 240, 300, 360)
prefix1 <- "sl_"
prefix2 <- "sl_"
unit1 <- "one"
unit2 <- "funds"
legend.text <- c("Mean Funds", "Mean VYP", "Stdv Funds", "Stdv VYP")
xlab.text <- "Max Months"

if (FALSE) {
  prefix2 <- "sl_single_date"
  unit2 <-  "one"
  legend.text <- c("Mean VYP: Single Date", "Mean VYP: Average Dates", 
                   "Stdv VYP: Single Date", "Stdv VYP: Average Datess")
  xlab.text <- "Single Date / Max Months"
}


df1 <- data.frame(
  months = months,
  VYP_mean = sapply(months, function(x) (l[[paste0(prefix1, x)]][[unit1]]$MKT$MKT$Mean)),
  VYP_sd = sapply(months, function(x) (l[[paste0(prefix1, x)]][[unit1]]$MKT$MKT$SD)),
  Funds_mean = sapply(months, function(x) (l[[paste0(prefix2, x)]][[unit2]]$MKT$MKT$Mean)),
  Funds_sd = sapply(months, function(x) (l[[paste0(prefix2, x)]][[unit2]]$MKT$MKT$SD))
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
     xlab = xlab.text, ylab = "Mean and Stdv of Estimates")
lines(df1$months, df1$Funds_mean, type = "b", col = "blue")
lines(df1$months, df1$VYP_sd, type = "b", lty = 3)
lines(df1$months, df1$Funds_sd, type = "b", lty = 3, col = "blue")
abline(h=c(0, 1), col="grey")
legend("topright", bty = "n", 
       legend = legend.text, 
       col = c("blue", "black"), lty = c(1,1,3,3), cex = 0.8)

if(do.eps){ 
  dev.off() 
}

## Shorter fund lifetime max(investing=2 + max.holding=4) ----

run <- "cache_q_factors_20250819_221615_simulated_cashflows_EW_VYP/"
l$sl_1$short_i2_h4 <-   bivar(paste0(run, "2025-08-20_053500_cached_res")) # 1
l$sl_30$short_i2_h4 <-  bivar(paste0(run, "2025-08-20_054522_cached_res")) # 30
l$sl_60$short_i2_h4 <-  bivar(paste0(run, "2025-08-20_055554_cached_res")) # 60
l$sl_90$short_i2_h4 <-  bivar(paste0(run, "2025-08-20_060642_cached_res")) # 90
l$sl_120$short_i2_h4 <- bivar(paste0(run, "2025-08-20_061757_cached_res")) # 120
l$sl_150$short_i2_h4 <- bivar(paste0(run, "2025-08-20_062941_cached_res")) # 150
l$sl_180$short_i2_h4 <- bivar(paste0(run, "2025-08-20_064157_cached_res")) # 180
l$sl_210$short_i2_h4 <- bivar(paste0(run, "2025-08-20_065429_cached_res")) # 210
l$sl_240$short_i2_h4 <- bivar(paste0(run, "2025-08-20_070716_cached_res")) # 240
l$sl_300$short_i2_h4 <- bivar(paste0(run, "2025-08-20_072016_cached_res")) # 300
l$sl_360$short_i2_h4 <- bivar(paste0(run, "2025-08-20_073400_cached_res")) # 360

## Shorter fund lifetime max(investing=4 + max.holding=2) ----

run <- "cache_q_factors_20250819_222708_simulated_cashflows_EW_VYP/"
l$sl_1$short_i4_h2 <-   bivar(paste0(run, "2025-08-20_074532_cached_res")) # 1
l$sl_30$short_i4_h2 <-  bivar(paste0(run, "2025-08-20_075514_cached_res")) # 30
l$sl_60$short_i4_h2 <-  bivar(paste0(run, "2025-08-20_080504_cached_res")) # 60
l$sl_90$short_i4_h2 <-  bivar(paste0(run, "2025-08-20_081515_cached_res")) # 90
l$sl_120$short_i4_h2 <- bivar(paste0(run, "2025-08-20_082557_cached_res")) # 120
l$sl_150$short_i4_h2 <- bivar(paste0(run, "2025-08-20_083720_cached_res")) # 150
l$sl_180$short_i4_h2 <- bivar(paste0(run, "2025-08-20_084913_cached_res")) # 180
l$sl_210$short_i4_h2 <- bivar(paste0(run, "2025-08-20_090124_cached_res")) # 210
l$sl_240$short_i4_h2 <- bivar(paste0(run, "2025-08-20_091350_cached_res")) # 240
l$sl_300$short_i4_h2 <- bivar(paste0(run, "2025-08-20_092632_cached_res")) # 300
l$sl_360$short_i4_h2 <- bivar(paste0(run, "2025-08-20_094023_cached_res")) # 360

## Shorter fund lifetime max(investing=1 + max.holding=5) ----

run <- "cache_q_factors_20250820_135255_simulated_cashflows_EW_VYP/"
l$sl_1$short_i1_h5 <-   bivar(paste0(run, "2025-08-20_162902_cached_res")) # 1
l$sl_30$short_i1_h5 <-  bivar(paste0(run, "2025-08-20_164006_cached_res")) # 30
l$sl_60$short_i1_h5 <-  bivar(paste0(run, "2025-08-20_165129_cached_res")) # 60
l$sl_90$short_i1_h5 <-  bivar(paste0(run, "2025-08-20_170256_cached_res")) # 90
l$sl_120$short_i1_h5 <- bivar(paste0(run, "2025-08-20_171436_cached_res")) # 120
l$sl_150$short_i1_h5 <- bivar(paste0(run, "2025-08-20_172646_cached_res")) # 150
l$sl_180$short_i1_h5 <- bivar(paste0(run, "2025-08-20_173918_cached_res")) # 180
l$sl_210$short_i1_h5 <- bivar(paste0(run, "2025-08-20_175209_cached_res")) # 210
l$sl_240$short_i1_h5 <- bivar(paste0(run, "2025-08-20_180507_cached_res")) # 240
l$sl_300$short_i1_h5 <- bivar(paste0(run, "2025-08-20_181818_cached_res")) # 300
l$sl_360$short_i1_h5 <- bivar(paste0(run, "2025-08-20_183213_cached_res")) # 360

## Shorter fund lifetime max(investing=5 + max.holding=1) ----

run <- "cache_q_factors_20250820_140352_simulated_cashflows_EW_VYP/"
l$sl_1$short_i5_h1 <-   bivar(paste0(run, "2025-08-20_184337_cached_res")) # 1
l$sl_30$short_i5_h1 <-  bivar(paste0(run, "2025-08-20_185301_cached_res")) # 30
l$sl_60$short_i5_h1 <-  bivar(paste0(run, "2025-08-20_190229_cached_res")) # 60
l$sl_90$short_i5_h1 <-  bivar(paste0(run, "2025-08-20_191225_cached_res")) # 90
l$sl_120$short_i5_h1 <- bivar(paste0(run, "2025-08-20_192304_cached_res")) # 120
l$sl_150$short_i5_h1 <- bivar(paste0(run, "2025-08-20_193431_cached_res")) # 150
l$sl_180$short_i5_h1 <- bivar(paste0(run, "2025-08-20_194656_cached_res")) # 180
l$sl_210$short_i5_h1 <- bivar(paste0(run, "2025-08-20_195944_cached_res")) # 210
l$sl_240$short_i5_h1 <- bivar(paste0(run, "2025-08-20_201251_cached_res")) # 240
l$sl_300$short_i5_h1 <- bivar(paste0(run, "2025-08-20_202625_cached_res")) # 300
l$sl_360$short_i5_h1 <- bivar(paste0(run, "2025-08-20_204054_cached_res")) # 360

## Shorter fund lifetime max(investing=1 + max.holding=1) ----

run <- "cache_q_factors_20250821_153327_simulated_cashflows_EW_VYP/"
l$sl_1$short_i1_h1 <-   bivar(paste0(run, "2025-08-21_210534_cached_res")) # 1
l$sl_30$short_i1_h1 <-  bivar(paste0(run, "2025-08-21_211433_cached_res")) # 30
l$sl_60$short_i1_h1 <-  bivar(paste0(run, "2025-08-21_212347_cached_res")) # 60
l$sl_90$short_i1_h1 <-  bivar(paste0(run, "2025-08-21_213327_cached_res")) # 90
l$sl_120$short_i1_h1 <- bivar(paste0(run, "2025-08-21_214346_cached_res")) # 120
l$sl_150$short_i1_h1 <- bivar(paste0(run, "2025-08-21_215449_cached_res")) # 150
l$sl_180$short_i1_h1 <- bivar(paste0(run, "2025-08-21_220620_cached_res")) # 180
l$sl_210$short_i1_h1 <- bivar(paste0(run, "2025-08-21_221759_cached_res")) # 210
l$sl_240$short_i1_h1 <- bivar(paste0(run, "2025-08-21_222952_cached_res")) # 240
l$sl_300$short_i1_h1 <- bivar(paste0(run, "2025-08-21_224228_cached_res")) # 300
l$sl_360$short_i1_h1 <- bivar(paste0(run, "2025-08-21_225617_cached_res")) # 360

## Shorter fund lifetime max(investing=2 + max.holding=2) ----

run <- "cache_q_factors_20250819_215442_simulated_cashflows_EW_VYP/"
l$sl_1$short_i2_h2 <-   bivar(paste0(run, "2025-08-20_011502_cached_res")) # 1
l$sl_30$short_i2_h2 <-  bivar(paste0(run, "2025-08-20_012502_cached_res")) # 30
l$sl_60$short_i2_h2 <-  bivar(paste0(run, "2025-08-20_013518_cached_res")) # 60
l$sl_90$short_i2_h2 <-  bivar(paste0(run, "2025-08-20_014548_cached_res")) # 90
l$sl_120$short_i2_h2 <- bivar(paste0(run, "2025-08-20_015638_cached_res")) # 120
l$sl_150$short_i2_h2 <- bivar(paste0(run, "2025-08-20_020808_cached_res")) # 150
l$sl_180$short_i2_h2 <- bivar(paste0(run, "2025-08-20_022012_cached_res")) # 180
l$sl_210$short_i2_h2 <- bivar(paste0(run, "2025-08-20_023242_cached_res")) # 210
l$sl_240$short_i2_h2 <- bivar(paste0(run, "2025-08-20_024522_cached_res")) # 240
l$sl_300$short_i2_h2 <- bivar(paste0(run, "2025-08-20_025825_cached_res")) # 300
l$sl_360$short_i2_h2 <- bivar(paste0(run, "2025-08-20_031225_cached_res")) # 360

## Shorter fund lifetime max(investing=3 + max.holding=3) ----

run <- "cache_q_factors_20250820_192412_simulated_cashflows_EW_VYP/"
l$sl_1$short_i3_h3 <-   bivar(paste0(run, "2025-08-20_235911_cached_res")) # 1
l$sl_30$short_i3_h3 <-  bivar(paste0(run, "2025-08-21_000928_cached_res")) # 30
l$sl_60$short_i3_h3 <-  bivar(paste0(run, "2025-08-21_001952_cached_res")) # 60
l$sl_90$short_i3_h3 <-  bivar(paste0(run, "2025-08-21_003037_cached_res")) # 90
l$sl_120$short_i3_h3 <- bivar(paste0(run, "2025-08-21_004147_cached_res")) # 120
l$sl_150$short_i3_h3 <- bivar(paste0(run, "2025-08-21_005332_cached_res")) # 150
l$sl_180$short_i3_h3 <- bivar(paste0(run, "2025-08-21_010553_cached_res")) # 180
l$sl_210$short_i3_h3 <- bivar(paste0(run, "2025-08-21_011830_cached_res")) # 210
l$sl_240$short_i3_h3 <- bivar(paste0(run, "2025-08-21_013122_cached_res")) # 240
l$sl_300$short_i3_h3 <- bivar(paste0(run, "2025-08-21_014427_cached_res")) # 300
l$sl_360$short_i3_h3 <- bivar(paste0(run, "2025-08-21_015822_cached_res")) # 360


## Shorter fund lifetime max(investing=4 + max.holding=4) ----

run <- "cache_q_factors_20250819_220526_simulated_cashflows_EW_VYP/"
l$sl_1$short_i4_h4 <-   bivar(paste0(run, "2025-08-20_032440_cached_res")) # 1
l$sl_30$short_i4_h4 <-  bivar(paste0(run, "2025-08-20_033500_cached_res")) # 30
l$sl_60$short_i4_h4 <-  bivar(paste0(run, "2025-08-20_034524_cached_res")) # 60
l$sl_90$short_i4_h4 <-  bivar(paste0(run, "2025-08-20_035607_cached_res")) # 90
l$sl_120$short_i4_h4 <- bivar(paste0(run, "2025-08-20_040711_cached_res")) # 120
l$sl_150$short_i4_h4 <- bivar(paste0(run, "2025-08-20_041851_cached_res")) # 150
l$sl_180$short_i4_h4 <- bivar(paste0(run, "2025-08-20_043102_cached_res")) # 180
l$sl_210$short_i4_h4 <- bivar(paste0(run, "2025-08-20_044329_cached_res")) # 210
l$sl_240$short_i4_h4 <- bivar(paste0(run, "2025-08-20_045613_cached_res")) # 240
l$sl_300$short_i4_h4 <- bivar(paste0(run, "2025-08-20_050909_cached_res")) # 300
l$sl_360$short_i4_h4 <- bivar(paste0(run, "2025-08-20_052245_cached_res")) # 360



## Shorter fund lifetime max(investing=5 + max.holding=5) ----
run <- "cache_q_factors_20250808_223358_simulated_cashflows_EW_VYP/"

l$sl_1$short_i5_h5 <-   bivar(paste0(run, "2025-08-18_195614_cached_res")) # 1
l$sl_30$short_i5_h5 <-  bivar(paste0(run, "2025-08-21_125236_cached_res")) # 30
l$sl_60$short_i5_h5 <-  bivar(paste0(run, "2025-08-18_200649_cached_res")) # 60
l$sl_90$short_i5_h5 <-  bivar(paste0(run, "2025-08-21_130325_cached_res")) # 90
l$sl_120$short_i5_h5 <- bivar(paste0(run, "2025-08-18_201811_cached_res")) # 120
l$sl_150$short_i5_h5 <- bivar(paste0(run, "2025-08-18_203005_cached_res")) # 150
l$sl_180$short_i5_h5 <- bivar(paste0(run, "2025-08-18_204205_cached_res")) # 180
l$sl_210$short_i5_h5 <- bivar(paste0(run, "2025-08-18_205414_cached_res")) # 210
l$sl_240$short_i5_h5 <- bivar(paste0(run, "2025-08-18_210639_cached_res")) # 240
l$sl_300$short_i5_h5 <- bivar(paste0(run, "2025-08-18_211921_cached_res")) # 300
l$sl_360$short_i5_h5 <- bivar(paste0(run, "2025-08-18_213309_cached_res")) # 360

## Shorter fund lifetime max(investing=6 + max.holding=6) ----

run <- "cache_q_factors_20250821_154353_simulated_cashflows_EW_VYP/"
l$sl_1$short_i6_h6 <-   bivar(paste0(run, "2025-08-21_230836_cached_res")) # 1
l$sl_30$short_i6_h6 <-  bivar(paste0(run, "2025-08-21_231918_cached_res")) # 30
l$sl_60$short_i6_h6 <-  bivar(paste0(run, "2025-08-21_233021_cached_res")) # 60
l$sl_90$short_i6_h6 <-  bivar(paste0(run, "2025-08-21_234141_cached_res")) # 90
l$sl_120$short_i6_h6 <- bivar(paste0(run, "2025-08-21_235310_cached_res")) # 120
l$sl_150$short_i6_h6 <- bivar(paste0(run, "2025-08-22_000505_cached_res")) # 150
l$sl_180$short_i6_h6 <- bivar(paste0(run, "2025-08-22_001721_cached_res")) # 180
l$sl_210$short_i6_h6 <- bivar(paste0(run, "2025-08-22_002954_cached_res")) # 210
l$sl_240$short_i6_h6 <- bivar(paste0(run, "2025-08-22_004246_cached_res")) # 240
l$sl_300$short_i6_h6 <- bivar(paste0(run, "2025-08-22_005548_cached_res")) # 300
l$sl_360$short_i6_h6 <- bivar(paste0(run, "2025-08-22_010921_cached_res")) # 360

## Shorter fund lifetime max(investing=7 + max.holding=7) ----

run <- "cache_q_factors_20250821_155431_simulated_cashflows_EW_VYP/"
l$sl_1$short_i7_h7 <-   bivar(paste0(run, "2025-08-22_012154_cached_res")) # 1
l$sl_30$short_i7_h7 <-  bivar(paste0(run, "2025-08-22_013243_cached_res")) # 30
l$sl_60$short_i7_h7 <-  bivar(paste0(run, "2025-08-22_014334_cached_res")) # 60
l$sl_90$short_i7_h7 <-  bivar(paste0(run, "2025-08-22_015442_cached_res")) # 90
l$sl_120$short_i7_h7 <- bivar(paste0(run, "2025-08-22_020603_cached_res")) # 120
l$sl_150$short_i7_h7 <- bivar(paste0(run, "2025-08-22_021758_cached_res")) # 150
l$sl_180$short_i7_h7 <- bivar(paste0(run, "2025-08-22_023018_cached_res")) # 180
l$sl_210$short_i7_h7 <- bivar(paste0(run, "2025-08-22_024244_cached_res")) # 210
l$sl_240$short_i7_h7 <- bivar(paste0(run, "2025-08-22_025530_cached_res")) # 240
l$sl_300$short_i7_h7 <- bivar(paste0(run, "2025-08-22_030839_cached_res")) # 300
l$sl_360$short_i7_h7 <- bivar(paste0(run, "2025-08-22_032211_cached_res")) # 360

## Shorter fund lifetime max(investing=8 + max.holding=8) ----

run <- "cache_q_factors_20250821_160511_simulated_cashflows_EW_VYP/"
l$sl_1$short_i8_h8 <-   bivar(paste0(run, "2025-08-22_035803_cached_res")) # 1
l$sl_30$short_i8_h8 <-  bivar(paste0(run, "2025-08-22_040832_cached_res")) # 30
l$sl_60$short_i8_h8 <-  bivar(paste0(run, "2025-08-22_041909_cached_res")) # 60
l$sl_90$short_i8_h8 <-  bivar(paste0(run, "2025-08-22_042954_cached_res")) # 90
l$sl_120$short_i8_h8 <- bivar(paste0(run, "2025-08-22_044054_cached_res")) # 120
l$sl_150$short_i8_h8 <- bivar(paste0(run, "2025-08-22_045219_cached_res")) # 150
l$sl_180$short_i8_h8 <- bivar(paste0(run, "2025-08-22_050407_cached_res")) # 180
l$sl_210$short_i8_h8 <- bivar(paste0(run, "2025-08-22_051609_cached_res")) # 210
l$sl_240$short_i8_h8 <- bivar(paste0(run, "2025-08-22_052827_cached_res")) # 240
l$sl_300$short_i8_h8 <- bivar(paste0(run, "2025-08-22_054106_cached_res")) # 300
l$sl_360$short_i8_h8 <- bivar(paste0(run, "2025-08-22_055417_cached_res")) # 360

## Shorter fund lifetime max(investing=9 + max.holding=9) ----

run <- "cache_q_factors_20250821_161553_simulated_cashflows_EW_VYP/"
l$sl_1$short_i9_h9 <-   bivar(paste0(run, "2025-08-22_060641_cached_res")) # 1
l$sl_30$short_i9_h9 <-  bivar(paste0(run, "2025-08-22_061705_cached_res")) # 30
l$sl_60$short_i9_h9 <-  bivar(paste0(run, "2025-08-22_062734_cached_res")) # 60
l$sl_90$short_i9_h9 <-  bivar(paste0(run, "2025-08-22_063816_cached_res")) # 90
l$sl_120$short_i9_h9 <- bivar(paste0(run, "2025-08-22_064912_cached_res")) # 120
l$sl_150$short_i9_h9 <- bivar(paste0(run, "2025-08-22_070029_cached_res")) # 150
l$sl_180$short_i9_h9 <- bivar(paste0(run, "2025-08-22_071208_cached_res")) # 180
l$sl_210$short_i9_h9 <- bivar(paste0(run, "2025-08-22_072410_cached_res")) # 210
l$sl_240$short_i9_h9 <- bivar(paste0(run, "2025-08-22_073630_cached_res")) # 240
l$sl_300$short_i9_h9 <- bivar(paste0(run, "2025-08-22_074918_cached_res")) # 300
l$sl_360$short_i9_h9 <- bivar(paste0(run, "2025-08-22_080240_cached_res")) # 360

## Normal fund lifetime max(investing=5 + max.holding=10) ------
l$sl_1$short_i5_h10 <-l$sl_1$one
l$sl_30$short_i5_h10 < NA
l$sl_60$short_i5_h10 <-l$sl_60$one 
l$sl_90$short_i5_h10 <-NA
l$sl_120$short_i5_h10 <-l$sl_120$one 
l$sl_150$short_i5_h10 <-l$sl_150$one 
l$sl_180$short_i5_h10 <-l$sl_180$one 
l$sl_210$short_i5_h10 <-l$sl_210$one 
l$sl_240$short_i5_h10 <-l$sl_240$one 
l$sl_300$short_i5_h10 <-l$sl_300$one
l$sl_360$short_i5_h10 <-l$sl_360$one

## Plot shorties (x+y=6) -----

# Build + plot "shorter fund lifetime" comparison in one chart
plot_shorter_lifetime <- function(l,
                                  ylim = c(0, 2),
                                  do.eps = FALSE,
                                  eps_file = "chart/Simulation_shorter_fund_lifetime.eps",
                                  xlab.text = "Max Months") {
  # Expect the 3+3 series to be named short_i3_h3 (fix your inputs if duplicated)
  if (is.null(l$sl_1$short_i3_h3)) {
    stop("Rename your 3+3 block to 'short_i3_h3' (your paste repeats 'short_i5_h1').")
  }
  
  months <- c(1, 30, 60, 90, 120, 150, 180, 210, 240, 300, 360)
  
  # Keys in the list and their human-readable labels (legend)
  series_keys <- c("short_i1_h5", "short_i2_h4", "short_i3_h3", "short_i4_h2", "short_i5_h1")
  series_labels <- c("Invest 1 / Hold 5",
                     "Invest 2 / Hold 4",
                     "Invest 3 / Hold 3",
                     "Invest 4 / Hold 2",
                     "Invest 5 / Hold 1")
  stopifnot(length(series_keys) == length(series_labels))
  
  # Colors for each scenario (feel free to tweak)
  series_cols <- c("black", "blue", "darkgreen", "red", "purple")
  
  # Helper to extract a vector for Mean/SD across months
  get_vec <- function(key, stat) {
    sapply(months, function(x) l[[paste0("sl_", x)]][[key]]$MKT$MKT[[stat]])
  }
  
  # Build a tidy-ish data.frame you can also inspect if you like
  df <- data.frame(months = months)
  for (k in series_keys) {
    df[[paste0(k, "_Mean")]] <- get_vec(k, "Mean")
    df[[paste0(k, "_SD")]]   <- get_vec(k, "SD")
  }
  print(df)  # optional: see the numbers
  
  # Optional EPS export
  if (do.eps) {
    dir.create(dirname(eps_file), showWarnings = FALSE, recursive = TRUE)
    setEPS()
    postscript(eps_file, width = 6, height = 3.5, family = "Helvetica", pointsize = 10)
    on.exit(dev.off(), add = TRUE)
  }
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(4.2, 4.6, 1.2, 1), xpd = FALSE)
  
  pch.mean <- 19
  pch.stdv <- 1

  # Start with the first series (Mean)
  y0 <- df[[paste0(series_keys[1], "_Mean")]]
  plot(months, y0, type = "b", pch=pch.mean, col = series_cols[1],
       ylim = ylim, xlab = xlab.text, ylab = "Mean and Stdv of Estimates")
  
  # Add the rest (Means as solid; SDs as dashed, same colors)
  # First, overlay SD for the first series so it matches loop behavior
  lines(months, df[[paste0(series_keys[1], "_SD")]], type = "b", pch=pch.stdv, lty = 3, col = series_cols[1])
  
  if (length(series_keys) > 1) {
    for (i in 2:length(series_keys)) {
      k <- series_keys[i]; col_i <- series_cols[i]
      lines(months, df[[paste0(k, "_Mean")]], type = "b", pch=pch.mean, col = col_i)
      lines(months, df[[paste0(k, "_SD")]],   type = "b", pch=pch.stdv, lty = 3, col = col_i)
    }
  }
  
  abline(h = c(0, 1), col = "grey")
  
  # Clean, unambiguous legends: one for Means, one for Std devs
  legend("topleft", bty = "n", title = "Mean (solid)",
         legend = series_labels, col = series_cols, lty = 1, pch=pch.mean, cex = 0.6)
  legend("topright", bty = "n", inset = c(0, 0), title = "Std dev (dashed)",
         legend = series_labels, col = series_cols, lty = 3, pch=pch.stdv, cex = 0.6)
}

# --- Usage ---
# After you load your l$sl_*$short_i*_h* objects (with 3+3 named short_i3_h3):
plot_shorter_lifetime(l, ylim = c(0, 2))
# or to save EPS:
# plot_shorter_lifetime(l, ylim = c(0, 2), do.eps = TRUE,
#                       eps_file = "chart/Simulation_shorter_fund_lifetime.eps")


## Plot shorties equal invest & holding periods ----

# Plot "equal investing/holding" shorter-lifetime scenarios in one chart
plot_short_equal_lifetime <- function(
    l,
    # keys present in your 'l$sl_*' lists:
    invest.holding.periods = c(2, 4, 6, 8),
    # candidate months; function will keep only those that actually exist in 'l'
    months_all = c(1, 30, 60, 90, 120, 150, 180, 210, 240, 300, 360),
    ylim = NULL,                # auto if NULL, otherwise e.g. c(0, 2)
    do.eps = FALSE,
    eps_file = "chart/Simulation_short_equal_lifetime.eps",
    xlab.text = "Max Months"
) {
  
  series <- c()
  labels <- c()
  for (h in invest.holding.periods) {
    series <- c(series, paste0("short_i", h, "_h", h))
    labels <- c(labels, paste("Invest", h, "/", "Hold", h))
  }
  
  stopifnot(length(series) == length(labels))
  
  # Keep only months that exist for at least one of the requested series
  has_series_at <- function(m) {
    node <- l[[paste0("sl_", m)]]
    if (is.null(node)) return(FALSE)
    any(sapply(series, function(k) !is.null(node[[k]])))
  }
  months <- months_all[sapply(months_all, has_series_at)]
  if (length(months) == 0) stop("No months found in 'l' for the given series.")
  
  # Safe extractor (returns NA if a path is missing)
  get_val <- function(m, key, stat) {
    node <- l[[paste0("sl_", m)]]
    if (is.null(node) || is.null(node[[key]])) return(NA_real_)
    x <- try(node[[key]]$MKT$MKT[[stat]], silent = TRUE)
    if (inherits(x, "try-error") || is.null(x)) return(NA_real_)
    x
  }
  
  # Build a compact data frame (also handy for inspection)
  df <- data.frame(months = months)
  for (k in series) {
    df[[paste0(k, "_Mean")]] <- sapply(months, get_val, key = k, stat = "Mean")
    df[[paste0(k, "_SD")]]   <- sapply(months, get_val, key = k, stat = "SD")
  }
  print(df)  # optional
  
  # Auto y-limits if not provided
  if (is.null(ylim)) {
    rng <- range(as.matrix(df[,-1]), na.rm = TRUE)
    pad <- diff(rng) * 0.06
    if (!is.finite(pad)) pad <- 0.1
    ylim <- c(max(0, rng[1] - pad), rng[2] + pad)
  }
  
  # Optional EPS export
  if (do.eps) {
    dir.create(dirname(eps_file), showWarnings = FALSE, recursive = TRUE)
    setEPS()
    postscript(eps_file, width = 6, height = 3.5, family = "Helvetica", pointsize = 10)
    on.exit(dev.off(), add = TRUE)
  }
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(4.2, 4.6, 1.2, 1), xpd = FALSE)
  
  cols <- c("black", "blue", "darkgreen", "red", "orange")[seq_along(series)]
  pch.mean <- 19
  pch.stdv <- 1
  
  # Start plot with the series that has the most available points
  idx_start <- which.max(sapply(series, function(k) sum(!is.na(df[[paste0(k, "_Mean")]]))))
  y0 <- df[[paste0(series[idx_start], "_Mean")]]
  plot(df$months, y0, type = "b", col = cols[idx_start], pch = pch.mean,
       ylim = ylim, xlab = xlab.text, ylab = "Mean and Stdv of Estimates")
  lines(df$months, df[[paste0(series[idx_start], "_SD")]], type = "b", lty = 3, col = cols[idx_start], pch = pch.stdv)
  
  # Add remaining series
  for (i in setdiff(seq_along(series), idx_start)) {
    lines(df$months, df[[paste0(series[i], "_Mean")]], type = "b", col = cols[i], pch = pch.mean)
    lines(df$months, df[[paste0(series[i], "_SD")]],   type = "b", lty = 3, col = cols[i], pch = pch.stdv)
  }
  
  abline(h = c(0, 1), col = "grey")
  
  # Clear, unambiguous legends
  legend("topleft", bty = "n", title = "Mean (solid)",
         legend = labels, col = cols, lty = 1, pch = pch.mean, cex = 0.6)
  legend("topright", bty = "n", inset = c(0, 0.0), title = "Std dev (dashed)",
         legend = labels, col = cols, lty = 3, pch = pch.stdv, cex = 0.6)
}

# --- Usage ---
# After loading your 'short', 'short_i2_h2', 'short_i3_h3', 'short_i4_h4' series:
# plot_short_equal_lifetime(l)                  # auto y-lims
plot_short_equal_lifetime(l, ylim = c(0, 2))  # fixed y-lims like your earlier chart
# plot_short_equal_lifetime(l, do.eps = TRUE,
#   eps_file = "chart/Simulation_short_equal_lifetime.eps")


## Plot investment period = 5 & holding periods 1,5,10 ----

# Plot "equal investing/holding" shorter-lifetime scenarios in one chart
plot_investment_peri_5 <- function(
    l,
    # keys present in your 'l$sl_*' lists:
    invest.holding.periods = c(1, 5, 10),
    # candidate months; function will keep only those that actually exist in 'l'
    months_all = c(1, 30, 60, 90, 120, 150, 180, 210, 240, 300, 360),
    ylim = NULL,                # auto if NULL, otherwise e.g. c(0, 2)
    do.eps = FALSE,
    eps_file = "chart/Simulation_short_equal_lifetime.eps",
    xlab.text = "Max Months"
) {
  
  series <- c()
  labels <- c()
  for (h in invest.holding.periods) {
    series <- c(series, paste0("short_i", 5, "_h", h))
    labels <- c(labels, paste("Invest", 5, "/", "Hold", h))
  }
  
  stopifnot(length(series) == length(labels))
  
  # Keep only months that exist for at least one of the requested series
  has_series_at <- function(m) {
    node <- l[[paste0("sl_", m)]]
    if (is.null(node)) return(FALSE)
    any(sapply(series, function(k) !is.null(node[[k]])))
  }
  months <- months_all[sapply(months_all, has_series_at)]
  if (length(months) == 0) stop("No months found in 'l' for the given series.")
  
  # Safe extractor (returns NA if a path is missing)
  get_val <- function(m, key, stat) {
    node <- l[[paste0("sl_", m)]]
    if (is.null(node) || is.null(node[[key]])) return(NA_real_)
    x <- try(node[[key]]$MKT$MKT[[stat]], silent = TRUE)
    if (inherits(x, "try-error") || is.null(x)) return(NA_real_)
    x
  }
  
  # Build a compact data frame (also handy for inspection)
  df <- data.frame(months = months)
  for (k in series) {
    df[[paste0(k, "_Mean")]] <- sapply(months, get_val, key = k, stat = "Mean")
    df[[paste0(k, "_SD")]]   <- sapply(months, get_val, key = k, stat = "SD")
  }
  print(df)  # optional
  
  # Auto y-limits if not provided
  if (is.null(ylim)) {
    rng <- range(as.matrix(df[,-1]), na.rm = TRUE)
    pad <- diff(rng) * 0.06
    if (!is.finite(pad)) pad <- 0.1
    ylim <- c(max(0, rng[1] - pad), rng[2] + pad)
  }
  
  # Optional EPS export
  if (do.eps) {
    dir.create(dirname(eps_file), showWarnings = FALSE, recursive = TRUE)
    setEPS()
    postscript(eps_file, width = 6, height = 3.5, family = "Helvetica", pointsize = 10)
    on.exit(dev.off(), add = TRUE)
  }
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(4.2, 4.6, 1.2, 1), xpd = FALSE)
  
  cols <- c("black", "blue", "darkgreen", "red", "orange")[seq_along(series)]
  pch.mean <- 19
  pch.stdv <- 1
  
  make.xy <- function(dfs, series, i, stats) {
    y.col <- paste0(series[i], stats)
    dfs <- dfs[, c("months", y.col)]
    dfs <- dfs[complete.cases(dfs), ]
    return(list(x = dfs$months, y = dfs[, y.col]))
  }
  
  # Start plot with the series that has the most available points
  idx_start <- which.max(sapply(series, function(k) sum(!is.na(df[[paste0(k, "_Mean")]]))))
  vals <- make.xy(df, series, idx_start, "_Mean")
  plot(vals$x, vals$y, type = "b", col = cols[idx_start], pch = pch.mean,
       ylim = ylim, xlab = xlab.text, ylab = "Mean and Stdv of Estimates")
  vals <- make.xy(df, series, idx_start, "_SD")
  lines(vals$x, vals$y, type = "b", lty = 3, col = cols[idx_start], pch = pch.stdv)
  
  # Add remaining series
  for (i in setdiff(seq_along(series), idx_start)) {
    vals <- make.xy(df, series, i, "_Mean")
    lines(vals$x, vals$y, type = "b", col = cols[i], pch = pch.mean)
    vals <- make.xy(df, series, i, "_SD")
    lines(vals$x, vals$y, type = "b", lty = 3, col = cols[i], pch = pch.stdv)
  }
  
  abline(h = c(0, 1), col = "grey")
  
  # Clear, unambiguous legends
  legend("topleft", bty = "n", title = "Mean (solid)",
         legend = labels, col = cols, lty = 1, pch = pch.mean, cex = 0.6)
  legend("topright", bty = "n", inset = c(0, 0.0), title = "Std dev (dashed)",
         legend = labels, col = cols, lty = 3, pch = pch.stdv, cex = 0.6)
}

# --- Usage ---
# After loading your 'short', 'short_i2_h2', 'short_i3_h3', 'short_i4_h4' series:
# plot_investment_peri_5(l)                  # auto y-lims
plot_investment_peri_5(l, ylim = c(0, 2))  # fixed y-lims like your earlier chart
# plot_investment_peri_5(l, do.eps = TRUE,
#   eps_file = "chart/Simulation_plot_investment_peri_5.eps")



## No investment period (investing=0 + max.holding=1) ----

run <- "cache_q_factors_20250825_143009_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h1 <-   bivar(paste0(run, "2025-08-25_180049_cached_res")) # 1
l$sl_30$no_i0_h1 <-  bivar(paste0(run, "2025-08-25_180534_cached_res")) # 30
l$sl_60$no_i0_h1 <-  bivar(paste0(run, "2025-08-25_181023_cached_res")) # 60
l$sl_90$no_i0_h1 <-  bivar(paste0(run, "2025-08-25_181526_cached_res")) # 90
l$sl_120$no_i0_h1 <- bivar(paste0(run, "2025-08-25_182049_cached_res")) # 120
l$sl_150$no_i0_h1 <- bivar(paste0(run, "2025-08-25_182704_cached_res")) # 150
l$sl_180$no_i0_h1 <- bivar(paste0(run, "2025-08-25_183341_cached_res")) # 180
l$sl_210$no_i0_h1 <- bivar(paste0(run, "2025-08-25_184031_cached_res")) # 210
l$sl_240$no_i0_h1 <- bivar(paste0(run, "2025-08-25_184737_cached_res")) # 240
l$sl_300$no_i0_h1 <- bivar(paste0(run, "2025-08-25_185510_cached_res")) # 300
l$sl_360$no_i0_h1 <- bivar(paste0(run, "2025-08-25_190300_cached_res")) # 360

## No investment period (investing=0 + max.holding=3) ----

run <- "cache_q_factors_20250825_143734_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h3 <-   bivar(paste0(run, "2025-08-25_190925_cached_res")) # 1
l$sl_30$no_i0_h3 <-  bivar(paste0(run, "2025-08-25_191418_cached_res")) # 30
l$sl_60$no_i0_h3 <-  bivar(paste0(run, "2025-08-25_191918_cached_res")) # 60
l$sl_90$no_i0_h3 <-  bivar(paste0(run, "2025-08-25_192427_cached_res")) # 90
l$sl_120$no_i0_h3 <- bivar(paste0(run, "2025-08-25_192952_cached_res")) # 120
l$sl_150$no_i0_h3 <- bivar(paste0(run, "2025-08-25_193542_cached_res")) # 150
l$sl_180$no_i0_h3 <- bivar(paste0(run, "2025-08-25_194150_cached_res")) # 180
l$sl_210$no_i0_h3 <- bivar(paste0(run, "2025-08-25_194804_cached_res")) # 210
l$sl_240$no_i0_h3 <- bivar(paste0(run, "2025-08-25_195421_cached_res")) # 240
l$sl_300$no_i0_h3 <- bivar(paste0(run, "2025-08-25_200045_cached_res")) # 300
l$sl_360$no_i0_h3 <- bivar(paste0(run, "2025-08-25_200738_cached_res")) # 360

## No investment period (investing=0 + max.holding=5) ----

run <- "cache_q_factors_20250825_144508_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h5 <-   bivar(paste0(run, "2025-08-25_201416_cached_res")) # 1
l$sl_30$no_i0_h5 <-  bivar(paste0(run, "2025-08-25_201923_cached_res")) # 30
l$sl_60$no_i0_h5 <-  bivar(paste0(run, "2025-08-25_202436_cached_res")) # 60
l$sl_90$no_i0_h5 <-  bivar(paste0(run, "2025-08-25_202954_cached_res")) # 90
l$sl_120$no_i0_h5 <- bivar(paste0(run, "2025-08-25_203524_cached_res")) # 120
l$sl_150$no_i0_h5 <- bivar(paste0(run, "2025-08-25_204111_cached_res")) # 150
l$sl_180$no_i0_h5 <- bivar(paste0(run, "2025-08-25_204713_cached_res")) # 180
l$sl_210$no_i0_h5 <- bivar(paste0(run, "2025-08-25_205321_cached_res")) # 210
l$sl_240$no_i0_h5 <- bivar(paste0(run, "2025-08-25_205934_cached_res")) # 240
l$sl_300$no_i0_h5 <- bivar(paste0(run, "2025-08-25_210550_cached_res")) # 300
l$sl_360$no_i0_h5 <- bivar(paste0(run, "2025-08-25_211230_cached_res")) # 360

## No investment period (investing=0 + max.holding=7) ----

run <- "cache_q_factors_20250825_145244_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h7 <-   bivar(paste0(run, "2025-08-25_211901_cached_res")) # 1
l$sl_30$no_i0_h7 <-  bivar(paste0(run, "2025-08-25_212403_cached_res")) # 30
l$sl_60$no_i0_h7 <-  bivar(paste0(run, "2025-08-25_212913_cached_res")) # 60
l$sl_90$no_i0_h7 <-  bivar(paste0(run, "2025-08-25_213427_cached_res")) # 90
l$sl_120$no_i0_h7 <- bivar(paste0(run, "2025-08-25_213951_cached_res")) # 120
l$sl_150$no_i0_h7 <- bivar(paste0(run, "2025-08-25_214531_cached_res")) # 150
l$sl_180$no_i0_h7 <- bivar(paste0(run, "2025-08-25_215119_cached_res")) # 180
l$sl_210$no_i0_h7 <- bivar(paste0(run, "2025-08-25_215714_cached_res")) # 210
l$sl_240$no_i0_h7 <- bivar(paste0(run, "2025-08-25_220313_cached_res")) # 240
l$sl_300$no_i0_h7 <- bivar(paste0(run, "2025-08-25_220922_cached_res")) # 300
l$sl_360$no_i0_h7 <- bivar(paste0(run, "2025-08-25_221552_cached_res")) # 360

## No investment period (investing=0 + max.holding=10) ----

run <- "cache_q_factors_20250825_150023_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h10 <-   bivar(paste0(run, "2025-08-25_222241_cached_res")) # 1
l$sl_30$no_i0_h10 <-  bivar(paste0(run, "2025-08-25_222755_cached_res")) # 30
l$sl_60$no_i0_h10 <-  bivar(paste0(run, "2025-08-25_223315_cached_res")) # 60
l$sl_90$no_i0_h10 <-  bivar(paste0(run, "2025-08-25_223840_cached_res")) # 90
l$sl_120$no_i0_h10 <- bivar(paste0(run, "2025-08-25_224413_cached_res")) # 120
l$sl_150$no_i0_h10 <- bivar(paste0(run, "2025-08-25_224955_cached_res")) # 150
l$sl_180$no_i0_h10 <- bivar(paste0(run, "2025-08-25_225546_cached_res")) # 180
l$sl_210$no_i0_h10 <- bivar(paste0(run, "2025-08-25_230141_cached_res")) # 210
l$sl_240$no_i0_h10 <- bivar(paste0(run, "2025-08-25_230740_cached_res")) # 240
l$sl_300$no_i0_h10 <- bivar(paste0(run, "2025-08-25_231352_cached_res")) # 300
l$sl_360$no_i0_h10 <- bivar(paste0(run, "2025-08-25_232022_cached_res")) # 360

# 13) No investment period (investing=0) - alpha & beta -----

### 13.a) No investment period  (investing=0 + max.holding=1) 

run <- "cache_q_factors_20250825_143009_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h1_alpha <-   bivar(paste0(run, "2025-08-26_143511_cached_res")) # 1
l$sl_30$no_i0_h1_alpha <-  bivar(paste0(run, "2025-08-26_144425_cached_res")) # 30
l$sl_60$no_i0_h1_alpha <-  bivar(paste0(run, "2025-08-26_145336_cached_res")) # 60
l$sl_90$no_i0_h1_alpha <-  bivar(paste0(run, "2025-08-26_150232_cached_res")) # 90
l$sl_120$no_i0_h1_alpha <- bivar(paste0(run, "2025-08-26_151214_cached_res")) # 120
l$sl_150$no_i0_h1_alpha <- bivar(paste0(run, "2025-08-26_152245_cached_res")) # 150
l$sl_180$no_i0_h1_alpha <- bivar(paste0(run, "2025-08-26_153324_cached_res")) # 180
l$sl_210$no_i0_h1_alpha <- bivar(paste0(run, "2025-08-26_154340_cached_res")) # 210
l$sl_240$no_i0_h1_alpha <- bivar(paste0(run, "2025-08-26_155358_cached_res")) # 240
l$sl_300$no_i0_h1_alpha <- bivar(paste0(run, "2025-08-26_160437_cached_res")) # 300
l$sl_360$no_i0_h1_alpha <- bivar(paste0(run, "2025-08-26_161521_cached_res")) # 360

### 13.b) No investment period  (investing=0 + max.holding=3) 

run <- "cache_q_factors_20250825_143734_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h3_alpha <-   bivar(paste0(run, "2025-08-26_162732_cached_res")) # 1
l$sl_30$no_i0_h3_alpha <-  bivar(paste0(run, "2025-08-26_163814_cached_res")) # 30
l$sl_60$no_i0_h3_alpha <-  bivar(paste0(run, "2025-08-26_164858_cached_res")) # 60
l$sl_90$no_i0_h3_alpha <-  bivar(paste0(run, "2025-08-26_165944_cached_res")) # 90
l$sl_120$no_i0_h3_alpha <- bivar(paste0(run, "2025-08-26_171033_cached_res")) # 120
l$sl_150$no_i0_h3_alpha <- bivar(paste0(run, "2025-08-26_172154_cached_res")) # 150
l$sl_180$no_i0_h3_alpha <- bivar(paste0(run, "2025-08-26_173317_cached_res")) # 180
l$sl_210$no_i0_h3_alpha <- bivar(paste0(run, "2025-08-26_174259_cached_res")) # 210
l$sl_240$no_i0_h3_alpha <- bivar(paste0(run, "2025-08-26_175207_cached_res")) # 240
l$sl_300$no_i0_h3_alpha <- bivar(paste0(run, "2025-08-26_180107_cached_res")) # 300
l$sl_360$no_i0_h3_alpha <- bivar(paste0(run, "2025-08-26_181015_cached_res")) # 360

### 13.c) No investment period  (investing=0 + max.holding=5) 

run <- "cache_q_factors_20250825_144508_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h5_alpha <-   bivar(paste0(run, "2025-08-26_182212_cached_res")) # 1
l$sl_30$no_i0_h5_alpha <-  bivar(paste0(run, "2025-08-26_183255_cached_res")) # 30
l$sl_60$no_i0_h5_alpha <-  bivar(paste0(run, "2025-08-26_184342_cached_res")) # 60
l$sl_90$no_i0_h5_alpha <-  bivar(paste0(run, "2025-08-26_185426_cached_res")) # 90
l$sl_120$no_i0_h5_alpha <- bivar(paste0(run, "2025-08-26_190517_cached_res")) # 120
l$sl_150$no_i0_h5_alpha <- bivar(paste0(run, "2025-08-26_191640_cached_res")) # 150
l$sl_180$no_i0_h5_alpha <- bivar(paste0(run, "2025-08-26_192811_cached_res")) # 180
l$sl_210$no_i0_h5_alpha <- bivar(paste0(run, "2025-08-26_193924_cached_res")) # 210
l$sl_240$no_i0_h5_alpha <- bivar(paste0(run, "2025-08-26_195033_cached_res")) # 240
l$sl_300$no_i0_h5_alpha <- bivar(paste0(run, "2025-08-26_200106_cached_res")) # 300
l$sl_360$no_i0_h5_alpha <- bivar(paste0(run, "2025-08-26_201104_cached_res")) # 360

### 13.d) No investment period  (investing=0 + max.holding=7) 

run <- "cache_q_factors_20250825_145244_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h7_alpha <-   bivar(paste0(run, "2025-08-26_202332_cached_res")) # 1
l$sl_30$no_i0_h7_alpha <-  bivar(paste0(run, "2025-08-26_203429_cached_res")) # 30
l$sl_60$no_i0_h7_alpha <-  bivar(paste0(run, "2025-08-26_204536_cached_res")) # 60
l$sl_90$no_i0_h7_alpha <-  bivar(paste0(run, "2025-08-26_205637_cached_res")) # 90
l$sl_120$no_i0_h7_alpha <- bivar(paste0(run, "2025-08-26_210744_cached_res")) # 120
l$sl_150$no_i0_h7_alpha <- bivar(paste0(run, "2025-08-26_211925_cached_res")) # 150
l$sl_180$no_i0_h7_alpha <- bivar(paste0(run, "2025-08-26_213043_cached_res")) # 180
l$sl_210$no_i0_h7_alpha <- bivar(paste0(run, "2025-08-26_214146_cached_res")) # 210
l$sl_240$no_i0_h7_alpha <- bivar(paste0(run, "2025-08-26_215234_cached_res")) # 240
l$sl_300$no_i0_h7_alpha <- bivar(paste0(run, "2025-08-26_220301_cached_res")) # 300
l$sl_360$no_i0_h7_alpha <- bivar(paste0(run, "2025-08-26_221354_cached_res")) # 360

### 13.e) No investment period  (investing=0 + max.holding=10) 

run <- "cache_q_factors_20250825_150023_simulated_cashflows_EW_VYP/"
l$sl_1$no_i0_h10_alpha <-   bivar(paste0(run, "2025-08-26_222612_cached_res")) # 1
l$sl_30$no_i0_h10_alpha <-  bivar(paste0(run, "2025-08-26_223747_cached_res")) # 30
l$sl_60$no_i0_h10_alpha <-  bivar(paste0(run, "2025-08-26_224936_cached_res")) # 60
l$sl_90$no_i0_h10_alpha <-  bivar(paste0(run, "2025-08-26_230142_cached_res")) # 90
l$sl_120$no_i0_h10_alpha <- bivar(paste0(run, "2025-08-26_231349_cached_res")) # 120
l$sl_150$no_i0_h10_alpha <- bivar(paste0(run, "2025-08-26_232550_cached_res")) # 150
l$sl_180$no_i0_h10_alpha <- bivar(paste0(run, "2025-08-26_233735_cached_res")) # 180
l$sl_210$no_i0_h10_alpha <- bivar(paste0(run, "2025-08-26_234850_cached_res")) # 210
l$sl_240$no_i0_h10_alpha <- bivar(paste0(run, "2025-08-26_235954_cached_res")) # 240
l$sl_300$no_i0_h10_alpha <- bivar(paste0(run, "2025-08-27_001018_cached_res")) # 300
l$sl_360$no_i0_h10_alpha <- bivar(paste0(run, "2025-08-27_002036_cached_res")) # 360

## Plot shorties 1:9 -----

# Old-style plot for equal investing/holding short_iX_hX (X = 1..9)
years1.9 <- 1:9
years1.9 <- c(1, 5, 9)
plot_short_equal_lifetime_oldstyle_1to9 <- function(
    l,
    series     = sprintf("short_i%d_h%d", years1.9, years1.9),
    labels     = sprintf("Invest %d / Hold %d", years1.9, years1.9),
    months_all = c(1, 30, 60, 90, 120, 150, 180, 210, 240, 300, 360),
    ylim       = NULL,                 # auto if NULL, else e.g. c(0, 2)
    do.eps     = FALSE,
    eps_file   = "chart/Simulation_short_equal_lifetime_i_eq_h_1to9.eps",
    xlab.text  = "Max Months",
    legend_cex = 0.60,                 # smaller legend to fit many entries
    legend_ncol= 2,                    # multi-column legend
    point_cex  = 0.6,                  # smaller points for readability
    pch_point  = 1                     # open circles (like base default)
) {
  stopifnot(length(series) == length(labels))
  
  # Keep only months that exist for at least one of the requested series
  has_series_at <- function(m) {
    node <- l[[paste0("sl_", m)]]
    if (is.null(node)) return(FALSE)
    any(sapply(series, function(k) !is.null(node[[k]])))
  }
  months <- months_all[sapply(months_all, has_series_at)]
  if (length(months) == 0) stop("No months found in 'l' for the given series.")
  
  # Safe extractor
  get_val <- function(m, key, stat) {
    node <- l[[paste0("sl_", m)]]
    if (is.null(node) || is.null(node[[key]])) return(NA_real_)
    x <- try(node[[key]]$MKT$MKT[[stat]], silent = TRUE)
    if (inherits(x, "try-error") || is.null(x)) return(NA_real_)
    x
  }
  
  # Build data frame
  df <- data.frame(months = months)
  for (k in series) {
    df[[paste0(k, "_Mean")]] <- sapply(months, get_val, key = k, stat = "Mean")
    df[[paste0(k, "_SD")]]   <- sapply(months, get_val, key = k, stat = "SD")
  }
  
  # Drop any series that are completely NA
  keep <- sapply(series, function(k) {
    any(!is.na(df[[paste0(k, "_Mean")]])) || any(!is.na(df[[paste0(k, "_SD")]]))
  })
  series <- series[keep]
  labels <- labels[keep]
  if (length(series) == 0) stop("All requested series are empty.")
  
  # Auto y-limits (like your old chart, with a small pad)
  if (is.null(ylim)) {
    cols_mat <- unlist(lapply(series, function(k) c(paste0(k,"_Mean"), paste0(k,"_SD"))))
    rng <- range(as.matrix(df[, cols_mat, drop = FALSE]), na.rm = TRUE)
    pad <- diff(rng) * 0.06
    if (!is.finite(pad)) pad <- 0.1
    ylim <- c(max(0, rng[1] - pad), rng[2] + pad)
  }
  
  # Optional EPS export
  if (do.eps) {
    dir.create(dirname(eps_file), showWarnings = FALSE, recursive = TRUE)
    setEPS()
    postscript(eps_file, width = 5.5, height = 3.0, family = "Helvetica", pointsize = 10)
    on.exit(dev.off(), add = TRUE)
  }
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(4.2, 4.2, 1, 1), xpd = NA)
  
  # Fixed, readable palette (mean & sd share same color)
  base_cols <- c("black","royalblue","firebrick","darkgreen",
                 "darkorange3","purple","brown","goldenrod3","deeppink3")
  cols <- base_cols[seq_along(series)]
  
  # Start with the series that has the most points
  npts <- sapply(series, function(k) sum(!is.na(df[[paste0(k, "_Mean")]])))
  i0 <- which.max(npts)
  k0 <- series[i0]
  
  # Plot first series (Mean solid, SD dashed), with points (type="b")
  plot(df$months, df[[paste0(k0, "_Mean")]], type = "b",
       col = cols[i0], pch = pch_point, cex = point_cex, lwd = 1,
       ylim = ylim, xlab = xlab.text, ylab = "Mean and Stdv of Estimates")
  lines(df$months, df[[paste0(k0, "_SD")]], type = "b",
        col = cols[i0], pch = pch_point, cex = point_cex, lwd = 1, lty = 3)
  
  # Add remaining series
  for (i in setdiff(seq_along(series), i0)) {
    k <- series[i]
    lines(df$months, df[[paste0(k, "_Mean")]], type = "b",
          col = cols[i], pch = pch_point, cex = point_cex, lwd = 1)
    lines(df$months, df[[paste0(k, "_SD")]], type = "b",
          col = cols[i], pch = pch_point, cex = point_cex, lwd = 1, lty = 3)
  }
  
  abline(h = c(0, 1), col = "grey")
  
  # Single, old-style combined legend: all Means first (solid), then SDs (dashed)
  legend("topright", bty = "n",
         legend = c(paste(labels, "- Mean"), paste(labels, "- Stdv")),
         col    = rep(cols, times = 2),
         lty    = c(rep(1, length(series)), rep(3, length(series))),
         cex    = legend_cex,
         ncol   = legend_ncol)
}

# --- Usage ---
# After loading l$sl_*$short_iX_hX for X=1..9:
# plot_short_equal_lifetime_oldstyle_1to9(l)                 # auto y-lims
plot_short_equal_lifetime_oldstyle_1to9(l, ylim = c(0, 2)) # like your earlier chart
# plot_short_equal_lifetime_oldstyle_1to9(l, do.eps = TRUE,
#   eps_file = "chart/Simulation_short_equal_lifetime_i_eq_h_1to9.eps")

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



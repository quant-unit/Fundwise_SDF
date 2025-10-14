# Analyze simulated data

if(sys.nframe() == 0L) rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


df.11 <- read.csv("data_prepared_sim/20250821_153327/20250821_153327_simulated_cashflows_EW_VYP.csv")
df.99 <- read.csv("data_prepared_sim/20250821_161553/20250821_161553_simulated_cashflows_EW_VYP.csv")

unique(df.11$Fund.ID)[1]
unique(df.11$Fund.ID)[2]

i <- 1
df.x <- df.99

dates <- as.Date(df.x$Date[df.x$Fund.ID == unique(df.x$Fund.ID)[i]])
as.numeric(max(dates) - min(dates)) / 365.25

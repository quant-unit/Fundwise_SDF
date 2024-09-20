######################
### load data part of: estim_model_optimized.R
######################
# 0) Prologue -----
if(FALSE) {
  files <- c("prep_preqin.R", "prep_public.R")
  for(file in files) {
    source(file)
    if(sys.nframe() == 0L) rm(list = ls())
  }
}

if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 1.1) PARAMETERS ----
export.data <- FALSE

use.vintage.year.pfs <- TRUE
use.simulation <- FALSE
do.cross.validation <- FALSE
do.cache <- FALSE
do.parallel <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)

private.source <- "pitchbook"
# private.source <- "preqin"

#public.filename <- "public_returns" # outdated
public.filename <- "msci_market_factors"
#public.filename <- "q_factors"
#public.filename <- "DebtFactorsEURUSD" # outdated
#public.filename <- "iBoxxFactorsMIX"
#public.filename <- "iBoxxFactorsUSD" 
#public.filename <- "iBoxxFactorsEUR"


# CHOICES
weighting <- "EW"
weighting <- "FW"
error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

max.months <- c(120, 150, 180, 210, 240) # c(10, 20) * 12 # c(12.5, 15, 17.5) * 12

include.alpha.term <- FALSE
lambdas <- 0
kernel.bandwidth <- 12
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0("pitchbook_2023", ifelse(do.cross.validation, "_cv_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, ifelse(include.alpha.term, "alpha_", ""))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# 1.2) load data -----

# load public data
if (public.filename == "q_factors") {
  df.public <- read.csv(paste0("data_prepared/", public.filename, ".csv"))
  
} else {
  df.public <- read.csv2(paste0("data_prepared/", public.filename, ".csv"))
}
colnames(df.public) <- gsub("_World", "", colnames(df.public))
df.public$Date <- as.Date(df.public$Date)
bond.files <- c(
  "DebtFactorsUSD", "DebtFactorsEUR", "DebtFactorsEURUSD",
  "iBoxxFactorsEUR", "iBoxxFactorsUSD", "iBoxxFactorsMIX"
)
if (public.filename %in% bond.files ) {
  df.public <- df.public[complete.cases(df.public[, 2:5]), ]
  df.public[is.na(df.public)] <- 0
  public.equity <- "msci_market_factors"
  df.pubeq <- read.csv2(paste0("data_prepared/", public.equity, ".csv"))
  colnames(df.pubeq) <- gsub("_World", "", colnames(df.pubeq))
  df.pubeq$Date <- as.Date(df.pubeq$Date)
  
  df.public$RF <- NULL
  df.public <- merge(df.public, df.pubeq[, c("Date", "RF", "MKT")], by = "Date", all.x = TRUE)
  remove(df.pubeq)
  df.public <- df.public[df.public$Date > as.Date("2000-01-31"), ]
}

df.public$Alpha <- 1

if (export.data) {
  write.csv2(df.public, paste0("data_private_public/public_", public.filename, ".csv"))
}

# Load private data
if(!use.simulation) {
  
  if (private.source == "preqin") {
    year.tag <- ""
    df.private.cfs <- read.csv(paste0("data_prepared/preqin_cashflows_", weighting, year.tag, ".csv"))
    colnames(df.private.cfs)
  }
  
  if (private.source == "pitchbook") {
    year.tag <- "_2023"
    df.private.cfs <- read.csv(paste0("data_prepared/pitchbook_cashflows_", weighting, year.tag, ".csv"))
    colnames(df.private.cfs)
  }
  
} else {
  df.private.cfs <- read.csv(paste0("data_prepared/simulated_cashflows_", weighting, ".csv"))
}
df.private.cfs$Date <- as.Date(df.private.cfs$Date)
df.private.cfs$Fund.ID <- as.factor(paste(df.private.cfs$Fund.ID, df.private.cfs$type, sep = "_"))

to.monthly <- function(df.ss) {
  # fill zero cash flows (necessary for estimation)
  max.date <- max(df.ss$Date) + 1
  #max.date <- as.Date("2020-01-01")
  df.m <- data.frame(Date = seq(min(df.ss$Date) + 1, max.date, by = "month") - 1)
  df.ss <- merge(df.ss, df.m, by = "Date", all = TRUE)
  df.ss$CF[is.na(df.ss$CF)] <- 0
  for(col in c("type", "Vintage", "Fund.ID")) df.ss[, col] <- df.ss[1, col]
  return(df.ss)
}

df.private.cfs <- as.data.frame(data.table::rbindlist(lapply(split(df.private.cfs, df.private.cfs$Fund.ID), to.monthly)))
df.private.cfs$type <- as.factor(as.character(df.private.cfs$type))
df.private.cfs$Fund.ID <- as.factor(as.character(df.private.cfs$Fund.ID))
length(levels(df.private.cfs$type))



# merge private and public data
df0 <- merge(df.private.cfs, df.public, by="Date")
df0$Fund.ID <- as.factor(df0$Fund.ID)
rm(df.private.cfs)

# check if we have enough public data
min(df0$Vintage)
min(df0$Date)

if (public.filename %in% c("DebtFactorsUSD", "DebtFactorsEUR", "DebtFactorsEURUSD") ) {
  df0 <- df0[df0$Vintage > 2000, ]  
}
#df0 <- df0[df0$Vintage <= 2011, ]


# Export df0
if (export.data) {
  write.csv2(df0, paste0("data_private_public/df0_",
                         private.source, "_",
                         public.filename,"_",
                         weighting,".csv"), 
             row.names = FALSE)
}


# Summary statistics for df0
aggregate(Vintage ~ type , df0, min)
number.of <- function(x) length(table(x))
aggregate(Vintage ~ type , df0, number.of)
aggregate(as.character(Fund.ID) ~ type , df0, number.of)

# print summary: funds per vintage
df1 <- df0
df1 <- df1[!duplicated(df1$Fund.ID), ]
df1 <- as.data.frame.matrix(table(df1$Vintage, df1$type))
df1 <- df1[, colnames(df1) %in% c("BO", "DD", "INF", "MEZZ", "NATRES", "PD", "RE", "VC")]
df1["Total", ] <- as.integer(colSums(df1))
print(xtable::xtable(df1, caption = "Number of funds per vintage year.", label = "tab:pitchbook_data"), include.rownames = TRUE)
rm(df1)

## prepare pitchbook
# prep ----
library(readxl)

if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
if(!dir.exists("data_prepared")) dir.create("data_prepared")


path <- "data_in/pitchbook_fc_estimation_validation_2023-04-19.csv"
df.xl <- read.csv2(path)
colnames(df.xl)

df.xl$Fund.ID <- df.xl$FundID
df.xl$Vintage <- as.integer(df.xl$Vintage)
df.xl$date <- as.Date(df.xl$date)
cols.numeric <- c("FundSize", "Distributed_investor_filled", "Contributed_investor_filled", "NAV_investor_filled")
for (col in cols.numeric) df.xl[, col] <- as.numeric(df.xl[, col])

table(df.xl$region)
table(df.xl$asset_class[!duplicated(df.xl$Fund.ID)])
table(df.xl$asset_class_segment[!duplicated(df.xl$Fund.ID)])
table(df.xl$Vintage[!duplicated(df.xl$Fund.ID)])


make.strategy.summary <- function() {
  l <- list()
  for(asset.class in levels(as.factor(df.xl$asset_class))) {
    #print(asset.class)
    #print(table(df.xl$asset_class_segment[df.xl$ == asset.class]))
    df <- data.frame(table(df.xl$asset_class_segment[(df.xl$asset_class == asset.class) & (!duplicated(df.xl$Fund.ID))]))
    colnames(df) <- c("asset_class_segment", "Fund Count")
    df["Asset Class"] <- asset.class
    df <- df[, c(3,1,2)]
    l[[asset.class]] <- df
  }
  df <- data.frame(do.call(rbind, l))
  
  print(
    xtable::xtable(df, caption = "Summary of asset classes in Preqin dataset.", 
                   digits=3, label = "tab:summary_preqin_strategy"), 
    include.rownames=FALSE
    #, format.args = list(big.mark = ",")
  )
  return(df)
}
make.strategy.summary()



make.preqin.df <- function(
    df = df.xl,
    fund.size.weighting=TRUE,
    acs.filter = "PE_FOF",
    out.name = "FOF",
    region.filter = NA,
    vin.year.pfs = FALSE) {
  col.to.keep <- c("Fund.ID", "FundSize", "Vintage",
                   "date", 
                   "NAV_investor_filled", "Contributed_investor_filled", "Distributed_investor_filled"
                   )
  
  if(length(acs.filter) > 0) {
    
    test <- df[df$asset_class %in% acs.filter, ]
    if(nrow(test) > 0) df <- test
    
    test <- df[df$asset_class_segment %in% acs.filter, ]
    if(nrow(test) > 0) df <- test
    
  }
  
  if(!is.na(region.filter)) {
    df <- df[df$region %in% region.filter, ]
  }
  
  df <- df[, col.to.keep]
  
  for(col in c("Fund.ID")) {
    df[, col] <- as.factor(df[, col])
  }
  df$FundSize <- ifelse(is.na(df$FundSize), 
                        mean(df$FundSize, na.rm = TRUE), 
                        df$FundSize)
  df <- df[order(df$Fund.ID, df$date), ]
  
  list.df <- list()
  for(fund.id in levels(df$Fund.ID)) {
    df.ss <- df[df$Fund.ID == fund.id, ]
    
    # turn cumulative cash flows to normal cash flows
    df.ss$contributions <- c(df.ss$Contributed_investor_filled[1], diff(df.ss$Contributed_investor_filled))
    df.ss$distributions <- c(df.ss$Distributed_investor_filled[1], diff(df.ss$Distributed_investor_filled))
    df.ss$CF <- df.ss$distributions - df.ss$contributions
    
    # for last date: CF = CF + NAV, i.e., regard final NAV as final distribution
    df.ss$CF[nrow(df.ss)] <- df.ss$CF[nrow(df.ss)] + df.ss$NAV_investor_filled[nrow(df.ss)]
    
    # sanity check if cash flows start around the vintage
    min.year <- min(as.numeric(format(df.ss$date, "%Y")))
    vintage <- as.numeric(df.ss$Vintage[1])
    year.diff <-  min.year - vintage
    if (!is.na(year.diff)) {
      if(abs(year.diff) < 2) {
        list.df[[fund.id]] <- df.ss
      }
    }

  }
  df <- data.frame(do.call(rbind, list.df))
  
  if(fund.size.weighting) {
    # df$CF <- df$CF * df$FundSize # PREQIN
    df$CF <- df$CF # PITCHBOOK
  } else {
    # df$CF <- df$CF * mean(df$FundSize) # PREQIN
    df$CF <- df$CF / df$FundSize # PITCHBOOK
  }
  
  df$CF <- df$CF/ 1000 / 1000
  df <- df[, c("Fund.ID", "date", "CF", "Vintage")]
  colnames(df) <- c("Fund.ID", "Date", "CF", "Vintage")
  df$type <- out.name
  
  # Vintage Year Porftolio Formation
  agg.vin.year.pf <- function(df) {
    Vintage <- df$Vintage[1]
    type <- df$type[1]
    df.out <- aggregate(CF ~ Date, data = df, sum)
    df.out$type <- type
    df.out$Vintage <- Vintage
    df.out$Fund.ID <- paste0(Vintage, "_", type)
    df.out <- df.out[order(df.out$Date), ]
    return(df.out)
  }
  if (vin.year.pfs) {
    dfx <- split(df, df$Vintage)
    df <- data.frame(do.call(rbind,lapply(dfx, agg.vin.year.pf))) 
  }
  
  return(df)
}
df <- make.preqin.df(acs.filter = "PE_PD")
length(df$Fund.ID[!duplicated(df$Fund.ID)])

# run ----

acs <- list(
  INF = "INF", 
  NATRES = "NATRES", 
  PD = "PE_PD", 
  PE = "PE", 
  RE = "RE", 
  VC = "PE_VC",
  MEZZ = "PE_MEZZ",
  DD = "PE_DD",
  BO = "PE_BO",
  FOF = "PE_FOF",
  SEC = "PE_SEC",
  ALL = c("ALL", "INF", "NATRES", "OA", "PE", "RE")
)


make.preqin.csv <- function(fund.size.weighting, region.filter=NA, vin.year.pfs=TRUE) {
  l <- list()
  for(i in names(acs)) {
    print(acs[[i]])
    df0 = make.preqin.df(
      fund.size.weighting = fund.size.weighting,
      acs.filter = acs[[i]],
      region.filter = region.filter,
      vin.year.pfs = vin.year.pfs,
      out.name = i)
    l[[i]] <- df0
  }
  df.out <- data.frame(do.call(rbind, l))
  
  if(fund.size.weighting) {
    tag <- "FW"
  } else {
    tag <- "EW"
  }
  
  region.tag = ""
  if (!is.na(region.filter)) region.tag <- paste0("_", region.filter)
  if (vin.year.pfs) region.tag <- paste0("_VYP", region.tag)
    
  file = paste0("data_prepared/pitchbook_cashflows_", tag, region.tag, "_2023.csv")
  write.csv(df.out, file, row.names = FALSE)
  invisible(df.out)
}
vin.year.pfs <- FALSE
df.out <- make.preqin.csv(TRUE, "EU", vin.year.pfs)
df.out <- make.preqin.csv(FALSE, "EU", vin.year.pfs)
df.out <- make.preqin.csv(TRUE, NA, vin.year.pfs)
df.out <- make.preqin.csv(FALSE, NA, vin.year.pfs)

vin.year.pfs <- TRUE
df.out <- make.preqin.csv(TRUE, "EU", vin.year.pfs)
df.out <- make.preqin.csv(FALSE, "EU", vin.year.pfs)
df.out <- make.preqin.csv(TRUE, NA, vin.year.pfs)
df.out <- make.preqin.csv(FALSE, NA, vin.year.pfs)

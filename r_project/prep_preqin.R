## prepare preqin
# prep ----
library(readxl)

if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
if(!dir.exists("data_prepared")) dir.create("data_prepared")


path <- "data_in/Preqin_Cashflow_export-26_Feb_20e003d5aa-5a18-4c12-aba8-7586a7435ac9.xlsx"
sheet <- "Preqin_Export"
df.xl <- data.frame(readxl::read_excel(path = path, sheet = sheet))
table(df.xl$PRIMARY.GEOGRAPHIC.FOCUS) 
colnames(df.xl)

#df <- df.xl
table(df.xl$ASSET.CLASS[!duplicated(df.xl$FUND.ID)])
table(df.xl$STRATEGY[!duplicated(df.xl$FUND.ID)])
table(df.xl$VINTAGE...INCEPTION.YEAR[!duplicated(df.xl$FUND.ID)])




make.strategy.summary <- function() {
  l <- list()
  for(asset.class in levels(as.factor(df.xl$ASSET.CLASS))) {
    #print(asset.class)
    #print(table(df.xl$STRATEGY[df.xl$ASSET.CLASS == asset.class]))
    df <- data.frame(table(df.xl$STRATEGY[(df.xl$ASSET.CLASS == asset.class) & (!duplicated(df.xl$FUND.ID))]))
    colnames(df) <- c("Strategy", "Fund Count")
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
  acs.filter = "Fund of Funds",
  out.name = "FOF",
  region.filter = NA,
  vin.year.pfs = FALSE) {
  col.to.keep <- c("FUND.ID", "FUND.SIZE..USD.MN.", "VINTAGE...INCEPTION.YEAR",
                   "TRANSACTION.DATE", "TRANSACTION.AMOUNT", "NET.CASHFLOW")
  
  df <- df[df$TRANSACTION.TYPE == "Value", ]

  if(!is.na(acs.filter)) {
    
    test <- df[df$ASSET.CLASS %in% acs.filter, ]
    if(nrow(test) > 0) df <- test
    
    test <- df[df$STRATEGY %in% acs.filter, ]
    if(nrow(test) > 0) df <- test
    
  }
  
  if(!is.na(region.filter)) {
    df <- df[df$PRIMARY.GEOGRAPHIC.FOCUS %in% region.filter, ]
  }
  
  df <- df[, col.to.keep]
  
  for(col in c("FUND.ID")) {
    df[, col] <- as.factor(df[, col])
  }
  df$FUND.SIZE..USD.MN. <- ifelse(is.na(df$FUND.SIZE..USD.MN.), 
                                  mean(df$FUND.SIZE..USD.MN., na.rm = TRUE), 
                                  df$FUND.SIZE..USD.MN.)
  df <- df[order(df$FUND.ID, df$TRANSACTION.DATE), ]
  
  list.df <- list()
  for(fund.id in levels(df$FUND.ID)) {
    df.ss <- df[df$FUND.ID == fund.id, ]
    df.ss$CF <- c(df.ss$NET.CASHFLOW[1], diff(df.ss$NET.CASHFLOW))
    # for last date: CF = CF + NAV, i.e., regard final NAV as final distribution
    df.ss$CF[nrow(df.ss)] <- df.ss$CF[nrow(df.ss)] + df.ss$TRANSACTION.AMOUNT[nrow(df.ss)]
    
    year.diff <- min(as.numeric(format(df.ss$TRANSACTION.DATE, "%Y"))) - as.numeric(df.ss$VINTAGE...INCEPTION.YEAR[1])
    if(abs(year.diff) < 2) {
      list.df[[fund.id]] <- df.ss
    }
  }
  df <- data.frame(do.call(rbind, list.df))

  if(fund.size.weighting) {
    df$CF <- df$CF * df$FUND.SIZE..USD.MN.
  } else {
    df$CF <- df$CF * mean(df$FUND.SIZE..USD.MN.)
  }
  
  df$CF <- df$CF/ 1000 / 1000
  df <- df[, c("FUND.ID", "TRANSACTION.DATE", "CF", "VINTAGE...INCEPTION.YEAR")]
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
df <- make.preqin.df(acs.filter = "Natural Resources")
length(df$Fund.ID[!duplicated(df$Fund.ID)])

# run ----

acs <- list(
  INF = "Infrastructure", 
  NATRES = "Natural Resources", 
  PD = "Private Debt", 
  PE = "Private Equity", 
  RE = "Real Estate", 
  VC = "Venture Capital",
  MEZZ = "Mezzanine",
  DD = c("Special Situations", "Distressed Debt"),
  BO = c("Growth", "Buyout"),
  FOF = "Fund of Funds",
  SEC = "Secondaries"
)


make.preqin.csv <- function(fund.size.weighting) {
  l <- list()
  for(i in names(acs)) {
    print(acs[[i]])
    df0 = make.preqin.df(
      fund.size.weighting = fund.size.weighting,
      acs.filter = acs[[i]],
      region.filter = "Europe",
      vin.year.pfs = TRUE,
      out.name = i)
    l[[i]] <- df0
  }
  df.out <- data.frame(do.call(rbind, l))
  
  if(fund.size.weighting) {
    tag <- "FW"
  } else {
    tag <- "EW"
  }
  file = paste0("data_prepared/preqin_cashflows_", tag, "_Europe.csv")
  write.csv(df.out, file, row.names = FALSE)
  
}
make.preqin.csv(TRUE)
make.preqin.csv(FALSE)

## prepare preqin
# prep ----
if(sys.nframe() == 0L) rm(list = ls())

library(readxl)

data.prepared.folder <- "data_prepared_2026"


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
if(!dir.exists(data.prepared.folder)) dir.create(data.prepared.folder)


path <- "data_in/Preqin_Cashflow_export-26_Feb_20e003d5aa-5a18-4c12-aba8-7586a7435ac9.xlsx"
path <- "data_in/Preqin_Cashflow_export-14_Apr_22272f54e5-79fe-43c0-8bc0-9206381de0b7.xlsx"
sheet <- "Preqin_Export"
df.xl <- data.frame(readxl::read_excel(path = path, sheet = sheet))
if ("GEOGRAPHIC.FOCUS" %in% colnames(df.xl)) df.xl$PRIMARY.GEOGRAPHIC.FOCUS <- df.xl$GEOGRAPHIC.FOCUS
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

make.vintage.summary <- function(asset.class = "Private Equity") {
  df <- df.xl[df.xl$ASSET.CLASS == asset.class, ]
  df <- df[!duplicated(df$FUND.ID), ]
  
  # Create the frequency table
  tbl <- table(df$VINTAGE...INCEPTION.YEAR)
  
  # Convert table to the desired string format: "Year (Count)"
  # names(tbl) gets the years, and tbl itself contains the counts
  output_string <- paste(paste0(names(tbl), " (", tbl, ")"), collapse = ", ")
  
  print(paste0("The ", asset.class, " sample contains ", length(unique(df$FUND.ID)), 
               " distinct funds spreading over ",  length(tbl), " vintage years."))
  print(paste0("The region distribution is as follows: ", output_string, "."))
  
  # Region split
  tbl <- table(df$REGION)
  output_string <- paste(paste0(names(tbl), " (", tbl, ")"), collapse = ", ")
  print(paste0("The region distribution is as follows: ", output_string, "."))
}
make.vintage.summary("Private Equity")

make.preqin.df <- function(
  df = df.xl,
  fund.size.weighting=TRUE,
  acs.filter = "Fund of Funds",
  out.name = "FOF",
  region.filter = NA,
  vin.year.pfs = FALSE) {
  col.to.keep <- c("FUND.ID", "FUND.SIZE..USD.MN.", "VINTAGE...INCEPTION.YEAR",
                   "TRANSACTION.DATE", "TRANSACTION.AMOUNT", "NET.CASHFLOW",
                   "CUMULATIVE.CONTRIBUTION", "CUMULATIVE.DISTRIBUTION")
  
  df <- df[df$TRANSACTION.TYPE == "Value", ]

  if(length(acs.filter) > 0) {
    
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
    df.ss$CUMULATIVE.CONTRIBUTION <- c(df.ss$CUMULATIVE.CONTRIBUTION[1], diff(df.ss$CUMULATIVE.CONTRIBUTION))
    df.ss$CUMULATIVE.DISTRIBUTION <- c(df.ss$CUMULATIVE.DISTRIBUTION[1], diff(df.ss$CUMULATIVE.DISTRIBUTION))
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
    df$TRANSACTION.AMOUNT <- df$TRANSACTION.AMOUNT * df$FUND.SIZE..USD.MN. # aka NAV
    df$CUMULATIVE.CONTRIBUTION <- df$CUMULATIVE.CONTRIBUTION * df$FUND.SIZE..USD.MN.
    df$CUMULATIVE.DISTRIBUTION <- df$CUMULATIVE.DISTRIBUTION * df$FUND.SIZE..USD.MN.
  } else {
    # df$CF <- df$CF * mean(df$FUND.SIZE..USD.MN.) # why?
    df$CF <- df$CF
  }
  
  ten.mio <- 10 * 1000 * 1000
  df$CF <- df$CF/ ten.mio
  df$TRANSACTION.AMOUNT <- df$TRANSACTION.AMOUNT / ten.mio # aka NAV
  df$CUMULATIVE.CONTRIBUTION <- df$CUMULATIVE.CONTRIBUTION / ten.mio
  df$CUMULATIVE.DISTRIBUTION <- df$CUMULATIVE.DISTRIBUTION / ten.mio
  df <- df[, c("FUND.ID", "TRANSACTION.DATE", "CF", 
               "TRANSACTION.AMOUNT", "CUMULATIVE.CONTRIBUTION", "CUMULATIVE.DISTRIBUTION",
               "VINTAGE...INCEPTION.YEAR")]
  colnames(df) <- c("Fund.ID", "Date", "CF", "NAV", "CON", "DIS", "Vintage")
  df$type <- out.name
  
  # Vintage Year Porftolio Formation
  agg.vin.year.pf <- function(df) {
    Vintage <- df$Vintage[1]
    type <- df$type[1]
    df.out <- aggregate(CF ~ Date, data = df, mean)
    df.nav <- aggregate(NAV ~ Date, data = df, mean)
    df.out <- merge(df.out, df.nav, by = "Date")
    df.con <- aggregate(CON ~ Date, data = df, mean)
    df.out <- merge(df.out, df.con, by = "Date")
    df.dis <- aggregate(DIS ~ Date, data = df, mean)
    df.out <- merge(df.out, df.dis, by = "Date")
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
df <- make.preqin.df(acs.filter = "Private Equity")
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


make.preqin.csv <- function(fund.size.weighting, vin.year.pfs, region.filter) {
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
  
  tag <- ifelse(fund.size.weighting, "FW", "EW")
  tag <- ifelse(vin.year.pfs, paste0(tag, "_VYP"), tag)
  tag <- ifelse(is.na(region.filter), tag, paste0(tag, "_", region.filter))
  file = paste0(data.prepared.folder, "/preqin_cashflows_2022_", tag, "_NAV.csv")
  file = paste0(data.prepared.folder, "/preqin_cashflows_", tag, "_NAV.csv")
  write.csv(df.out, file, row.names = FALSE)
  invisible(df.out)
}
df.out <- make.preqin.csv(TRUE, TRUE, NA)
df.out <- make.preqin.csv(FALSE, TRUE, NA)
df.out <- make.preqin.csv(TRUE, FALSE, NA)
df.out <- make.preqin.csv(FALSE, FALSE, NA)
df.out <- make.preqin.csv(TRUE, TRUE, "Europe")
df.out <- make.preqin.csv(FALSE, TRUE, "Europe")

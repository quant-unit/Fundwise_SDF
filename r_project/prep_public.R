## prepare public market data
# load data -----
library(readxl)

if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
if(!dir.exists("data_prepared")) dir.create("data_prepared")

path <- "data_in/2020-02-26 150100 Update-51aec0-52aa6c.xlsx"
sheets <- c("msci_market_factors", "public_indices")
list.xl <- list()
for(sheet in sheets) {
  df.xl <- data.frame(readxl::read_excel(path = path, sheet = sheet))
  list.xl[[sheet]] <- df.xl
}

# public_indices ----
df <- list.xl[["public_indices"]]
colnames(df.xl)

# calc returns
for(col in colnames(df)[2:ncol(df)]) {
  df[, col] <- c(0, diff(df[, col]) / df[-1, col])
}

col.before <- colnames(df)

df$RF <- df$S.P.US.TREASURY.BILL
df$MKT <- df$MSCI.World.Net.Return.Daily - df$S.P.US.TREASURY.BILL
df$RE <- df$Wilshire.Global.Real.Estate.Securities.Total.Return.Index - df$MSCI.World.Net.Return.Daily
df$INF <- df$Thomson.Reuters.Global.Infrastructure.Total.Return.Index..USD. - df$MSCI.World.Net.Return.Daily
df$NATRES <- df$CS.Commodities.Benchmark.S.P.GSCI.Total.Return.Index - df$MSCI.World.Net.Return.Daily
df$VC <- df$NASDAQ.100.Total.Return - df$MSCI.World.Net.Return.Daily
df$PD <- df$iBoxx.USD.Liquid.Investment.Grade.Index.Index - df$MSCI.World.Net.Return.Daily
df$MEZZ <- df$iBoxx.USD.Liquid.High.Yield.Index - df$MSCI.World.Net.Return.Daily
df$DD <- df$iBoxx.USD.Liquid.High.Yield.Index - df$MSCI.World.Net.Return.Daily
df$BO <- df$Russell.2000.Total.Return.Index - df$MSCI.World.Net.Return.Daily
df$FOF <- df$Russell.2000.Total.Return.Index - df$MSCI.World.Net.Return.Daily
df$SEC <- df$Russell.2000.Total.Return.Index - df$MSCI.World.Net.Return.Daily
df$PE <- df$Russell.2000.Total.Return.Index - df$MSCI.World.Net.Return.Daily

df$Factor <- "Public"
  
df <- df[, c("Date", setdiff(colnames(df), col.before))]

write.csv(df, "data_prepared/public_returns.csv", row.names = FALSE)

# msci_market_factors ----

factor.df <- function(df = list.xl[["msci_market_factors"]]) {
  colnames(df)[1] <- "Date"
  
  cols <- grepl("World", colnames(df)) | grepl("Date", colnames(df))
  
  df <- df[, cols]
  for (col in colnames(df)) {
    if(col != "Date") {
      df[, col] <- as.numeric(df[, col])
      df[, col] <- as.numeric(ifelse(is.na(df[, col]), df$Market.World, df[, col]))
    }
  }
  col.before <- colnames(df)
  
  # Risk Free Rate
  df.RF <- read.csv("data_prepared/public_returns.csv")[, c("Date", "RF")]
  df <- merge(df, df.RF, by="Date")
  
  df$MKT <- df$Market.World - df$RF
  df$HML <- df$Value.World - df$Growth.World
  df$SMB <- df$Small.Cap.World - df$Large.Cap.World
  df$HDY <- df$High.Dividend.Yield.World - df$Market.World
  df$QLT <- df$Quality.World - df$Market.World
  df$MOM <- df$Momentum.World - df$Market.World
  df$ESG <- df$ESG.World - df$Market.World
  df$LOV <- df$Low.Volatility.World - df$Market.World
  
  df <- df[, c("Date", setdiff(colnames(df), col.before))]

  return(df)
}

df <- factor.df()

write.csv(df, "data_prepared/msci_market_factors.csv", row.names = FALSE)

# q-factor data ----
url <- "http://global-q.org/uploads/1/2/2/6/122679606/q5_factors_monthly_2019.csv"
url <- "data_in/q5_factors_monthly_2019.csv"
df.q <- read.csv(url)
cor(df.q[, 3:8])

cols.before <- colnames(df.q)
cols <- cols.before[3:ncol(df.q)]
for (col in cols) {
  new.col <- sub("R_", "", col)
  if (new.col == "F") new.col <- "RF"
  df.q[, new.col] <- df.q[, col] / 100
}
df.q$Date <- seq(as.Date("1967-02-01"),length=nrow(df.q),by="months")-1
df.q <- df.q[, setdiff(colnames(df.q), cols.before)]

apply(df.q[, 1:6], 2, function(x) prod(1+x))

write.csv(df.q, "data_prepared/q_factors.csv", row.names = FALSE)

# add future (yiels worse results)
if (FALSE) {
  df.q2 <- df.q
  df.q2$Date <- seq(as.Date("2020-02-01"), by = "month", length.out = nrow(df.q)) - 1
  df.q <- data.frame(rbind(df.q, df.q2))
  write.csv(df.q, "data_prepared/q_factors.csv", row.names = FALSE)
}



###### BOND FACTORS ------
# load Excel ------
library(readxl)
library(tidyr)

filename <- "data_in/Bond Indices iBoxx - large coverage.xlsx"
df.eur <- data.frame(readxl::read_excel(filename, sheet = "EuroData"))
df.usd <- data.frame(readxl::read_excel(filename, sheet = "USDData"))

# calculate returns -----------
make.return.df <- function(df) {
  # Extend date axis (to all daily dates)
  colnames(df)[colnames(df) == "Name"] <- "Date"
  df$Date <- as.Date(df$Date)
  min.date <- min(df$Date)
  # min.date <- as.Date("1969-12-31")
  df.date <- data.frame(Date = seq(min.date, max(df$Date),by ="day"))
  df <- merge(df, df.date, by = "Date", all.y=TRUE)
  
  # Forward fill all columns
  for (col in colnames(df)){
    df <-  df %>% fill(col, .direction = "down")
  }
  
  # select month end dates
  df$YearMonth <- format(df$Date, "%Y-%m")
  df <- df[!duplicated(df$YearMonth, fromLast=TRUE), ]
  
  # Return calculation
  calc.return <- function(x) {
    x <- as.numeric(x)
    y <- c(NA, x[-1] / x[1:(length(x)- 1)])
    # y <- y - 1
    return(y)
  }
  
  cols <- colnames(df)
  cols <- cols[!(cols %in% c("Date", "YearMonth"))]
  for (col in cols) {
    try(df[, col] <- calc.return(df[, col]))
  }
  
  # sort cols by number of missing values
  cols <- names(sort(apply(df, 2, function(x) sum(is.na(x)))))
  df <- df[, cols]
  
  ####
  return(df)
}

df.eur <- make.return.df(df.eur)
df.usd <- make.return.df(df.usd)

# calculate factors ------
colnames(df.usd)
apply(df.usd, 2, function(x) sum(is.na(x)))

fac.eur <- list(
  RF.long = "IBOXX.EURO.SOVEREIGNS.1.3...Tot..Rtn.Idx.Today",
  Term.long = "IBOXX.EURO.SOVEREIGNS.7.10...Tot..Rtn.Idx.Today",
  Corp.long = "IBOXX.EURO.CORPORATES.1.3...Tot..Rtn.Idx.Today",
  Corp.BBB.long = "IBOXX.EURO.CORPORATES.BBB.1.3...Tot..Rtn.Idx.Today",
  Liq.BBB.long =  "IBOXX.EURO.CORPORATES.BBB...Tot..Rtn.Idx.Today",
  Liq.BBB.short = "IBOXX.EURO.LIQUID.CORPORATES.BBB...Tot..Rtn.Idx.Today"
)

fac.usd <- list(
  RF.long = "IBOXX...SOVEREIGNS.1.3Y...Tot..Rtn.Idx.Today",
  Term.long = "IBOXX...SOVEREIGNS.7.10Y...Tot..Rtn.Idx.Today",
  Corp.long = "IBOXX...CORPORATES.1.3Y...Tot..Rtn.Idx.Today",
  Corp.BBB.long = "IBOXX...CORPORATES.BBB.1.3Y...Tot..Rtn.Idx.Today",
  Liq.BBB.long =  "IBOXX...CORPORATES.BBB...Tot..Rtn.Idx.Today" ,
  Liq.BBB.short = "IBOXX...LIQUID.INVESTMENT.GRADE.BBB.INDEX...Tot..Rtn.Idx.Today"
)


make.factors <- function(df, fac) {
  df$RF <- df[, fac$RF.long] - 1
  df$TERM <- df[, fac$Term.long] - df[, fac$RF.long]
  df$CORP <- df[, fac$Corp.long] - df[, fac$RF.long]
  df$HY <- df[, fac$Corp.BBB.long] - df[, fac$Corp.long]
  df$LIQ <- df[, fac$Liq.BBB.long] - df[, fac$Liq.BBB.short]
  
  rownames(df) <- df$Date
  cols <- c("RF", "TERM", "CORP", "HY", "LIQ")
  df <- df[, cols]
  
  return(df)
}

df.eur.fac <- make.factors(df.eur, fac.eur)
df.usd.fac <- make.factors(df.usd, fac.usd)

rownames2col <- function(df) {
  df <- cbind(Date = row.names(df), df)
  return(df)
}

write.csv(rownames2col(df.eur.fac), file = "data_prepared/DebtFactorsEUR.csv", row.names = FALSE)
write.csv(rownames2col(df.usd.fac), file = "data_prepared/DebtFactorsUSD.csv", row.names = FALSE)

# analyze factors ------
factor.summary <- function(df) {
  y <- apply(df, 2, function(x) {
    x[is.na(x)] <- 0
    y <- prod(1 + x)
    return(y)
  })
  
  return(y)
}

factor.summary(df.eur.fac)
factor.summary(df.usd.fac)

var(df.eur.fac, na.rm = TRUE)
var(df.usd.fac, na.rm = TRUE)

cor(df.eur.fac, use = "pairwise.complete.obs")
cor(df.usd.fac, use = "pairwise.complete.obs")


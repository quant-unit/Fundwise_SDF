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



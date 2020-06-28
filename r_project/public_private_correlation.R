# public private correlation
# load data -----
library(readxl)

if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
if(!dir.exists("data_prepared")) dir.create("data_prepared")

path <- "data_in/2020-02-26 150100 Update-51aec0-52aa6c.xlsx"
sheets <- c("msci_market_factors", "public_indices", "ca_index_100")
list.xl <- list()
for(sheet in sheets) {
  df.xl <- data.frame(readxl::read_excel(path = path, sheet = sheet))
  if("Name" %in% colnames(df.xl)) {
    df.xl$Date <- df.xl$Name
    df.xl$Name <- NULL
  }
  df.xl$Date <- as.Date(df.xl$Date)

  if (sheet == "public_indices") {
    for (col in colnames(df.xl)) {
      if (col != "Date") {
        df.xl[, col] <- c(NA, diff(df.xl[, col]) / df.xl[-nrow(df.xl), col])
      }
    }
    print(df.xl$Date)
  }

  for(col in colnames(df.xl)) {
    if (col != "Date") {
      df.xl[, col] <- as.numeric(as.character(df.xl[, col]))
    }
  }
  list.xl[[sheet]] <- df.xl
}

df <- list.xl[["msci_market_factors"]]
df <- Reduce(function(x,y) merge(x = x, y = y, by = "Date", all=TRUE),  list.xl)

# analyze ----

colnames(df)
cols <- c("MSCI.World.Net.Return.Daily", "Market.World", "Market.NA", "iBoxx.USD.Liquid.High.Yield.Index", "REXP.Total.Return")
df0 <- df[df$Date > as.Date("1998-12-31"), cols]
cor(df0, use="pairwise.complete.obs")

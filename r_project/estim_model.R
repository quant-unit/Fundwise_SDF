#### estimated model
# 1) load data -----

if(TRUE) {
  files <- c("prep_preqin.R", "prep_public.R")
  for(file in files) {
    source(file)
    if(sys.nframe() == 0L) rm(list = ls())
  }
}

if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

df.preqin <- read.csv("data_prepared/preqin_cashflows.csv")
public.filename <- "public_returns"
public.filename <- "msci_market_factors"
# public.filename <- "q_factors"
df.public <- read.csv(paste0("data_prepared/", public.filename, ".csv"))

df.preqin$Date <- as.Date(df.preqin$Date)
df.public$Date <- as.Date(df.public$Date)

df0 <- merge(df.preqin, df.public, by="Date")
df0$Fund.ID <- as.factor(df0$Fund.ID)
df0$Alpha <- 1

# 2) loss function ----
library(optimx)

npv.calc <- function(par, max.quarter = 60) {
  par0 <-  c("RF" = 1, par)
  
  npvs <- c()
  for(fund.id in levels(df$Fund.ID)) {
    df.ss <- df[df$Fund.ID == fund.id, ]
    df.ss$Return <- cumprod(1 + (as.matrix(df.ss[, names(par0)]) %*% par0))
  
    end <- min(max.quarter, nrow(df.ss))
    for(i in 1:end) {
      DCF <- df.ss$CF / df.ss$Return * df.ss$Return[i]
      npvs <- c(npvs, sum(DCF))
    }
    
  }
  
  y <- sum(npvs^2)
  return(y)
}

par <- c("MKT" = 1, "Alpha" = 0.05)
df <- df0[df0$type == "PD", ]
npv.calc(par)

# 3.0) rbind all ----
rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}

# 3) iterative run ----

max.quarters <- c(40)
types <- levels(df0$type)
#types <- c("DD")
factors <- colnames(df.public)[4:ncol(df.public)] # msci_factors
#factors <- colnames(df.public)[3:6] # q_factors
#factors <- c( "ALL", factors)

l <- list()
for (factor in factors) {
  for(max.quarter in max.quarters) {
    for (type in types) {
      print(paste(factor, max.quarter, type))
      par <- c("MKT" = 0, "Alpha"=0)
      
      if(factor == "ALL") {
        for(fac in setdiff(factors, "ALL")) {
          par[fac] <- 0
        }
      } else {
        par[factor] <- 0
        if ("ALL" %in% factors) break
      }
      
      df <- df0[df0$type == type, ]
      
      res <- optimx::optimx(par, npv.calc, max.quarter = max.quarter,
                            #lower = c("MKT" = -100, "RE" = -100, "Alpha" = -0.01), 
                            #upper = c("MKT" = 100, "RE" = 100, "Alpha" = 0.01),
                            method = "L-BFGS-B"
                            #method = c("Nelder-Mead", "L-BFGS-B", "nlminb", "nlm", "ucminf")
                            #control = list(all.methods=TRUE)
      )
      res$Factor <- factor
      res$Type <- type
      res$max.quarter <- max.quarter
      l[[paste0(factor, type)]] <- res
    }
  }

}

df.res <- data.frame(Reduce(rbind.all.columns, l))

file.out <- paste0("data_out/result_", public.filename, ".csv")
write.csv(df.res, file.out, row.names = FALSE)

df.res <- read.csv("data_out/result_msci_market_factors_60_2010.csv")
df.res1 <- df.res[order(df.res$value), ]
df.res1 <- df.res1[!(duplicated(df.res1$Type)), ]
df.res1

# 4) Grid Search ----

if (FALSE) {
  grid <- expand.grid(Alpha = seq(-0.005, 0.005, by=0.0025),
                      MKT = seq(-2, 2, 0.2),
                      RE = seq(-2, 2, 0.2))
  print(nrow(grid))
  
  
  val <- c()
  for(i in 1:nrow(grid)) {
    print(i)
    par <- as.numeric(grid[i, ])
    names(par) <- names(grid)
    
    val <- c(val, npv.calc(par))
  }
  grid <- data.frame(grid)
  grid$res <- val
  
  g1 <- grid
  
  head(g1[order(g1$res), ])
}



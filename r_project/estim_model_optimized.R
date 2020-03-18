#### estimated model
# 1) load data -----

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

weighting <- "FW"
df.preqin <- read.csv(paste0("data_prepared/preqin_cashflows_", weighting, ".csv"))
public.filename <- "public_returns"
public.filename <- "msci_market_factors"
public.filename <- "q_factors"
df.public <- read.csv(paste0("data_prepared/", public.filename, ".csv"))

df.preqin$Date <- as.Date(df.preqin$Date)
df.public$Date <- as.Date(df.public$Date)

df0 <- merge(df.preqin, df.public, by="Date")
df0$Fund.ID <- as.factor(df0$Fund.ID)
df0$Alpha <- 1

# 2.1) getNPVs function ----
library(Rcpp)

Rcpp::cppFunction("
    NumericVector getNPVs(NumericVector cf, NumericVector ret, int max_quarter){
    int n = cf.length();
    int i_max = std::min(max_quarter, n);
    NumericVector VectorOut(i_max, 0.0);

    for(int i = 0; i < i_max; ++i) {
      VectorOut[i] = sum(cf / ret * ret[i]);
    }
      return VectorOut;
    }
                      ")
getNPVs(df0[1:1000, 7], df0[1:1000, 5], 10)

# 2.2) err.sqr.calc function ----

error.function <- "log_L2_Lasso"

err.sqr.calc <- function(par, max.quarter = 60, lambda = 0, df) {
  par0 <-  c("RF" = 1, par)
  
  dfx <- split(df, df$Fund.ID)
  
  f1 <- function(df.ss) {
    Return <- exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0))))
    return(getNPVs(df.ss$CF, Return, max.quarter))
  }
  
  npvs <- lapply(dfx, f1)
  y <- sum(unlist(npvs)^2) / length(npvs)
  y <- log(y) + lambda * sum(abs(par[-1])) # lasso ignoring MKT
  #y <- log(sum(abs(unlist(npvs))) / length(npvs))
  return(y)
}

par <- c("MKT" = 1, "Alpha" = 0.05)
df.in <- df0[df0$type == "SEC", ]
system.time(
  print(err.sqr.calc(par, df=df.in, lambda = 0)) # 42181603
)

# 3.0) rbind all ----
rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}

# 3.1) cross validation -------

vintage.blocks <- list()
for(year in seq(1988, 2015, by=3)) {
  validate <- year:(year + 2)
  block <- (year - 3):(year - 1)
  block <- c(block, (year + 3):(year + 5))
  
  estimate <- 1980:2020
  estimate <- estimate[!(estimate %in% validate)]
  estimate <- estimate[!(estimate %in% block)]
  
  vintage.blocks[[paste(year)]] <- list(validate=validate, block=block, estimate=estimate)
}

for(year in names(vintage.blocks)) {
  for(key in c("estimate", "validate")) {
    vintage.blocks[[year]][[paste0("df.", key)]] <- df0[df0$Vintage %in% vintage.blocks[[year]][[key]], ]
    vintage.blocks[[year]][["cv.year"]] <- year
  }
}

viblo.list <- list()
for(block in vintage.blocks) {
  viblo.list[[paste(block$validate[1])]] <- data.frame(
    training.before = paste0("start-", block$block[1]-1),
    h.block.before = paste(block$block[1:3], collapse = ","), 
    validation = paste(block$validate, collapse = ","), 
    h.block.after =paste(block$block[4:6], collapse = ","),
    training.after = paste0(block$block[6]+1, "-end")
    
  )
}
print(xtable::xtable(data.frame(do.call(rbind, viblo.list))), include.rownames = FALSE)


vintage.blocks[["ALL"]][["df.estimate"]] <- df0
vintage.blocks[["ALL"]][["df.validate"]] <- df0
vintage.blocks[["ALL"]][["cv.year"]] <- "ALL"



# 3.2) iterative run ----
library(optimx)
library(parallel)
library(doParallel)
library(foreach)

iter.run <- function(input.list) {
  
  # 1) INITIALIZE
  df.estimate <- input.list$df.estimate
  df.validate <- input.list$df.validate
  
  lambdas <- seq(0,0.1,0.01)
  max.quarters <- c(40, 60)
  types <- levels(df0$type)
  types <- c("BO", "VC", "PE")
  if (public.filename == "msci_market_factors") {
    factors <- c("HML", "SMB", "HDY", "QLT", "MOM", "ESG", "LOV")
  }
  if (public.filename == "q_factors") {
    factors <- colnames(df.public)[3:6] # q_factors
    factors <- c( "ALL", factors)
  }
  
  # 2) RUNNER
  l <- list()
  for (lambda in lambdas) {
    for(max.quarter in max.quarters) {
      print(paste("Max.Q", max.quarter, "Lambda", lambda))
      
      cl <- parallel::makeForkCluster(3)
      doParallel::registerDoParallel(cl)
      
      output <- foreach (type = types, .combine = "rbind") %dopar% {
        factor.list <- list()
        for (factor in factors) {
          print(paste(type, factor))
          
          par <- c("MKT" = 1, "Alpha"=0)
          par <- c("MKT" = 1)
          
          if(factor == "ALL") {
            for(fac in setdiff(factors, "ALL")) {
              par[fac] <- 0
            }
          } else {
            par[factor] <- 0
            if ("ALL" %in% factors) break
          }
          
          if (TRUE) {
            res <- optimx::optimx(par, err.sqr.calc, 
                                  lambda = lambda,
                                  max.quarter = max.quarter,
                                  df = df.estimate[df.estimate$type == type, ],
                                  method = "nlminb"
                                  #method = c("Nelder-Mead", "L-BFGS-B", "nlminb", "nlm", "ucminf")
                                  #control = list(all.methods=TRUE)
            )
          } else {
            res <- data.frame(MKT=1, HML=1, SMB=1, HDY=1, QLT=1, MOM=1, ESG=1, LOV=1)
          }

          res$Factor <- factor
          res$Type <- type
          res$max.quarter <- max.quarter
          res$lambda <- lambda
          res$CV.key <- input.list$cv.year
          res$weighting <- weighting
          res$error.fun <- error.function
          res$datetime <- Sys.time()
          
          # validation error
          for(factor in names(par)) {
            par[factor] <- res[1, factor]
          }
          res$validation.error <- err.sqr.calc(par, max.quarter = max.quarter, lambda = lambda, df = df.validate[df.validate$type == type, ])
          
          factor.list[[factor]] <- res
        }
        return(Reduce(rbind.all.columns, factor.list))
      }
      output <- data.frame(output)
      print(output)
      l[[paste0(max.quarter, lambda)]] <- output
      
      parallel::stopCluster(cl)
    }
  }
  df.res <- data.frame(Reduce(rbind.all.columns, l))
  return(df.res)
}

#system.time(df.res <- iter.run(vintage.blocks$ALL))
#tag <- paste0("_50_2019_L1_", weighting)
#file.out <- paste0("data_out/result_", public.filename, tag, ".csv")
#write.csv(df.res, file.out, row.names = FALSE)

system.time(out.list <- lapply(vintage.blocks, iter.run))
saveRDS(out.list, "data_out/out_list_60_EW.RDS")

# x) Grid Search ----

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



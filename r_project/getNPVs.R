###############################
# getNPVs function & rbind.all
##############################
# 1) getNPVs function ----

use.cpp <- TRUE

if (use.cpp) {
  library(Rcpp)
  
  Rcpp::cppFunction("
                    double getNPVs(NumericVector cf, NumericVector ret, int max_month){
                    double DCF = sum(cf / ret);
                    
                    int n = cf.length();
                    int i_max = std::min(max_month, n);
                    NumericVector VectorOut(i_max, 0.0);
                    
                    for(int i = 0; i < i_max; i++) {
                    VectorOut[i] = DCF * ret[i];
                    }
                    
                    double y = mean(VectorOut);
                    return y;
                    }")
  n <- 13
  system.time(
    replicate(10000, getNPVs(rep(1,n), seq(1,n,1), n))
  )
  
} else {
  # use R function instead
  
  getNPVs <- function(cf, ret, max_month) {
    DCF = sum(cf / ret)
    
    n = length(cf)
    i_max = min(max_month, n)
    VectorOut = c()
    
    for(i in 0:i_max) {
      VectorOut <- c(VectorOut, DCF * ret[i])
    }
    y <- mean(VectorOut)
    return(y)
  }
  
  n <- 13
  system.time(
    replicate(10000, getNPVs(rep(1,n), seq(1,n,1), n))
  )
}

# 2) rbind.all ----
rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}
# 3) getNPVsDuration function ----

getNPVsDuration <- function(cf, ret, max_month) {
  
  duration <- sum(cf * 1:length(cf)) / sum(cf)
  floor.duration <- floor(duration)
    
  DCF = sum(cf / ret)
  
  n = length(cf)
  i_min = min(floor.duration, n)
  i_max = min(floor.duration + 1, n)
  VectorOut = c()
  
  for(i in i_min:i_max) {
    VectorOut <- c(VectorOut, DCF * ret[i])
  }
  y <- mean(VectorOut)
  return(y)
}

n <- 13
system.time(
  replicate(1, getNPVsDuration(rep(1,n), seq(1,n,1), n))
)

# 4) getNPVsSingleDate function ----

getNPVsSingleDate <- function(cf, ret, max_month) {
  
  DCF = sum(cf / ret)
  
  n = length(cf)
  i_max = min(max_month, n)
  
  y <- DCF * ret[i_max]
  return(y)
}

n <- 13
system.time(
  replicate(1, getNPVsSingleDate(rep(1,n), seq(1,n,1), n))
)


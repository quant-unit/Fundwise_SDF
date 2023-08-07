###############################
# getNPVs function & rbind.all
##############################
# 1) getNPVs function ----

use.cpp <- FALSE

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
##################################
### HAC estimator implemented in R
### Heberle Sattarhoff - A Fast Algorithm for the Computation of HAC Covariance Matrix Estimators (2017)
### https://doi.org/10.3390/econometrics5010009
# algorithmns -----
require(fftwtools)

HAC.new <- function ( mcond , method , bw ) {
   dimmcond <- dim ( mcond )
   Nlen <- dimmcond [1]
   qlen <- dimmcond [2]
   ww <- kweightsHAC ( kernel = method , Nlen , bw )
   ww <- c (1 , ww [1:( Nlen -1) ] , 0 , ww [( Nlen -1) :1])
   ww <- Re ( fftw ( ww ) )
   FF <- rbind ( mcond , matrix (0 , Nlen , qlen ) )
   FF <- mvfftw ( FF )
   FF <- FF * matrix ( rep ( ww , qlen ) , ncol = qlen )
   FF <- Re ( mvfftw ( FF , inverse = TRUE ) ) / (2 * Nlen )
   FF <- FF [1: Nlen ,]
   return (( t ( mcond ) %*% FF ) / Nlen )
   }

HAC.ron <- function ( mcond , method , bw ) {
   dimmcond <- dim ( mcond )
   Nlen <- dimmcond [1]
   qlen <- dimmcond [2]
   ww <- kweightsHAC ( kernel = method , Nlen , bw )
   LL <- ( crossprod ( mcond , mcond ) ) / Nlen
   for ( i in 1: bw ) {
     GG <- rbind ( matrix (0 ,i , qlen ) , mcond [1:( Nlen - i ) ,])
     GG <- ( crossprod ( mcond , GG ) ) / Nlen
     LL <- LL + ww [ i ] * ( GG + t ( GG ) )
     }
   return ( LL )
}

HAC.zei <- function ( mcond , method , bw ) {
   Nlen <- dim ( mcond ) [1]
   ww <- kweightsHAC ( kernel = method , Nlen , bw )
   LL <- crossprod ( mcond , mcond ) / Nlen
   for ( i in 1: bw ) {
     GG <- ( crossprod ( mcond [( i +1) : Nlen ,] , mcond [1:( Nlen - i ) ,]) ) / Nlen
     LL <- LL + ww [ i ] * ( GG + t ( GG ) )
     }
   return ( LL )
   }

HAC.kyr <- function ( mcond , method , bw ) {
   Nlen <- dim ( mcond ) [1]
   ww <- kweightsHAC ( kernel = method , Nlen , bw )
   ww <- c (1 , ww [1:( Nlen -1) ])
   TT <- matrix (0 , Nlen , Nlen )
   TT <- matrix ( ww [ abs ( col ( TT ) - row ( TT ) ) + 1] , Nlen , Nlen )
   return ( crossprod ( t ( crossprod ( mcond , TT )) , mcond ) / Nlen )
   }

kweightsHAC <- function (kernel = c ("Truncated" , "Bartlett" , "Parzen", "Tukey - Hanning", "Quadratic Spectral"),
                          dimN , bw ) {
  ww <- numeric ( dimN )
  switch ( kernel ,
             Truncated = {
               ww [1: bw ] <- 1
               } ,
             Bartlett = {
               ww [1: bw ] <- 1 - ( seq (1 , bw ) / ( bw +1) )
               } , Parzen = {
                 seq1 <- ( seq (1 , floor ( bw / 2) ) ) / bw
                 seq2 <- ( seq ( floor ( bw / 2) +1 , bw ) ) / bw
                 ww [1: length ( seq1 ) ] <- 1 - 6 * seq1 ^2 + 6 * seq1 ^3
                 ww [( length ( seq1 ) +1) : bw ] <- 2 * (1 - seq2 ) ^3
                 } ,  "Tukey - Hanning" = {
                   ww [1: bw ] <- (1 + cos ( pi * (( seq (1 , bw ) ) / bw ) ) ) / 2
                   } , "Quadratic Spectral" = {
                     aa <- pi * (( seq (1 , dimN ) ) / bw ) / 5
                     ww <- 1 / (12 * aa ^2) * ( sin (6 * aa ) / (6 * aa ) - cos (6 * aa ) )
                     })
  return ( ww )
  }


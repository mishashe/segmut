#' Chi statistics
#'
#' @description calculates sum of Chi statistics given breakpoints and penalize.
#'
#' @param par vector of breakpoints
#' @param muts vector of mutations locations
#' @param L length of locus
#' @param Kmin minimal length of a segment
#'
#' @return penalized chi square
#' @export
#'
getChiSquare <- function(par,muts,L=max(muts)-min(muts)+1,Kmin=0)
{
  averageDivergence <- length(muts)/L
  xB <- c(0,sort(par),L)
  Ks <- diff(xB)
  if (any(Ks<Kmin)) return(Inf)
  D <- 0
  for (i in 1:(length(xB)-1))
  {
    nmuts <- sum(muts>=xB[i] & muts<=xB[i+1])
    D <- D  + (nmuts - Ks[i]*averageDivergence)^2/(Ks[i]*averageDivergence)
  }
  return(length(par)*log(length(muts))-D)
}



#' finds optimal breaks
#'
#' @description fits optimal breaks given number of breaks
#'
#' @param muts vector of mutations locations
#' @param L length of locus
#' @param Kmin minimal length of a segment
#' @param n number of breaks
#'
#' @return results of the DE optimization
#' @export
#'
getBreaksChiSquare <- function(muts,L=max(muts)-min(muts)+1,Kmin=0,n=1)
{
  res <- DEoptim(getChiSquare, lower = rep(0,n), upper = rep(L,n),
                 control = DEoptim.control(trace = 0),
                 muts=muts,L=L,Kmin=Kmin)
  return(res)
}



#' finds optimal number of breaks
#'
#' @description fits optimal breaks NOT given number of breaks
#'
#' @param muts vector of mutations locations
#' @param L length of locus
#' @param Kmin minimal length of a segment
#' @param pThreshold significance threshold of the uniform KS test
#'
#' @return breaks locations for optimal number of breaks
#' @export
#'
getNumberBreaksChiSquare <- function(muts,L=max(muts)-min(muts)+1,Kmin=0)
{
  n <- 0
  parPrev <- c()
  ChiSquarePrev <- getChiSquare(parPrev,muts,L=L,Kmin=Kmin)
  repeat
  {
    n <- n + 1
    par <- getBreaksChiSquare(muts = muts, L = L, Kmin=Kmin, n=n)$optim$bestmem
    ChiSquare <- getChiSquare(par,muts,L=L,Kmin=Kmin)
    if(ChiSquare>ChiSquarePrev){return(parPrev)}
    parPrev <- par
    ChiSquarePrev <- ChiSquare
  }
}




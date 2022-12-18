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
# Data driven chi-square test for
# uniformity with unequal cells
# Tadeusz Inglot a
# & Alicja Janic-Wróblewska
getChiSquareGreedy <- function(par,parPrev,muts,L=max(muts)-min(muts)+1,Kmin=0)
{
  # par <- sort(par)
  xB <- sort(c(0,c(par,parPrev),L+1))
  Ks <- diff(xB) + 1
  if (any(Ks<Kmin)) return(Inf)
  expected <- (length(muts)/L)*Ks
  nmuts <- .Call(graphics:::C_BinCount, muts, xB, TRUE, TRUE) # a bit faster than nmuts <- hist(muts,breaks=xB,plot=FALSE)$counts
  chiSquare <- sum((nmuts - expected)^2/expected)
  return(-chiSquare)
}


#' finds optimal breaks ChiSquare
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
getBreaksChiSquareGreedy <- function(parPrev,muts,L=max(muts)-min(muts)+1,Kmin=0,n=1)
{
  res <- DEoptim(getChiSquareGreedy, lower = rep(Kmin,1), upper = rep(L-Kmin,1),
                 control = DEoptim.control(trace = FALSE,
                                           parallelType='none',
                                           strategy = 1,
                                           NP=100,
                                           storepopfreq = 100000,
                                           reltol=1e-7),
                 parPrev,muts=muts,L=L,Kmin=Kmin)$optim$bestmem
  return(res)
}



#' finds optimal number of breaks ChiSquare
#'
#' @description fits optimal breaks NOT given number of breaks
#'
#' @param muts vector of mutations locations
#' @param L length of locus
#' @param Kmin minimal length of a segment
#'
#' @return breaks locations for optimal number of breaks
#' @export
#'
getNumberBreaksChiSquareGreedy <- function(muts,L=max(muts)-min(muts)+1,Kmin=0, nmax=0)
{
  n <- 0
  parPrev <- c()
  ChiSquareVector <- c(getChiSquareGreedy(c(),c(),muts,L=L,Kmin=Kmin))
  if (length(muts)<2 | (Kmin*2+2)>=L) return(list(parPrev,ChiSquareVector))
  ChiSquarePrev <- 0*log(length(muts)) + ChiSquareVector[1]
  repeat
  {
    n <- n + 1
    par <- getBreaksChiSquareGreedy(parPrev=parPrev,muts = muts, L = L, Kmin=Kmin, n=n)
    ChiSquare <- n*log(length(muts)) + getChiSquareGreedy(par, parPrev,muts, L=L, Kmin=Kmin)
    ChiSquareVector <- c(ChiSquareVector,ChiSquare)
    if(ChiSquare>=ChiSquarePrev & n>=nmax){return(list(parPrev,ChiSquareVector))}
    if(Kmin*(n+1)>=L){return(list(parPrev,ChiSquareVector))}
    parPrev <- c(parPrev,par)
    ChiSquarePrev <- ChiSquare
  }
}

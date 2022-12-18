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
getChiSquare <- function(par,muts,L=max(muts)-min(muts)+1,Kmin=0,parPrev=c())
{
  xB <- c(0,sort(c(par,parPrev)),L+1)
  Ks <- diff(xB) + 1
  if (any(Ks<Kmin)) return(-Inf)
  expected <- (length(muts)/L)*Ks
  nmuts <- .Call(graphics:::C_BinCount, muts, xB, TRUE, TRUE) # a bit faster than nmuts <- hist(muts,breaks=xB,plot=FALSE)$counts
  chiSquare <- sum((nmuts - expected)^2/expected)
  return(chiSquare)
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
getBreaksChiSquare <- function(muts,L=max(muts)-min(muts)+1,Kmin=0,n=1,parPrev=c())
{
  res <- ga(type = "real-valued", fitness = getChiSquare,optim=FALSE,keepBest=FALSE,
            lower = rep(Kmin,1),
            upper = rep(L-Kmin,1),
            parallel=FALSE,
            monitor = FALSE,
            muts=muts,L=L,Kmin=Kmin,parPrev=parPrev)@solution
  res <- sort(c(parPrev,res))
  suggestions <- matrix(res,nrow=1)
  res <- ga(type = "real-valued", fitness = getChiSquare,optim=FALSE,keepBest=FALSE,
            lower = res-diff(c(0,res))/2,
            upper = res+diff(c(res,L))/2,
            parallel=FALSE,
            suggestions=suggestions,monitor = FALSE,
           muts=muts,L=L,Kmin=Kmin,parPrev=c())@solution
  return(sort(c(res)))
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
getNumberBreaksChiSquare <- function(muts,L=max(muts)-min(muts)+1,Kmin=0, nmax=0)
{
  n <- 0
  parPrev <- c()
  ChiSquareVector <- c(-getChiSquare(parPrev,muts,L=L,Kmin=Kmin))
  parBest <- parPrev
  if (length(muts)<2 | (Kmin*2+2)>=L) return((parBest))
  ChiSquarePrev <- 0*log(length(muts)) - ChiSquareVector[1]
  repeat
  {
    n <- n + 1
    par <- getBreaksChiSquare(muts = muts, L = L, Kmin=Kmin, n=n, parPrev)
    ChiSquare <- n*log(length(muts)) - getChiSquare(par, muts, L=L, Kmin=Kmin)
    ChiSquareVector <- c(ChiSquareVector,ChiSquare)
    if(ChiSquare<ChiSquarePrev) parBest <- par
    if(ChiSquare>=ChiSquarePrev & ChiSquare>=ChiSquarePrev){return((parBest))}
    if(Kmin*(n+1)>=L){return((parBest))}
    parPrev <- par
    ChiSquarePrev <- ChiSquare
  }
}


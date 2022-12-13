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
getChiSquare <- function(par,muts,L=max(muts)-min(muts)+1,Kmin=0)
{
  # par <- sort(par)
  xB <- c(0,par,L)
  Ks <- diff(xB)
  if (any(Ks<Kmin)) return(10)
  expected <- (length(muts)/L)*Ks
  # nmuts <- .Internal(tabulate(.Internal(findInterval(vec=xB,x=muts, FALSE,FALSE,FALSE)), nbins=length(xB)-1))
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
getBreaksChiSquare <- function(muts,L=max(muts)-min(muts)+1,Kmin=0,n=1)
{
  initialpop <- Kmin + 1 + matrix(runif(10*1*n)*(L-Kmin-1),ncol=n)
  if (n>1) initialpop=t(apply(initialpop,1,sort))
  res <- DEoptim(getChiSquare, lower = rep(Kmin,n), upper = rep(L-Kmin,n),
                 control = DEoptim.control(trace = FALSE,
                                           parallelType='none',
                                           strategy = 1,
                                           initialpop=initialpop,
                                           storepopfreq = 10000),
                 muts=muts,L=L,Kmin=Kmin)$optim$bestmem
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
getNumberBreaksChiSquare <- function(muts,L=max(muts)-min(muts)+1,Kmin=0)
{
  n <- 0
  parPrev <- c()
  if (length(muts)<2 | Kmin>=L) return(parPrev)
  ChiSquarePrev <- 0*log(length(muts)) + getChiSquare(parPrev,muts,L=L,Kmin=Kmin)
  repeat
  {
    n <- n + 1
    par <- getBreaksChiSquare(muts = muts, L = L, Kmin=Kmin, n=n)
    ChiSquare <- n*log(length(muts)) + getChiSquare(par,muts,L=L,Kmin=Kmin)
    if(ChiSquare>ChiSquarePrev | Kmin*(n+1)>=L){return(parPrev)}
    parPrev <- par
    ChiSquarePrev <- ChiSquare
  }
}




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
  par <- sort(par)
  xB <- c(0,par,L)
  Ks <- diff(xB)
  if (any(Ks<Kmin)) return(Inf)
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
  res <- DEoptim(getChiSquare, lower = seq(1,n,1)*Kmin, upper = L-seq(n,1,-1)*Kmin,
                 control = DEoptim.control(trace = 0,parallelType='none'),
                 muts=muts,L=L,Kmin=Kmin)
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
  ChiSquarePrev <- 0*log(length(muts)) + getChiSquare(parPrev,muts,L=L,Kmin=Kmin)
  if (Kmin>=L) return(parPrev)
  repeat
  {
    n <- n + 1
    par <- getBreaksChiSquare(muts = muts, L = L, Kmin=Kmin, n=n)$optim$bestmem
    ChiSquare <- n*log(length(muts)) + getChiSquare(par,muts,L=L,Kmin=Kmin)
    if(ChiSquare>ChiSquarePrev | Kmin*(n+1)>=L){return(parPrev)}
    parPrev <- par
    ChiSquarePrev <- ChiSquare
  }
}




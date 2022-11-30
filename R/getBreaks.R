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
getBreaks <- function(muts,L=max(muts)-min(muts)+1,Kmin=0,n=1)
{
  res <- DEoptim(D_KS, lower = rep(0,n), upper = rep(L,n),
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
getNumberBreaks <- function(muts,L=max(muts)-min(muts)+1,Kmin=0,pThreshold=0.05)
{
  p <- min(p_KS(c(),muts,L=L,Kmin=Kmin))
  if (p>pThreshold) return(c())
  n <- 1
  repeat
  {
    res <- getBreaks(muts = muts, L = L, Kmin=Kmin, n=n)
    p <- min(p_KS(res$optim$bestmem,muts,L=L,Kmin=Kmin))
    if(min(p)>pThreshold){return(res$optim$bestmem)}
    n <- n + 1
  }
}




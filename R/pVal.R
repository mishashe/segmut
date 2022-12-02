#' KS statistics
#'
#' @description calculates sum of KS statistics given breakpoints.
#'
#' @param par vector of breakpoints
#' @param muts vector of mutations locations
#' @param L length of locus
#' @param Kmin minimal length of a segment
#'
#' @return sum of KS statistics
#' @export
#'
#' @examples
#' D_KS(par=c(3000,5000),muts =sort(c(sample(1:3000,400),sample(3001:8000,250),sample(8001:10000,400))),L=10000,Kmin=0)
D_KS <- function(par,muts,L=max(muts)-min(muts)+1,Kmin=0)
{
  xB <- c(0,sort(par),L)
  if (any(diff(xB)<Kmin)) return(Inf)
  D <- 0
  for (i in 1:(length(xB)-1))
  {
    Ind <- which(muts>=xB[i] & muts<=xB[i+1])
    if (length(Ind)>3)
    {
      tes <- ks.test(muts[Ind],"punif",xB[i],xB[i+1],exact=FALSE)
      D <- D  + 1/abs(tes$p.value) + abs(tes$statistic)
    }
  }
  return(D)
}


#' KS p values
#'
#' @description calculates vector of KS p-values given breakpoints.
#'
#' @param par vector of breakpoints
#' @param muts vector of mutations locations
#' @param L length of locus
#' @param Kmin minimal length of a segment
#'
#' @return vector of KS p-values
#' @export
#'
#' @examples
#' p_KS(par=c(3000,5000),muts =sort(c(sample(1:3000,400),sample(3001:8000,250),sample(8001:10000,400))),L=10000,Kmin=0)

p_KS <- function(par,muts,L=max(muts)-min(muts)+1,Kmin=0)
{
  xB <- c(0,sort(par),L)
  if (any(diff(xB)<Kmin)) return(c(1))
  p <- c()
  for (i in 1:(length(xB)-1))
  {
    Ind <- which(muts>=xB[i] & muts<=xB[i+1])
    if (length(Ind)>3)
    {
      p <- c(p,ks.test(muts[Ind],"punif",xB[i],xB[i+1],exact=FALSE)$p.value)
    }
  }
  return(p)
}



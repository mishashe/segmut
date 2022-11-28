#' calculates geometrically averaged p-value
#'
#' @description calculates geometrically averaged p-value given breakpoints.
#'
#' @param par vector of breakpoints
#' @param muts vector of mutations locations
#' @param L length of locus
#' @param Kmin minimal length of a segment
#'
#' @return p-value of all breakpoints after harmonic mean
#' @export
#'
#' @examples
#' pVal(par=c(3000,5000),muts =sort(c(sample(1:3000,400),sample(3001:8000,250),sample(8001:10000,400))),L=10000,Kmin=0)
pVal <- function(par,muts,L=max(muts)-min(muts)+1,Kmin=0)
{
  xB <- c(0,sort(par),L)
  if (any(diff(xB)<Kmin)) return(0)
  pVal <- 0
  n1 <- sum(muts>=xB[1] & muts<=xB[2])
  K1 <- xB[2]-xB[1]
  for (i in 2:(length(xB)-1))
  {
    n2 <- sum(muts>=xB[i] & muts<=xB[i+1])
    K2 <- xB[i+1]-xB[i]
    if (n1<10 & n2<10) return(0)
    tau <- (n1+n2)/(K1+K2)
    pVal <- pVal + stats::pnorm(abs(n1/K1-n2/K2), mean = 0, sd = sqrt(tau/K1+tau/K2), lower.tail=FALSE, log.p = TRUE)
    n1 <- n2
    K1 <- K2
  }
  return(pVal/length(par))
}



#' finds optimal breaks
#'
#' @description fits optimal breaks given number of breaks
#'
#' @param muts vector of mutations locations
#' @param L length of locus
#' @param Kmin minimal length of a segment
#'
#' @return results of the DE optimization
#' @export
getBreaks <- function(muts,L=max(muts)-min(muts)+1,Kmin=0,n)
{
  res <- DEoptim(pVal, lower = rep(0,n), upper = rep(L,n),
                 control = DEoptim.control(trace = 0),
                 muts=muts,L=L)
  return(res)
}







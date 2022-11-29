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
  res <- DEoptim(pVal, lower = rep(0,n), upper = rep(L,n),
                 control = DEoptim.control(trace = 0),
                 muts=muts,L=L,Kmin=Kmin)
  return(res)
}







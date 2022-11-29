#' p-value of the null model
#'
#' @description estimates p-value of the null model
#'
#' @param muts vector of mutations locations
#' @param L length of locus
#' @param Kmin minimal length of a segment
#' @param A number of realizations
#'
#' @return p-value of the null model
#' @export
#'
get_p0 <- function(muts,L=max(muts)-min(muts)+1,Kmin=0,A=20)
{
  p.values <- foreach(a=1:A, .combine='c') %dopar%
    {return(getBreaks(muts = runif(length(muts), 0, L), L = L, Kmin=0, n=1)$optim$bestval)}
  p0 <- quantile(p.values,0.05)
  return(p0)
}





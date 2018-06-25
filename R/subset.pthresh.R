#' @title Subsetting a matrix of p-value threshold betas
#' @description Subsetting a matrix of p-value threshold betas
#' @export
`[.pthresh` <- function(obj, row=1:length(obj$beta), col) {
  stopifnot(!missing(col))
  stopifnot(col <= obj$nbin)
  beta <- obj$beta[row] * 0
  include <- obj$pbin[row] <= col
  beta[row][include] <- obj$beta[row][include]
  return(beta)
}

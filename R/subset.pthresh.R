#' @title Subsetting a matrix of p-value threshold betas
#' @description Subsetting a matrix of p-value threshold betas
#' @export
`[.pthresh` <- function(obj, row=1:length(obj$beta), col=1:obj$nbin) {
  if(is.logical(col)) col <- which(col)
  if(length(col) > 1) {
    return(do.call(cbind, lapply(col, function(i) obj[row,i])))
  }
  stopifnot(col <= obj$nbin)
  beta <- obj$beta[row] * 0
  include <- obj$pbin[row] <= col
  beta[row][include] <- obj$beta[row][include]
  return(matrix(beta, ncol=1))
}

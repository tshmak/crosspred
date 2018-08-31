#' @title Subsetting a matrix of p-value threshold betas
#' @description Subsetting a matrix of p-value threshold betas
#' @export
`[.pthresh` <- function(obj, row=1:length(obj$beta), col=1:obj$nbin, pthresh=FALSE) {
  if(is.logical(col)) col <- which(col)
  
  if(pthresh) {
    obj2 <- obj
    obj2$beta <- obj2$beta[row]
    obj2$pbin <- obj2$pbin[row]
    obj2$nbin <- length(col)
    obj2$p.thresholds <- obj2$p.thresholds[col]
    class(obj2) <- "pthresh"
    return(obj2)
  } else {
    if(length(col) > 1) {
      return(do.call(cbind, lapply(col, function(i) obj[row,i])))
    }
    stopifnot(col <= obj$nbin)
    beta <- obj$beta[row] * 0
    include <- obj$pbin[row] <= col
    beta[row][include] <- obj$beta[row][include]
    return(matrix(beta, ncol=1))
  }
}

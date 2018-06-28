cor2p <- function(cor, n) {
  #' Reverse of p2cor
  t <- abs(cor) * sqrt((n - 2)/(1-cor^2))
  p <- pt(t, df=n-2, lower.tail=F) * 2
  return(p)
}
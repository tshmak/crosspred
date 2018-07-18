subset.lassosum.pipeline <- function(lassosum.pipeline, validate.lassosum, s, lambda) {
  #' subsetting a lassosum.pipeline object, by s and lambda. Alternatively, 
  #' provide a validate.lassosum object. 
  if(!missing(validate.lassosum)) {
    s <- validate.lassosum$best.s
    lambda <- validate.lassosum$best.lambda
  } else {
    stopifnot(!(missing(s) || missing(lambda)) )
  }
  
  lp <- lassosum.pipeline
  select.s <- lp$s %in% s
  select.lambda <- lp$lambda %in% lambda
  lp$lambda <- lp$lambda[select.lambda]
  for(i in 1:length(lp$beta)) {
    lp$beta[[i]] <- lp$beta[[i]][,select.lambda]
    if(!is.null(lp$pgs)) {
      lp$pgs[[i]] <- lp$pgs[[i]][,select.lambda] 
    }
  }
  lp$beta <- lp$beta[select.s]
  lp$pgs <- lp$pgs[select.s]
  
  return(lp)
}
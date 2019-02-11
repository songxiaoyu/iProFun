trigammaInverse <- function(x) {
  #	Solve trigamma(y) = x for y
  #	Gordon Smyth
  #	8 Sept 2002.  Last revised 12 March 2004.

  #	Non-numeric or zero length input
  if(!is.numeric(x)) stop("Non-numeric argument to mathematical function")
  if(length(x)==0) return(numeric(0))

  #	Treat out-of-range values as special cases
  omit <- is.na(x)
  if(any(omit)) {
    y <- x
    if(any(!omit)) y[!omit] <- Recall(x[!omit])
    return(y)
  }
  omit <- (x < 0)
  if(any(omit)) {
    y <- x
    y[omit] <- NaN
    warning("NaNs produced")
    if(any(!omit)) y[!omit] <- Recall(x[!omit])
    return(y)
  }
  omit <- (x > 1e7)
  if(any(omit)) {
    y <- x
    y[omit] <- 1/sqrt(x[omit])
    if(any(!omit)) y[!omit] <- Recall(x[!omit])
    return(y)
  }
  omit <- (x < 1e-6)
  if(any(omit)) {
    y <- x
    y[omit] <- 1/x[omit]
    if(any(!omit)) y[!omit] <- Recall(x[!omit])
    return(y)
  }

  #	Newton's method
  #	1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
  #	so iteration to solve 1/x = 1/trigamma is monotonically convergent
  y <- 0.5+1/x
  iter <- 0
  repeat {
    iter <- iter+1
    tri <- trigamma(y)
    dif <- tri*(1-tri/x)/psigamma(y,deriv=2)
    y <- y+dif
    if(max(-dif/y) < 1e-8) break
    if(iter > 50) {
      warning("Iteration limit exceeded")
      break
    }
  }
  y
}

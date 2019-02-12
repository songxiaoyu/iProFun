##fitFDist is the function to obtain the empirical bayes estimate of n0 and s0^2 described on the top of the file. This function is extracted from limma package.
##The input x is the vector of s_g_square's (all trios in one tissue), and df1 is the vector of df_g's (all trios in one tissue).
fitFDist <- function(x,df1,covariate=NULL)
  #	Moment estimation of the parameters of a scaled F-distribution.
  #	The numerator degrees of freedom are given, the denominator is to be estimated.
  #	Gordon Smyth and Belinda Phipson
  #	8 Sept 2002.  Last revised 25 Jan 2017.
{
  #	Check x
  n <- length(x)
  if(n <= 1L) return(list(scale=NA,df2=NA))

  #	Check df1
  ok <- is.finite(df1) & df1 > 1e-15
  if(length(df1)==1L) {
    if(!ok) {
      return(list(scale=NA,df2=NA))
    } else {
      ok <- rep_len(TRUE,n)
    }
  } else {
    if(length(df1) != n) stop("x and df1 have different lengths")
  }

  #	Check covariate
  if(is.null(covariate)) {
    splinedf <- 1L
  } else {
    if(length(covariate) != n) stop("x and covariate must be of same length")
    if(anyNA(covariate)) stop("NA covariate values not allowed")
    isfin <- is.finite(covariate)
    if(!all(isfin)) {
      if(any(isfin)) {
        r <- range(covariate[isfin])
        covariate[covariate == -Inf] <- r[1]-1
        covariate[covariate == Inf] <- r[2]+1
      } else {
        covariate <- sign(covariate)
      }
    }
    splinedf <- min(4L,length(unique(covariate)))
    #		If covariate takes only one value, recall with NULL covariate
    if(splinedf < 2L) {
      out <- Recall(x=x,df1=df1)
      out$scale <- rep_len(out$scale,n)
      return(out)
    }
  }

  #	Remove missing or infinite values and zero degrees of freedom
  ok <- ok & is.finite(x) & (x > -1e-15)
  nok <- sum(ok)
  notallok <- !all(ok)
  if(notallok) {
    x <- x[ok]
    if(length(df1)>1L) df1 <- df1[ok]
    if(!is.null(covariate)) {
      covariate.notok <- covariate[!ok]
      covariate <- covariate[ok]
    }
  }

  #	Check whether enough observations to estimate variance around trend
  if(nok <= splinedf) {
    s20 <- NA
    if(!is.null(covariate)) s20 <- rep_len(s20,n)
    return(list(scale=s20,df2=NA))
  }

  #	Avoid exactly zero values
  x <- pmax(x,0)
  m <- median(x)
  if(m==0) {
    warning("More than half of residual variances are exactly zero: eBayes unreliable")
    m <- 1
  } else {
    if(any(x==0)) warning("Zero sample variances detected, have been offset away from zero",call.=FALSE)
  }
  x <- pmax(x, 1e-5 * m)

  #	Better to work on with log(F)
  z <- log(x)
  e <- z-digamma(df1/2)+log(df1/2)

  if(is.null(covariate)) {
    emean <- mean(e)
    evar <- sum((e-emean)^2)/(nok-1L)
  } else {
    if(!requireNamespace("splines",quietly=TRUE)) stop("splines package required but is not available")
    design <- try(splines::ns(covariate,df=splinedf,intercept=TRUE),silent=TRUE)
    if(is(design,"try-error")) stop("Problem with covariate")
    fit <- lm.fit(design,e)
    if(notallok) {
      design2 <- predict(design,newx=covariate.notok)
      emean <- rep_len(0,n)
      emean[ok] <- fit$fitted
      emean[!ok] <- design2 %*% fit$coefficients
    } else {
      emean <- fit$fitted
    }
    evar <- mean(fit$effects[-(1:fit$rank)]^2)
  }

  #	Estimate scale and df2
  evar <- evar - mean(trigamma(df1/2))
  if(evar > 0) {
    df2 <- 2*trigammaInverse(evar)
    s20 <- exp(emean+digamma(df2/2)-log(df2/2))
  } else {
    df2 <- Inf
    if(is.null(covariate))
      #			Use simple pooled variance, which is MLE of the scale in this case.
      #			Versions of limma before Jan 2017 returned the limiting value of the evar>0 estimate, which is larger.
      s20 <- mean(x)
    else
      s20 <- exp(emean)
  }

  list(scale=s20,df2=df2)
}

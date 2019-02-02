##Once obtaining the s0^2 and n0, for each trio, we create the moderated t-statistics (a_2g_hat)/((s_g^2_tilde)*sqrt(v_g)) where s_g^2_tilde=(n0*s0^2+df_g*s_g^2)/(n0+df_g)
##This function gives the estimate of v0. The imput tstat is the vector of moderated t-statistics for all trios in one tissue
##The input stdev.unscaled is the vector of v_g for all trios in one tissue, taking sqrt ##lin
##The input df shoud be the vector of the degree of freedom for the moderated t-statistics (for all trios in one tissue) which is df_g+n0.
tmixture.vector <- function(tstat,stdev.unscaled,df,proportion=0.01,v0.lim=NULL)
  #	Estimate scale factor in mixture of two t-distributions
  #	tstat is assumed to follow sqrt(1+v0/v1)*t(df) with probability proportion and t(df) otherwise
  #	v1 is stdev.unscaled^2 and v0 is to be estimated
  #	Gordon Smyth
  #	18 Nov 2002.  Last modified 15 April 2016.
{
  #	Remove missing values
  if(anyNA(tstat)) {
    o <- !is.na(tstat)
    tstat <- tstat[o]
    stdev.unscaled <- stdev.unscaled[o]
    df <- df[o]
  }

  #	ntarget t-statistics will be used for estimation
  ngenes <- length(tstat)
  ntarget <- ceiling(proportion/2*ngenes)
  if(ntarget < 1) return(NA)

  #	If ntarget is v small, ensure p at least matches selected proportion
  #	This ensures ptarget < 1
  p <- max(ntarget/ngenes,proportion)

  #	Method requires that df be equal
  tstat <- abs(tstat)
  MaxDF <- max(df)
  i <- df < MaxDF
  if(any(i)) {
    TailP <- pt(tstat[i],df=df[i],lower.tail=FALSE,log.p=TRUE)
    tstat[i] <- qt(TailP,df=MaxDF,lower.tail=FALSE,log.p=TRUE)
    df[i] <- MaxDF
  }

  #	Select top statistics
  o <- order(tstat,decreasing=TRUE)[1:ntarget]
  tstat <- tstat[o]
  v1 <- stdev.unscaled[o]^2

  #	Compare to order statistics
  r <- 1:ntarget
  p0 <- 2*pt(tstat,df=MaxDF,lower.tail=FALSE)
  ptarget <- ( (r-0.5)/ngenes - (1-p)*p0 ) / p
  v0 <- rep.int(0,ntarget)
  pos <- ptarget > p0
  if(any(pos)) {
    qtarget <- qt(ptarget[pos]/2,df=MaxDF,lower.tail=FALSE)
    v0[pos] <- v1[pos]*((tstat[pos]/qtarget)^2-1)
  }
  if(!is.null(v0.lim)) v0 <- pmin(pmax(v0,v0.lim[1]),v0.lim[2])
  mean(v0)
}

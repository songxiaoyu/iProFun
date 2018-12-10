####For each trio (trans, cis, snp)and the corresponding covariates, we run the following linear model for the gth trio:
## trans =a_0g+a_1g*snp+a_2g*cis+a_3g*covariates+epsilon_g, epsilon~N(0,sigma_g^2)
## Then a_2g_hat | a_2g, sigma_g^2 ~ N(a_2, sigma_g^2*v_g) where v_g is the (3,3) entry in the matrix (X^T X)^{-1}, X is the design matrix (1,snp, cis, covariates)
## mean squared error s_g^2 | sigma_g^2 ~ [sigma_g^2/(n_g-m_g-2)]Chi^2_{df=n_g-m_g-2} where n_g is sample size, m_g is the number of confounders in covariates matrix.
## priors: 1/sigma_g^2 ~ [1/(n0+s0^2)]Chi^2_{df=n0}, a_2g|a_2g!=0, sigma_g^2 ~ N(0,v0*sigma_g^2). n0, s0^2 and v0 are hyperparameters to be estimated.


# install.packages("metRology")
# install.packages("matrixStats")
require(metRology)
require(matrixStats)

##estimates is a function that gives the needed output from each trio's regression
my.solve <- function(X) {
  if (!is.matrix(X))
    X <- matrix(X, nrow = sqrt(length(X)))
  ss <- svd(X)
  Xinv <- ss$u %*% diag(1/ss$d, nrow = nrow(X), ncol = nrow(X)) %*% t(ss$v)
  return(Xinv)
}
estimates<-function(trans,snp,cis,covariates){
	n_g<-length(snp)
	x<-cbind(1,snp,cis,covariates)
	p_x<-dim(x)[2]
	inverse_xx<-my.solve(t(x) %*% x)
	alpha<-inverse_xx %*% t(x) %*% trans
	s_g_square<-as.numeric(1/(n_g - p_x) * (sum(trans^2) - t(trans) %*% x %*% inverse_xx %*% t(x) %*% trans)) 
  	return(list(a_2g_hat = alpha[3], s_g_square = s_g_square, v_g =inverse_xx[3,3], df_g=n_g-p_x))
##a_2g_hat is defined before, s_g_square is s_g^2, v_g is defined before, df_g is essentially n_g-m_g-2
}


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



Estep<-function(cl, oldpi, n_trio, n_pattern, q, density_0, density_1){
  Bmatrix<-matrix(NA,nrow=n_trio,ncol=n_pattern)
  
  if(!is.null(cl)){
    Bmatrix <- parSapply(cl, 1:n_pattern,function(j,q,density_0,density_1,oldpi){
      nT=ncol(q)
      m=nrow(density_0)
      Iq <- matrix(rep(q[j,],m),byrow=T,ncol=nT)
      Pi <- rep(log(oldpi[j]),m)
      oo <- Pi+rowSums((1-Iq)*log(density_0)+Iq*log(density_1))
      return(oo)
    },q=q, density_0=density_0, density_1=density_1, oldpi=oldpi)
  }else{
    nT=ncol(q)
    m=nrow(density_0)
    for(j in 1:n_pattern){
      Iq <- matrix(rep(q[j,],m),byrow=T,ncol=nT)
      Pi <- rep(log(oldpi[j]),m)
      Bmatrix[,j] <- Pi+rowSums((1-Iq)*log(density_0)+Iq*log(density_1))
    }
  }
  
  
  ##' consider matrixStats::rowMins() or rpgm::rowMins()
  Bmatrix<-Bmatrix-rowMins(Bmatrix)
  Bmatrix<-exp(Bmatrix)	
  Bmatrix<-Bmatrix/rowSums(Bmatrix)
  Lst<-list(Bmatrix=Bmatrix)
  Lst
}



##Mstep obtains the pi that maximizes the posterior expectation function in the expecation step.
Mstep<-function(oldB, n_trio, n_pattern){
#  newpi<-(apply(oldB,2,sum)+1)/(n_trio+n_pattern)
  newpi<-colSums(oldB)/(n_trio)
  Lst<-list(newpi=newpi)
  Lst
}


makeQ <- function(grp){
  nT <- length(grp)
  ngi=-sort(-table(grp))
  Q <- rep(0,nT)
  ng <- max(grp)
  ngi <- rep(1,ng)
  
  for (qi in 1:ng){
    oo <- combn(1:ng,qi)
    Q <- cbind(Q, apply(oo,2,function(x,ng,ngi) {
      vv <- rep(0,ng)
      vv[x] <- 1
      oo <- rep(vv,ngi)
      return(oo)
    }, ng=ng, ngi=ngi))
  }
  Q <- t(Q)
  return(Q)
}


## this is the final function to co-map QTLs.
## betas are a matrix of beta-coefficient each row is a testing feature/gene/locus, each column is a study/data-source/omics-data-type
## sds are the standard deviation corresponding to the betas
## for SNP level test, MAFs are the minor allele frequencies
## df is the degree of freedom for each study, sample size minus number of covariates minus 1
## pi1 is the specified proportion of statistics likely from alternative
## cl is using parallel computing, if the number of test is not large, you can ignore it
## tol is the covergence criteria


mQTLoc_func <- function(betas, sds, mafs, dfs, pi1, cl=NULL, tol=1e-3){
  
  m <- nrow(betas)
  d <- ncol(betas)
  vg = 1/(2*mafs*(1-mafs))
  sigma2 <- sds^2*matrix(rep(2*mafs*(1-mafs), each=d), byrow=T, ncol=d)
  
  Tstat_m <- D0 <- D1 <- NULL
  
  for (j in 1:d){
    d1=rep(dfs[j],m)
    xx <- fitFDist(sigma2[,j],d1)
    s02 <- xx$scale; n0 <- xx$df2
    sg_tilde <- sqrt((n0*s02+d1*sigma2[,j])/(n0+d1))
    moderate.t <- betas[,j]/(sg_tilde*sqrt(vg))
    Tstat_m <- cbind(Tstat_m, moderate.t)
    
    #pt <- 2*(1-pt(abs(moderate.t), df=d1+n0))
    #pi0 <- 1
    #try({pi0 <- qvalue(pt)$pi0}, silent=TRUE)
    
    D0 <- cbind(D0, dt(moderate.t, df=d1+n0))
    v0 <- tmixture.vector(moderate.t, sqrt(vg),d1+n0,proportion=pi1[j],v0.lim=NULL)  
    scaler=sqrt(1+v0/vg)
    D1 <- cbind(D1, dt.scaled(moderate.t,df=d1+n0,mean=0,sd=scaler))
  }
  
  
  Q <- makeQ(1:d)
  n_pattern <- nrow(Q)
  curpi<- c(0.80, rep((1-0.80)/(n_pattern-1),n_pattern-1))
  diff<-1
  numiters<-1
  itermat<-curpi
  while(diff>tol){
    numiters<-numiters+1
    curestep<-Estep(cl, curpi, n_trio=m, n_pattern=n_pattern, q=Q, density_0=D0, density_1=D1)
    curb<-curestep[[1]]
    curmstep<-Mstep(curb, n_trio=m, n_pattern=n_pattern)
    curpi<-curmstep[[1]]
    itermat<-rbind(itermat,curpi)
    diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
    # print(c(numiters,diff, curpi))
  }
  
  
  
  return(list(PostProb = curb, colocProb = curpi, Tstat_m = Tstat_m, D0=D0, D1=D1))
}


get.fdr <- function(pp){
  fdr <- 1-cumsum(pp[order(pp,decreasing=T)])/(1:length(pp))
  idx <- match(1:length(pp), order(pp,decreasing=T))
  fdr <- fdr[idx]
  return(fdr)
}


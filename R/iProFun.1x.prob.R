# a single x data type
iProFun.1x.prob <- function(input, pi1,
                             miss.platform.include=T, cl=NULL, tol=1e-3){
  betas_J=input$betas_J
  dfs_J=input$dfs_J
  betas_se_J=input$betas_se_J
  sigma2_J=input$sigma2_J
  v_g_J=input$v_g_J
  xName=input$xName_J
  yName=input$yName_J

  ylength=ncol(betas_J)

  # use genes with no missing data types to calculate the parameters
  Tstat_L <- D0 <- D1 <- NULL
  s02_J=n0_J=v0_J=NULL
  for (j in 1:ylength){
    d1=dfs_J[,j]
    density <- fitFDist(sigma2_J[,j],d1 )# obtain d_0 and s02
    s02 <- density$scale;
    n0 <- density$df2
    if (n0==Inf) {n0=1000000}
    sg_tilde <- sqrt((n0*s02+d1*sigma2_J[,j])/(n0+d1))
    moderate.t <- betas_J[,j]/(sg_tilde*sqrt(v_g_J[,j]))
    Tstat_L <- cbind(Tstat_L, moderate.t)

    D0 <- cbind(D0, dt(moderate.t, df=d1+n0))
    v0 <- tmixture.vector(moderate.t, sqrt(v_g_J[,j]),d1+n0,proportion=pi1[j],v0.lim=NULL)
    scaler=sqrt(1+v0/v_g_J[,j])
    D1 <- cbind(D1, dt.scaled(moderate.t,df=d1+n0,mean=0,sd=scaler))
    s02_J=c(s02_J, s02)
    n0_J=c(n0_J, n0)
    v0_J=c(v0_J, v0)
  }

  # estimate parameters using all non-missing data types
  Q <- makeQ(1:ylength)
  colnames(Q)= 1:ylength
  n_pattern <- nrow(Q)
  curpi<- c(0.80, rep((1-0.80)/(n_pattern-1),n_pattern-1))
  diff<-1
  numiters<-1
  itermat<-curpi


  while(diff>tol & numiters<200 ){
    #print(diff)
    numiters<-numiters+1
    curestep<-Estep(cl, oldpi=curpi, n_trio=nrow(Tstat_L), n_pattern=n_pattern, q=Q, density_0=D0, density_1=D1)
    curb<-curestep[[1]]
    curmstep<-Mstep(curb, n_trio=nrow(Tstat_L), n_pattern=n_pattern)
    curpi<-curmstep[[1]]
    itermat<-rbind(itermat,curpi)
    diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])

  }

  # adding genes with missing outcome platforms
  MissIndex=apply(betas_J, 1, function(f) which(is.na(f)))
  MissPattern=unique(MissIndex)
  MissPattern=MissPattern[(sapply(MissPattern, length)!=0 & sapply(MissPattern, length)!=ylength)]




  if (miss.platform.include==F | length(MissPattern)==0) {
    return(list(NoComputation=dim(betas_J)[1],
                Config=Q, Config.miss=NULL,
                PostProb = curb, PostProb.miss= NULL,
                xName.miss=NULL,
                colocProb = curpi, Tstat_L = Tstat_L, D0=D0, D1=D1))
  } else {

    # use above estimated parameters to estimate postprob for genes with missing in some data types

    D0t=D1t=Tstatt=Q1=n_pattern1=colocProb1=xName.miss=vector("list",  length(MissPattern))

    for ( k in 1:length(MissPattern)) {
      trait_idx=setdiff(1:ylength, MissPattern[[k]])
      t=which(sapply(MissIndex, function(f) identical(f, MissPattern[[k]])))
      xName.miss[[k]]=xName[t,]

      for (j in trait_idx ) {

        d1=dfs_J[t,j]
        sg_tilde <- sqrt((n0_J[j]*s02_J[j]+d1*sigma2_J[t,j])/(n0_J[j]+d1))
        moderate.t <- betas_J[t,j]/(sg_tilde*sqrt(v_g_J[t,j]))
        Tstatt[[k]] <- cbind(Tstatt[[k]], moderate.t)

        D0t[[k]] <- cbind(D0t[[k]], dt(moderate.t, df=d1+n0_J[j]))
        scaler=sqrt(1+v0_J[j]/v_g_J[t,j])
        D1t[[k]] <- cbind(D1t[[k]], dt.scaled(moderate.t,df=d1+n0_J[j],mean=0,sd=scaler))

      }

      Q1[[k]] <- makeQ(1:length(trait_idx))
      colnames(Q1[[k]])= trait_idx
      n_pattern1[[k]] <- nrow(Q1[[k]])
      pattern_index=row.match(data.frame(Q[,trait_idx]), data.frame(Q1[[k]]))

      curpi1<- tapply(curpi, pattern_index, sum)
      curestep1<-Estep(cl, oldpi=curpi1, n_trio=nrow(Tstatt[[k]]), n_pattern=n_pattern1[[k]], q=Q1[[k]], density_0=as.matrix(D0t[[k]]), density_1=as.matrix(D1t[[k]]))
      colocProb1[[k]]=curestep1[[1]]
    }

    return(list(NoComputation=dim(betas_J)[1],
                Config=Q, Config.miss=Q1,
                PostProb = curb,PostProb.miss= colocProb1,
                xName.miss=xName.miss,
                colocProb = curpi, Tstat_L = Tstat_L, D0=D0, D1=D1))

  }# end missing platfor codes



}

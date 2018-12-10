MultiOmics_Input <- function(input, pi1=c(0.0001, 0.0001, 0.0001), cl=NULL, tol=1e-3){
  betas_J=input$betas_J
  dfs_J=input$dfs_J
  betas_se_J=input$betas_se_J
  sigma2_J=input$sigma2_J
  v_g_J=input$v_g_J
  xName=input$xName
  yName=input$yName
  ylength=ncol(betas_J)
  platform = input$name

  Tstat_L <- D0 <- D1 <- NULL

  for (j in 1:ylength){
    d1=dfs_J[,j]
    density <- fitFDist(sigma2_J[,j],d1 )# obtain d_0 and s02
    s02 <- density$scale; n0 <- density$df2
    if (n0==Inf) {n0=100000}
    sg_tilde <- sqrt((n0*s02+d1*sigma2_J[,j])/(n0+d1))
    moderate.t <- betas_J[,j]/(sg_tilde*sqrt(v_g_J[,j]))
    Tstat_L <- cbind(Tstat_L, moderate.t)

    D0 <- cbind(D0, dt(moderate.t, df=d1+n0))
    v0 <- tmixture.vector(moderate.t, sqrt(v_g_J[,j]),d1+n0,proportion=pi1[j],v0.lim=NULL)
    scaler=sqrt(1+v0/v_g_J[,j])
    D1 <- cbind(D1, dt.scaled(moderate.t,df=d1+n0,mean=0,sd=scaler))
  }


  Q <- makeQ(1:ylength)
  n_pattern <- nrow(Q)
  curpi<- c(0.80, rep((1-0.80)/(n_pattern-1),n_pattern-1))
  diff<-1
  numiters<-1
  itermat<-curpi
  while(diff>tol){
    numiters<-numiters+1
    curestep<-Estep(cl, curpi, n_trio=nrow(Tstat_L), n_pattern=n_pattern, q=Q, density_0=D0, density_1=D1)
    curb<-curestep[[1]]
    curmstep<-Mstep(curb, n_trio=nrow(Tstat_L), n_pattern=n_pattern)
    curpi<-curmstep[[1]]
    itermat<-rbind(itermat,curpi)
    diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
    print(c(numiters,diff, curpi))
  }

  betas_J=cbind(xName, yName, betas_J)
  betas_se_J=cbind(xName, yName, betas_se_J)
  curb=cbind( xName, curb) ### with name here

  return(list(NoComputation=dim(betas_J)[1], Config=Q, PostProb = curb, colocProb = curpi, Tstat_L = Tstat_L, D0=D0, D1=D1))
}

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

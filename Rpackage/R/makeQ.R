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

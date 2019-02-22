inverse <- function(X) {
  if (!is.matrix(X))
    X <- matrix(X, nrow = sqrt(length(X)))
  ss <- svd(X)
  Xinv <- ss$u %*% diag(1/ss$d, nrow = nrow(X), ncol = nrow(X)) %*% t(ss$v)
  return(Xinv)
}

## Inputs ##
# N: number of faces to be simulated
# delta: size of signal of X~N(0,1)
# err.sd: error standard deviation

## Outputs ##
# Y: N x 3 x 7150 (7150 points of x,y,z-coordinate for N faces)
# X: N x 1

simulate_face <- function(N, delta = 100, err.sd = 0.002) {

  require(MASS)
  require(tensorA)

  load("./rundata/for_simulation.Rdata")
  # includes Ymean, Yvar, eval_eig, pt_mean, beta

  P <- 7150

  pred <- mvrnorm(N, Ymean, Yvar)

  #nparam = nrow(pred) ## would be N
  pred <- to.tensor(c(pred), c(nparam = N, lambda=ncol(pred)))
  ptwise <- mul.tensor(pred, "lambda", eval_eig, "lambda")

  xmat <- matrix(unlist(t(ptwise[,1,])), nrow = P, ncol = N) + 
		matrix(rep(unlist(pt_mean$x),N), nrow=P, ncol = N)
  ymat <- matrix(unlist(t(ptwise[,2,])), nrow = P, ncol = N) + 
		matrix(rep(unlist(pt_mean$y),N), nrow=P, ncol = N)
  zmat <- matrix(unlist(t(ptwise[,3,])), nrow = P, ncol = N) + 
		matrix(rep(unlist(pt_mean$z),N), nrow = P, ncol = N)

  Y <- matrix(0, nrow=(3*P), ncol=N)
  Y[(3*c(1:7150)-2),] <- xmat
  Y[(3*c(1:7150)-1),] <- ymat
  Y[(3*c(1:7150)-0),] <- zmat

  # simulation
  X <- rnorm(N, 0, 1)
  err <- matrix(rnorm(N*P*3, 0, err.sd), nrow=(3*P), ncol=N)
  tmp <- Y + delta * (as.matrix(beta) %*% t(as.matrix(X))) + err

  Y <- array(0, dim=c(N, 3, P))
  Y[,1,] <- t(tmp[(3*c(1:7150)-2),])
  Y[,2,] <- t(tmp[(3*c(1:7150)-1),])
  Y[,3,] <- t(tmp[(3*c(1:7150)-0),])

  return(list(Y=Y, X=X))
}
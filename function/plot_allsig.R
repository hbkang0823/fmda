## Inputs ##
# fitted_beta: the beta to be plotted
# sig.level: the significance level
# Ymat: coefficient matrix of basis functions for functional Y's (N x K matrix)
# Xmat: matrix of scalar predictors (N x R matrix)
# SVD: SVD result from 'fpca2' function
# ref: reference face for the plot to be on (P x 3 matrix)
# evalmat: P x nharms1 matrix
# orthVec: orthogonal vectors relative to reference face (P x 3 matrix)

## Outputs ##
# plot that gives the overall significance of the given beta

plot_allsig <- function(fitted_beta, sig.level = .95, Ymat, Xmat, SVD, ref, evalmat, orthVec=c()) {

  require(tensorA)

  if (length(orthVec) == 0) {
    orthVec <- get_orth_vec(ref)
  }
 
  eigmat <- SVD$v[,,1:nharms2]
  evalmat <- to.tensor(c(evalmat), c(P = dim(evalmat)[1], J=dim(evalmat)[2]))
  eval_eig <- mul.tensor(eigmat, "J", evalmat, "J")
  if (length(dim(fitted_beta)) == 0) {
    fitted_beta.t <- to.tensor(c(fitted_beta), c(ncoef = 1, lambda = nharms2))
    Ni = 1; i=1
  } else {
    fitted_beta.t <- to.tensor(c(fitted_beta), c(ncoef = dim(fitted_beta)[1], lambda = dim(fitted_beta)[2]))
  }
  eval_beta <- mul.tensor(fitted_beta.t, "lambda", eval_eig, "lambda")
  eval_beta <- t(matrix(unlist(eval_beta), nrow=dim(eval_beta)[2], ncol=dim(eval_beta)[3]))

  diff <- (Ymat - Xmat %*% t(matrix(fitted_beta)))
  diff.t <- to.tensor(c(diff), c(N = nrow(diff), lambda = ncol(diff)))
  eval_diff <- mul.tensor(diff.t, "lambda", eval_eig, "lambda")


  require(expm)

  N <- dim(eval_diff)[1]
  P <- dim(eval_diff)[3]
  half_corr_pt <- eval_diff

  t1 = proc.time()
  for (p in 1:P) {
    half_corr_pt[,,p] <- eval_diff[,,p] %*% sqrtm(solve( cov(eval_diff[,,p]) ))
  }
  proc.time()-t1

  eSVD <- svd.tensor(half_corr_pt/sqrt(N-1), i="N", j=c("I","P"))  


  ######
  invXX <- solve(t(Xmat)%*%(Xmat))


  betacoef_rotated <- eval_beta
   for (p in 1:P) {
     tmp <- matrix(c(t(eval_diff[,,p]) %*% eval_diff[,,p] / (N-Ni)), nrow=3)
     betacoef_rotated[p,] <- c(sqrtm(solve(N*invXX[i,i]*tmp)) %*% matrix(eval_beta[p,]))
  }

  beta_norm_sq <- matrix(NA, nrow=P, ncol=dim(Xmat)[2])
  for (i in 1:ncol(beta_norm_sq)) {
    for (p in 1:P) {
      beta_norm_sq[p,i] <- t(matrix(betacoef_rotated[p,],nrow=3))%*%matrix(betacoef_rotated[p,],nrow=3)
    }
  }

  inprod_vec_beta <- diag(eval_beta %*% t(orthVec))


  s = sum((eSVD$d)^2)*2 / sum((eSVD$d))	# scale parameter
  a = sum((eSVD$d))/s				# shape parameter
  critval <- qgamma(sig.level, shape=a, scale=s)
  #require('CompQuadForm'); print(imhof(critval, (eSVD$d)))

  r_ball_sq <- integer(P)
  for (p in 1:P) {
    r_ball_sq[p] <- critval * sum(  diag( eSVD$v[,p,] %*% diag( eSVD$d / sqrt(N) ) %*%
  						               t(eSVD$v[,p,]) )   )
  }

  norm_outside_ball <- ((N*beta_norm_sq) > r_ball_sq)

  ####################################################
  # plot

  conf_dir <- sign(inprod_vec_beta) * (norm_outside_ball)

  qbreaks = c(-1, -0.5, 0.5, 1)
  color = c("blue", "grey", "red")

  rglplot(conf_dir, ref, qbreaks, color)

}
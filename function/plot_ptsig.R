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
# plot that gives the pointwise significance of the given beta

plot_ptsig <- function(fitted_beta, sig.level=.95, Ymat, Xmat, SVD, ref, evalmat, orthVec=c()) {

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


  P <- dim(eval_diff)[3]
  diag_e_cov_pt <- to.tensor(rep(0, (3*3*7150)), c(I1 = 3, I2 = 3, P = 7150))
  for (p in 1:P) {
    diag_e_cov_pt[,,p] <- t(eval_diff[,,p]) %*% eval_diff[,,p] / (N-Ni)
  }
  for (p in 1:P) {
    diag_e_cov_pt[,,p] <- cov(eval_diff[,,p])
  }



  invXX <- solve(t(Xmat)%*%(Xmat))

  require(expm)

  betacoef_rotated <- eval_beta
   for (p in 1:P) {
     tmp <- matrix(c(diag_e_cov_pt[,,p]), nrow=dim(diag_e_cov_pt)[1])
     betacoef_rotated[p,] <- c(sqrtm(solve(N*invXX[i,i]*tmp)) %*% matrix(eval_beta[p,]))
  }

  beta_norm_sq <- matrix(NA, nrow=P, ncol=dim(Xmat)[2])
  for (i in 1:ncol(beta_norm_sq)) {
    for (p in 1:P) {
      beta_norm_sq[p,i] <- t(matrix(betacoef_rotated[p,],nrow=3))%*%matrix(betacoef_rotated[p,],nrow=3)
    }
  }

  inprod_vec_beta <- diag(eval_beta %*% t(orthVec))


  ####################################################
  # plot
  Pk_mat <- matrix(NA, nrow=P, ncol=Ni) 
  for (k in 1:Ni) { 				
    Pk_mat[,k] <- 1 - pchisq((N*beta_norm_sq[,k]), df=3)	#p-values
  }

  sig_dir <- sign(inprod_vec_beta) * (1-Pk_mat)

  qbreaks = c(-1, -sig.level, sig.level, 1)
  color = c("blue", "grey", "red")

  rglplot(sig_dir, ref, qbreaks, color)

}
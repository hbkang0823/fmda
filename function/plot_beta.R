## Inputs ##
# beta: the beta to be plotted. K x 1.
# SVD: SVD result from 'fpca2' function
# nharms2: same as K
# ref: reference face for the plot to be on (P x 3 matrix)
# evalmat: P x nharms1 matrix
# orthVec: orthogonal vectors relative to reference face (P x 3 matrix)

## Outputs ##
# plot

plot_beta <- function(fitted_beta, SVD, nharms2, ref, evalmat, orthVec=c()) {

  if (length(orthVec) == 0) {
    orthVec <- get_orth_vec(obj)
  }
 
  eigmat <- SVD$v[,,1:nharms2]
  evalmat <- to.tensor(c(evalmat), c(P = dim(evalmat)[1], J=dim(evalmat)[2]))
  eval_eig <- mul.tensor(eigmat, "J", evalmat, "J")
  fitted_beta.t <- to.tensor(c(fitted_beta), c(ncoef = 1, lambda = nharms2))
  eval_beta <- mul.tensor(fitted_beta.t, "lambda", eval_eig, "lambda")
  eval_beta <- t(matrix(unlist(eval_beta), nrow=dim(eval_beta)[2], ncol=dim(eval_beta)[3]))

  # get the inner product between the orthogonal vector and beta.
  inprod_vec_beta <- diag(eval_beta %*% t(orthVec))

  # plot
  rglplot(inprod_vec_beta, ref, qbreaks=c(-1,0,1), color=c("blue", "red"))

}
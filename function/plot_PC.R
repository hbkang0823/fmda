## Inputs ##
# pc: the principal component to be plotted. nharms1 x 3.
# ref: reference face for the plot to be on (P x 3 matrix)
# evalmat: P x nharms1 matrix
# orthVec: orthogonal vectors relative to reference face (P x 3 matrix)

## Outputs ##
# plot

plot_PC <- function(pc, ref, evalmat, orthVec=c()) {

  if (length(orthVec) == 0) {
    orthVec <- get_orth_vec(ref)
  }
 

  # get the inner product between the orthogonal vector and PC.
  pc <- t(matrix(unlist(pc), nrow=dim(pc)[1], ncol=dim(pc)[2]))
  evalmat <- matrix(unlist(evalmat),nrow=dim(evalmat)[1])
  eval_pc <- evalmat%*%pc ## PC * evalmat

  inprod_vec_pc <- diag(eval_pc %*% t(orthVec))


  # plot
  require(colorRamps)
  limcol = ceiling(max(abs(max(inprod_vec_pc)), abs(min(inprod_vec_pc))))
  limseq = seq(0, limcol, by = 4)
  qbreaks = c(-limseq, 0, limseq); qbreaks = unique(qbreaks); qbreaks = qbreaks[order(qbreaks)]
  color = blue2red(length(qbreaks))

  rglplot(inprod_vec_pc, ref, qbreaks, color)

}
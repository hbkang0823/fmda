## Inputs ##
# result: the result from fpca1.m
# SVD: SVD result from 'fpca2' function
# nharms2: number of 

## Outputs ##
# Laplacian penalty matrix

get_penalty <- function(result, SVD, nharms2) {
  require(tensorA)

  nharms1 <- dim(result$evalmat)[2]

  PsiLap <- result$PsiLap
  PsiLap <- to.tensor(c(PsiLap), c(Ji = nharms1, Jj=nharms1))
  
  tmp1 <- mul.tensor(SVD$v[1,,1:nharms2], "I1", PsiLap, "Ji")
  names(tmp1) <- c("K1", "Jj")
  tmp1 <- mul.tensor(tmp1, "Jj", SVD$v[1,,1:nharms2],"I1")
  names(tmp1) <- c("K1", "K2")

  tmp2 <- mul.tensor(SVD$v[2,,1:nharms2], "I1", PsiLap, "Ji")
  names(tmp2) <- c("K1", "Jj")
  tmp2 <- mul.tensor(tmp2, "Jj", SVD$v[2,,1:nharms2],"I1")
  names(tmp2) <- c("K1", "K2")

  tmp3 <- mul.tensor(SVD$v[3,,1:nharms2], "I1", PsiLap, "Ji")
  names(tmp3) <- c("K1", "Jj")
  tmp3 <- mul.tensor(tmp3, "Jj", SVD$v[3,,1:nharms2],"I1")
  names(tmp3) <- c("K1", "K2")

  Vpenmat <- tmp1 + tmp2 + tmp3

  return(Vpenmat)
}
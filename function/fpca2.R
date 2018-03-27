fpca2 <- function(coef) {

  require(tensorA)

  nface <- dim(coef)[1]/3
  nbasis <- dim(coef)[2]

  coef <- to.tensor(c(coef), c(N=nface, I = 3, J = nbasis))
  SVD <- svd.tensor(coef/sqrt(N-1), i="N", j=c("I","J"))

#  names(coef) <- c("N", "I1", "J1")
#  coef2 <- coef
#  names(coef2) <- c("N", "I2", "J2")

#  Sigmat <- mul.tensor(coef, "N", coef2, "N")/(N-1)

#  SVD <- svd.tensor(Sigmat, c("J1", "I1"), c("J2","I2"))
  
  return(SVD)
}
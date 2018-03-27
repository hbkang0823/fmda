## Inputs ##
# Ymat: coefficient matrix of basis functions for functional Y's (N x K matrix)
# Xmat: matrix of scalar predictors (N x R matrix)
# penmat: penalty matrix (K x K matrix)
#		if no penalty, use 
#			penmat = matrix(0,nrow=dim(Ymat)[2],ncol=dim(Ymat)[2])
# penlamb_all: possible values for penalty lambda 

## Outputs ##
# Bhat: coefficient matrix for functional beta that takes same basis as Y's
#		(R x K matrix)
# pvals: test statistics size and p-values based on L2-test, PC-test, 
#		and Choi-test


MS_penreg <- function(Ymat, Xmat, penmat, penlamb_all = c(0,10^(-5:10))) {
  
  	N = dim(Ymat)[1]
  	K = dim(Ymat)[2]
  	R = dim(Xmat)[2]

  	if (dim(Xmat)[1] != N) {
    	stop("Check the dimension of Ymat and Xmat")
  	}
  	if (dim(penmat)[1] != K | dim(penmat)[2] != K) {
    	stop("Check the dimension of penmat")
  	}
  	
  	penlamb<-rep(0,times=R)
	if (R == 1) {
		Lamb = penlamb
		A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb * penmat
	} else {
	  	Lamb <- diag(penlamb)
	  	A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
	}

  	# Penalized Estimation
  	V = c(t(Ymat)%*%Xmat)
  	A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
  	Ainv = chol2inv(chol(A))
  	Bhat_v = Ainv%*%V
  	Bhat_m_T = matrix(Bhat_v,nrow=K,ncol=R)
  	Bhat = t(Bhat_m_T)	### Penalized LSE
  	Eps_mat = Ymat - Xmat%*%Bhat
  	df<-N-sum(diag(Ainv%*%((t(Xmat)%*%Xmat)%x%diag(1,nrow=K))))/K
  	GCV<-N*sum((Eps_mat^2))/df^2
  	
  	max_iter<-10;cnt<-0
	lamb_ind<-rep(1,times=R)
	go<-TRUE
  	while(go==TRUE & cnt<max_iter){
  		go<-FALSE
  		for(j in 1:R){
  			if(lamb_ind[j]<length(penlamb_all)){
	  			penlamb_new<-penlamb
	  			penlamb_new[j]<-penlamb_all[lamb_ind[j]+1]
				if (R == 1) {
					Lamb = penlamb
					A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb * penmat
				} else {
				  	Lamb <- diag(penlamb)
				  	A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
				}
  				# Penalized Estimation
  				V = c(t(Ymat)%*%Xmat)
  				A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
  				Ainv = chol2inv(chol(A))
  				Bhat_v = Ainv%*%V
  				Bhat_m_T = matrix(Bhat_v,nrow=K,ncol=R)
  				Bhat = t(Bhat_m_T)	### Penalized LSE
  				Eps_mat = Ymat - Xmat%*%Bhat
  				df<-N-sum(diag(Ainv%*%((t(Xmat)%*%Xmat)%x%diag(1,nrow=K))))/K
  				GCV_new<-N*sum((Eps_mat^2))/df^2
  				if(GCV_new<=GCV){
  					GCV<-GCV_new
  					penlamb<-penlamb_new
  					lamb_ind[j]<-lamb_ind[j]+1
  					go<-TRUE
  				}
			}
  		}
  		cnt<-cnt+1
  		
  	}
  	
  	
	if (R == 1) {
		Lamb = penlamb
		A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb * penmat
	} else {
	  	Lamb <- diag(penlamb)
	  	A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
	}
  	# Penalized Estimation
  	V = c(t(Ymat)%*%Xmat)
  	A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
  	Ainv = chol2inv(chol(A))
  	Bhat_v = Ainv%*%V
  	Bhat_m_T = matrix(Bhat_v,nrow=K,ncol=R)
  	Bhat = t(Bhat_m_T)	### Penalized LSE
  	Eps_mat = Ymat - Xmat%*%Bhat
  	Sig_mat = (1/(N - R)) * t(Eps_mat)%*%Eps_mat
  	Big_Cov_mat = Ainv%*%((t(Xmat)%*%Xmat)%x%Sig_mat)%*%Ainv


  	# Testing
  	EV_PC<-.85
  	EV_Choi<-.99
  	EV_L2<-.999

  	pvalmat <- data.frame(betanum=1:nrow(Bhat), Ti_L2 = NA, pval_L2 = NA, 
						Ti_PC = NA, pval_PC = NA,
						Ti_Choi = NA, pval_Choi = NA)
  	require('CompQuadForm')
  	for (i in 1:nrow(Bhat)) {

    	covBi = Big_Cov_mat[(K*(i-1)+1):(K*i),(K*(i-1)+1):(K*i)]
    	Eig = eigen(covBi, symmetric=T)
    	eigvals_i <- Eig$values
    	eigvecs_i <- Eig$vectors

    	Scores<-Bhat%*%eigvecs_i

    	k_PC<-1+sum(cumsum(eigvals_i)/sum(eigvals_i)<EV_PC)
    	k_Choi<-1+sum(cumsum(eigvals_i)/sum(eigvals_i)<EV_Choi)
    	k_L2<-1+sum(cumsum(eigvals_i)/sum(eigvals_i)<EV_L2)
    	 
    	Ti_L2 = sum( Scores[i,1:k_L2]^2 )
    	Ti_PC = sum( Scores[i,1:k_PC]^2/eigvals_i[1:k_PC] )
    	Ti_Choi = sum( Scores[i,1:k_Choi]^2/sqrt(eigvals_i[1:k_Choi]) )
    	pvalmat$Ti_L2[i] = Ti_L2
    	pvalmat$Ti_PC[i] = Ti_PC
    	pvalmat$Ti_Choi[i] = Ti_Choi
    	pvalmat$pval_L2[i] = imhof(Ti_L2/eigvals_i[1], eigvals_i[1:k_L2]/eigvals_i[1])$Qq
    	pvalmat$pval_PC[i] = pchisq(Ti_PC, df = k_PC, lower.tail=FALSE)
    	pvalmat$pval_Choi[i] = imhof(Ti_Choi/sqrt(eigvals_i[1]), sqrt(eigvals_i[1:k_Choi])/sqrt(eigvals_i[1]),epsabs=10^(-16))$Qq

  	}

  	return(list(Bhat = Bhat, pvals = pvalmat, lambda=penlamb))
}

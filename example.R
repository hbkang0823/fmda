# Before running anything, download 'fdaM.zip' 
# from http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/ 
# and unzip it. 
# It is important to include the codes in './fdaM' for the code to run.

# load necessary packages and codes
source("setup.R")


####### Generate Faces #######
set.seed(1234)
N <- 100
sim <- simulate_face(N, delta = 200, err.sd = 0.001)
	# delta: size of effect
	# err.sd: err ~ N(0, err.sd)
ls(sim) # X and Y
dim(sim$Y) # N x 3 x P(=7150)

# plot the first face
n=1; plot_face(t(sim$Y[n,,]))


####################### Getting Mesh ##########################
# mean as refernece face, do manifold learning on reference face, and get mesh
meanface <- cbind(colMeans(sim$Y[,1,]), colMeans(sim$Y[,2,]), colMeans(sim$Y[,3,]))

wdir <- getwd() # should be the directory where example.R is.
filename <- paste0(wdir, "/rundata/for_dimred.mat")
writeMat(filename, data = meanface, d=2, type = 'ltsa')
# possible types: 'pca', 'ltsa', 'lle', 'hlle' 'laplacian', 'diffusemap'
run_matlab_script(paste0(wdir,'/function/dimred.m')) # this MATLAB code takes quite some time

# NEED TO WAIT until './rundata/from_dimred.mat' is created.
M0 <- readMat('./rundata/from_dimred.mat')$M0
mesh <- create_mesh(M0, ifplot=T)


######################## fdobj + FPCA #########################
# MATLAB: turning into fdobj using felsplines and first step of FPCA
# connection established through package 'matlabr' and 'R.matlab'
wdir <- getwd() # should be the directory where example.R is.
filename <- paste0(wdir, "/rundata/for_fpca1.mat")
writeMat(filename, pp=mesh$p, tt=mesh$t, data = sim$Y, M0 = M0, lamb_all = c(0, 10^c(-8:-5)))
run_matlab_script(paste0(wdir,'/function/fpca1.m')) # this MATLAB code takes quite some time

# NEED TO WAIT until './rundata/from_fpca1.mat' is created.
result <- readMat('./rundata/from_fpca1.mat')
ls(result)
# includes 'fpca1varprop', 'evalmat_x_mean', 'evalmat_y_mean', 'evalmat_z_mean', 
# 'harmcoef', 'coef', 'evalmat', 'PsiLap'

# R: second step of FPCA
SVD <- fpca2(result$coef)

# cutoff percentage for FPCA
EV_fpca <- 0.99	
varprop1 <- result$fpca1varprop
nharms1 <- dim(result$evalmat)[2]
nharms2 <- 1+sum((varprop1*(cumsum(SVD$d^2)/sum(SVD$d^2))) < EV_fpca)

# PC plot
orthVec <- get_orth_vec(meanface) # may take a while
k=1; plot_PC(SVD$v[,,k], meanface, result$evalmat, orthVec) # plot first PC


############# Penalized Manifold-on-Scalar Regression ############
U <- get_penalty(result, SVD, nharms2)

coef <- to.tensor(c(result$coef), c(N = N, I1 = 3, J1 = nharms1))
Ymat <- matrix(0, nrow=N, ncol=nharms2)
Ymat[1:N,1:nharms2] <- mul.tensor(coef, c("J1", "I1"), SVD$v, c("J", "I"))[1:N, 1:nharms2]
penmat <- matrix(c(U), nrow=dim(U)[1])
fit <- MS_penreg(Ymat, as.matrix(sim$X), penmat = penmat, penlamb_all = c(0,10^(-5:10)))

fit$pvals # shows the result based on Norm test, PC test, and Choi test.


############# Plots Showing the Significance of Beta #############
# Beta plot
plot_beta(fit$Bhat[1,], SVD, nharms2, meanface, result$evalmat, orthVec)
# compare this to the original beta
inprod_vec_beta <- diag(cbind(beta[(3*c(1:7150)-2),], beta[(3*c(1:7150)-1),], beta[(3*c(1:7150)-0),]) %*% t(orthVec))
rglplot(inprod_vec_beta, meanface, qbreaks=c(-1,0,1), color=c("blue", "red"))

# Pointwise significance plot
plot_ptsig(fit$Bhat[1,], sig.level=.95, Ymat, as.matrix(sim$X), SVD, meanface, result$evalmat, orthVec) # takes a bit of time
# Test is done pointwise.
# This plot includes the ovelapping area of fitted beta plot and original beta plot.

# Overall significance plot
plot_allsig(fit$Bhat[1,], sig.level=.95, Ymat, as.matrix(sim$X), SVD, meanface, result$evalmat, orthVec) # takes a bit of time
# Test is done for overall face -- alpha level (1-sig.level) is maintained over the whole face.
# This plot includes the colored part of pointwise significance plot.


####### Delete files in 'rundata' folder after analysis is done #######
file.remove("./rundata/for_dimred.mat")
file.remove("./rundata/from_dimred.mat")
file.remove("./rundata/for_fpca1.mat")
file.remove("./rundata/from_fpca1.mat")


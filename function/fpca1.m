addpath('..') 
addpath('../rundata')
addpath('../fdaM')
addpath('../function')

load('for_fpca1.mat')
% includes data, M0, pp, tt, lamb_all

datax = squeeze(data(:,1,:))';
datay = squeeze(data(:,2,:))';
dataz = squeeze(data(:,3,:))';

N = size(datax, 2);
P = size(datax, 1);
test = randsample(1:N, ceil(N/10));
testx = datax(:,test);
testy = datay(:,test);
testz = dataz(:,test);

basisobj = create_FEM_basis(pp, [], tt, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SSEvec = zeros(length(lamb_all), 1);
for j = 1:length(lamb_all)
	tmp_fdobj_x = smooth_FEM_basis(M0, testx, basisobj, lamb_all(j));
    tmp_fdobj_y = smooth_FEM_basis(M0, testy, basisobj, lamb_all(j));
	tmp_fdobj_z = smooth_FEM_basis(M0, testz, basisobj, lamb_all(j));

    evalmat_x = eval_FEM_fd(M0, tmp_fdobj_x);
	evalmat_y = eval_FEM_fd(M0, tmp_fdobj_y);
	evalmat_z = eval_FEM_fd(M0, tmp_fdobj_z);
    
    xdiff = evalmat_x - testx;
	ydiff = evalmat_y - testy;
	zdiff = evalmat_z - testz;

    SSEvec(j) = sum(sum(xdiff.^2 + ydiff.^2 + zdiff.^2));
end

[minval, whichmin] = min(SSEvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdobj_x = smooth_FEM_basis(M0, datax, basisobj, lamb_all(whichmin));
fdobj_y = smooth_FEM_basis(M0, datay, basisobj, lamb_all(whichmin));
fdobj_z = smooth_FEM_basis(M0, dataz, basisobj, lamb_all(whichmin));
    
if (3*N) < 200
    nharm = floor(3*N*0.9);
else
    nharm = 200;
end

pcastr = pca_fd_pooled(fdobj_x, fdobj_y, fdobj_z, nharm);

fpca1varprop = sum(pcastr.varprop);     %save
evalmat_x_mean = eval_FEM_fd(M0, pcastr.meanfd1); %save
evalmat_y_mean = eval_FEM_fd(M0, pcastr.meanfd2); %save
evalmat_z_mean = eval_FEM_fd(M0, pcastr.meanfd3); %save

fdobj1 = center(fdobj_x);
fdobj2 = center(fdobj_y);
fdobj3 = center(fdobj_z);

coef1 = getcoef(fdobj1);
coef2 = getcoef(fdobj2);
coef3 = getcoef(fdobj3);
ctemp = [coef1'; coef2'; coef3']';

Jmat = inprod_FEM_basis(basisobj, basisobj);
harmfd = pcastr.harmfd;

harmcoef = getcoef(harmfd); %save
coef = ctemp'*Jmat*getcoef(harmfd); %save
evalmat = eval_FEM_fd(M0, harmfd);    %save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params   = getbasispar(basisobj);
nodeStruct.order     = params.order;
nodeStruct.nodes     = params.nodes;
nodeStruct.nodeindex = params.nodeindex;
nodeStruct.J         = params.J;
nodeStruct.metric    = params.metric;

K0 = mass(nodeStruct);
%same as Jmat = inprod_FEM_basis(basisobj, basisobj);

K1 = stiff(nodeStruct);
Fmat = K1 * harmcoef;
PsiLap = Fmat' * K0 * Fmat;  %save

save('../rundata/from_fpca1', 'fpca1varprop', 'evalmat_x_mean', 'evalmat_y_mean', 'evalmat_z_mean', 'harmcoef', 'coef', 'evalmat', 'PsiLap');


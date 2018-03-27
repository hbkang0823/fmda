function pcastr = pca_fd_pooled(fdobj1, fdobj2, fdobj3, nharm, harmfdPar, centerfns)
%  PCA Functional principal components analysis with regularization
%
%  Arguments:
%  FDOBJ     ... Functional data object (a struct object)
%  NHARM     ... Number of principal components to be kept. Default 2
%  HARMFDPAR ... A functional parameter object specifying the
%                basis, differential operator, and level of smoothing
%                for eigenfunctions.
%  CENTERFNS ... If 1, the mean function is first subtracted from each 
%                function.  1 is the default.
%
%  Returns:
%  A struct object PCASTR with the fields:
%  HARMFD  ... A functional data object for the harmonics or eigenfunctions
%  VALUES  ... The complete set of eigenvalues
%  HARMSCR ... A matrix of scores on the principal components or harmonics
%  VARPROP ... A vector giving the proportion of variance explained
%                 by each eigenfunction
%  FDHATFD ... A functional data object for the approximation to the
%              FDOBJ based on NHARM principal components
%  MEANFD  ... A functional data object giving the mean function
%
%  If NHARM = 0, all fields except MEANFD are empty.

%  Last modified:  29 January 2013 by Jim Ramsay

%  check FDOBJ

if ~isa_fd(fdobj1) || ~isa_fd(fdobj2) || ~isa_fd(fdobj3)
    error ('One of the first three arguments is not a functional data object.');
end

%  get basis information for functional data

fdbasis1  = getbasis(fdobj1);
fdbasis2  = getbasis(fdobj2);
fdbasis3  = getbasis(fdobj3);

% if they are all same
fdbasis = fdbasis1;

%  set up default values

if nargin < 6
    centerfns = 1;   %  subtract mean from data before PCA
end

if nargin < 5
    %  default Lfd object: penalize 2nd deriv., lambda = 0
    Lfdobj    = int2Lfd(2);
    lambda    = 0;
    harmfdPar = fdPar(fdbasis, Lfdobj, lambda);
else
    %  check harmfdPar object
    if ~isa_fdPar(harmfdPar)
        if isa_fd(harmfdPar) || isa_basis(harmfdPar)
            harmfdPar = fdPar(harmfdPar);
        else
            error(['HARMFDPAR is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
        end
    end
end

if nargin < 4
    nharm = 2;  %  default to two harmonics
end

%  compute mean function

meanfd1 = mean(fdobj1);
meanfd2 = mean(fdobj2);
meanfd3 = mean(fdobj3);

if nharm == 0
    pcastr.harmfd  = [];
    pcastr.values  = [];
    pcastr.harmscr = [];
    pcastr.varprop = [];
    pcastr.fdhatfd = [];
    pcastr.meanfd  = meanfd;
    return
end

%  --------------------   begin principal components analysis  ------------

% center data if required

if centerfns ~= 0
    fdobj1 = center(fdobj1);
    fdobj2 = center(fdobj2);
    fdobj3 = center(fdobj3);
end

%  set up HARMBASIS

harmbasis = getbasis(getfd(harmfdPar));
nhbasis   = getnbasis(harmbasis);

%  set up LFDOBJ

Lfdobj = getLfd(harmfdPar);
Lfdobj = int2Lfd(Lfdobj);

%  set up LAMBDA

lambda = getlambda(harmfdPar);

%  get coefficient matrix for FDOBJ and its dimensions

coef1   = getcoef(fdobj1);
coef2   = getcoef(fdobj2);
coef3   = getcoef(fdobj3);

coefd1  = size(coef1);
coefd2  = size(coef2);
coefd3  = size(coef3);

nbasis1 = coefd1(1);
nrep1   = coefd1(2);
ndim1   = length(coefd1);

nbasis2 = coefd2(1);
nrep2   = coefd2(2);
ndim2   = length(coefd2);

nbasis3 = coefd3(1);
nrep3   = coefd3(2);
ndim3   = length(coefd3);

if nrep1 < 2 || nrep2 < 3 || nrep3 < 3
    error('PCA not possible without replications');
end

%  compute CTEMP whose cross product is needed

% if ndim == 3
%     nvar  = coefd(3);
%     ctemp = zeros(nvar*nbasis,nrep);
%     for j = 1:nvar
%         index = (1:nbasis) + (j-1)*nbasis;
%         ctemp(index,:) = coef(:,:,j);
%     end
% else
    nvar = 1;
    ctemp = [coef1'; coef2'; coef3']';
    nbasis = nbasis1;
    nrep = size(ctemp, 2);
    ndim = ndim3;

%  set up cross product Lmat for harmonic basis, 
%  roughness penalty matrix Rmat, and
%  penalized cross product matrix Lmat.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% saw until this line
%Lmat = eval_penalty(harmbasis, int2Lfd(0));
% with the assumption that we use 'FEM' type, use 'FEMpen_rv'.
Lmat = FEMpen_rv(harmbasis, int2Lfd(0));
if lambda > 0
%    Rmat = eval_penalty(harmbasis, Lfdobj);
    Rmat = FEMpen_rv(harmbasis, Lfdobj);
    Lmat = Lmat + lambda .* Rmat;
end
Lmat = (Lmat + Lmat')/2;

%  Choleski factor Mmat of Lmat = Mmat'*Mmat

Mmat    = chol(Lmat);
Mmatinv = inv(Mmat);

%  coefficient cross product matrix for covariance operator

%Wmat1 = ctemp1*ctemp1'./nrep1;
%Wmat2 = ctemp2*ctemp2'./nrep2;
%Wmat3 = ctemp3*ctemp3'./nrep3;

%Wmat = blkdiag(Wmat1, Wmat2, Wmat3);

% OR
% ctemp = [ctemp1; ctemp2; ctemp3];
% Wmat = ctemp*ctemp'./nrep;

Wmat = ctemp*ctemp'./nrep;

%  set up matrix for eigenanalysis depending on whether
%  a special basis was supplied for the eigenfunctions or not

Jmat = inprod_FEM_basis(harmbasis, fdbasis);
%Jmat = blkdiag(Jmat, Jmat, Jmat);

MIJW = Mmatinv'*Jmat;

if nvar == 1
    Cmat = MIJW*Wmat*MIJW';
else
    Cmat = zeros(nvar*nhbasis);
    for i = 1:nvar
        indexi =   (1:nbasis) + (i-1)*nbasis;
        for j = 1:nvar
            indexj = (1:nbasis) + (j-1)*nbasis;
            Cmat(indexi,indexj) = MIJW*Wmat(indexi,indexj)*MIJW';
        end
    end
end

% Eigenanalysis

Cmat = (Cmat + Cmat')./2;
[eigvecs, eigvals] = eig(Cmat);
[eigvals, indsrt ] = sort(diag(eigvals));
eigvecs = eigvecs(:,indsrt);
neig    = nvar*nhbasis;
indx    = neig + 1 - (1:nharm);
eigvals = eigvals(neig + 1 - (1:neig));
eigvecs = eigvecs(:,indx);
sumvecs = sum(eigvecs);
eigvecs(:,sumvecs < 0) = -eigvecs(:,sumvecs < 0);

varprop = eigvals(1:nharm)./sum(eigvals);

%  Set up fdnames for harmfd

harmnames = getnames(fdobj1);
%  Name and labels for harmonics
harmlabels = ['I   '; 'II  '; 'III '; 'IV  '; 'V   '; ...
              'VI  '; 'VII '; 'VIII'; 'IX  '; 'X   '];
if nharm <= 10
    harmnames2    = cell(1,2);
    harmnames2{1} = 'Harmonics';
    harmnames2{2} = harmlabels(1:nharm,:);
    harmnames{2}  = harmnames2;
else
    harmnames{2} = 'Harmonics';
end
%  Name and labels for variables
if iscell(harmnames{3})
    harmnames3    = harmnames{3};
    harmnames3{1} = ['Harmonics for ',harmnames3{1}];
    harmnames{3}  = harmnames3;
else
    if ischar(harmnames{3}) && size(harmnames{3},1) == 1
        harmnames{3} = ['Harmonics for ',harmnames{3}];
    else
        harmnames{3} = 'Harmonics';
    end
end

%  set up harmfd

if nvar == 1
    %  harmcoef = Lmat\eigvecs; changed 20 Nov 13
    %  see p 38 of AFDA
    harmcoef = Mmat\eigvecs;
else
    harmcoef = zeros(nbasis,nharm,nvar);
    for j = 1:nvar
        index     = (1:nbasis) + (j-1)*nbasis;
        eigvecsj  = eigvecs(index,:);
        %  harmcoef(:,:,j) = Lmat\eigvecsj; changed 20 Nov 13
        harmcoef(:,:,j) = Mmat\eigvecsj;
    end
end

harmfd = fd(harmcoef, harmbasis, harmnames);
% harmfd = fd(harmcoef, fdbasis);


%  set up harmscr

% if nvar == 1
%     if strcmp(getbasistype(fdbasis), 'FEM')     %added this part!!!!!!
        harmscr1 = inprod_FEM_fd(fdobj1, harmfd);
        harmscr2 = inprod_FEM_fd(fdobj2, harmfd);
        harmscr3 = inprod_FEM_fd(fdobj3, harmfd);
        
%     else
%         harmscr = inprod(fdobj, harmfd);
%     end
% else
%     harmscr       = zeros(nrep, nharm, nvar);
%     coefarray     = getcoef(fdobj);
%     harmcoefarray = getcoef(harmfd);
%     for j=1:nvar
%         coefj     = squeeze(coefarray(:,:,j));
%         harmcoefj = squeeze(harmcoefarray(:,:,j));
%         fdobjj    = fd(coefj, fdbasis);
%         harmfdj   = fd(harmcoefj, fdbasis);
%         if strcmp(getbasistype(fdobjj), 'FEM')     %added this part!!!!!!
%             harmscr(:,:,j) = inprod_FEM_fd(fdobjj, harmfdj);
%         else
%             harmscr(:,:,j) = inprod(fdobjj,harmfdj);
%         end
%     end
% end

%  set up functional data object for fit to the data

% if nvar == 1
    fdhatcoef1 = harmcoef*harmscr1';
    fdhatcoef2 = harmcoef*harmscr2';
    fdhatcoef3 = harmcoef*harmscr3';
% else
%     fdhatcoef = zeros(nbasis,nrep,nvar);
%     for j=1:nvar
%         fdhatcoef(:,:,j) = harmcoef(:,:,j)*harmscr(:,:,j)';
%     end
% end
fdhatfd1 = fd(fdhatcoef1, harmbasis);
fdhatfd2 = fd(fdhatcoef2, harmbasis);
fdhatfd3 = fd(fdhatcoef3, harmbasis);

%  set up structure object PCASTR

pcastr.harmfd  = harmfd;
pcastr.values  = eigvals;
pcastr.harmscr1 = harmscr1;
pcastr.harmscr2 = harmscr2;
pcastr.harmscr3 = harmscr3;
pcastr.varprop = varprop;
pcastr.fdhatfd1 = fdhatfd1;
pcastr.fdhatfd2 = fdhatfd2;
pcastr.fdhatfd3 = fdhatfd3;
pcastr.meanfd1  = meanfd1;
pcastr.meanfd2  = meanfd2;
pcastr.meanfd3  = meanfd3;


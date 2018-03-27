% reference:
%   mani: MANIfold learning demonstration GUI
%   by Todd Wittman, Department of Mathematics, University of Minnesota
%   E-mail wittman@math.ucla.edu with comments & questions.
%   MANI Website: http://www.math.ucla.edu/~wittman/mani/index.html
%   Last Modified by GUIDE v2.5 10-Apr-2005 13:28:36

%   Methods obtained from various authors.
%      MDS -- Michael Lee
%      ISOMAP -- J. Tenenbaum, de Silva, & Langford
%      LLE -- Sam Roweis & Lawrence Saul
%      Hessian LLE  -- D. Donoho & C. Grimes
%      Laplacian -- M. Belkin & P. Niyogi
%      Diffusion Map -- R. Coifman & S. Lafon
%      LTSA -- Zhenyue Zhang & Hongyuan Zha

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- HLLE function
% Written by David Donoho & Carrie Grimes, 2003.
function [Y, mse] = HLLE(X,k,d)
N = size(X,2);
if max(size(k)) ==1
    kvec = repmat(k,N,1);
elseif max(size(k)) == N
    kvec=k;
end;
%Compute Nearest neighbors
D1 = L2_distance(X,X,1);
dim = size(X,1);
nind = repmat(0, size(D1,1), size(D1,2));
dp = d*(d+1)/2;
W = repmat(0,dp*N,N);
if(mean(k)>d) 
  tol=1e-3; % regularlizer in case constrained fits are ill conditioned
else
  tol=0;
end;
for i=1:N
    tmp = D1(:,i);
    [ts, or] = sort(tmp);
    %take k nearest neighbors
    nind(or(2:kvec(i)+1),i) = 1;
    thisx = X(:,or(2:kvec(i)+1));
    %center using the mean 
    thisx = thisx - repmat(mean(thisx')',1,kvec(i));
    %compute local coordinates
    [U,D,Vpr] = svd(thisx);
    V = Vpr(:,1:d);
    %Neighborhood diagnostics
    vals = diag(D);
    mse(i) = sum(vals(d+1:end));
    %build Hessian estimator
    clear Yi; clear Pii;
    ct = 0;
    for mm=1:d
        startp = V(:,mm);
        for nn=1:length(mm:d)
            indles = mm:d;
            Yi(:,ct+nn) = startp.*(V(:,indles(nn)));
        end;
        ct = ct+length(mm:d);
    end;
    Yi = [repmat(1,kvec(i),1), V, Yi];
    %orthogonalize linear and quadratic forms
    [Yt, Orig] = mgs(Yi);
    Pii = Yt(:,d+2:end)';
    %double check weights sum to 1
    for j=1:dp
        if sum(Pii(j,:)) >0.0001
            tpp = Pii(j,:)./sum(Pii(j,:)); 
        else
            tpp = Pii(j,:);
        end;
        %fill weight matrix
       W((i-1)*dp+j, or(2:kvec(i)+1)) = tpp;
    end;
end;
%%%%%%%%%%%%%%%%%%%%Compute eigenanalysis of W
G=W'*W;
G = sparse(G);
options.disp = 0; 
options.isreal = 1; 
options.issym = 1;
tol=0;
[Yo,eigenvals] = eigs(G,d+1,tol,options);
Y = Yo(:,1:d)'*sqrt(N); % bottom evect is [1,1,1,1...] with eval 0
%compute final coordinate alignment
R = Y'*Y;
R2 = R^(-1/2);
Y = Y*R2;
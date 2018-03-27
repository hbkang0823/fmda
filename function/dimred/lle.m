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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- LLE ALGORITHM (using K nearest neighbors)
% Written by Sam Roweis & Lawrence Saul
function [Y] = lle(X,K,d)
warning off;
[D,N] = size(X);
% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);
% STEP2: SOLVE FOR RECONSTRUCTION WEIGHTS
if(K>D) 
  tol=1e-3; % regularlizer in case constrained fits are ill conditioned
else
  tol=0;
end
W = zeros(K,N);
for ii=1:N
   z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
   C = z'*z;                                        % local covariance
   C = C + eye(K,K)*tol*trace(C);                   % regularlization (K>D)
   W(:,ii) = C\ones(K,1);                           % solve Cw=1
   W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
end;
% STEP 3: COMPUTE EMBEDDING FROM EIGENVECTS OF COST MATRIX M=(I-W)'(I-W)
M = sparse(1:N,1:N,ones(1,N),N,N,4*K*N); 
for ii=1:N
   w = W(:,ii);
   jj = neighborhood(:,ii);
   M(ii,jj) = M(ii,jj) - w';
   M(jj,ii) = M(jj,ii) - w;
   M(jj,jj) = M(jj,jj) + w*w';
end;
% CALCULATION OF EMBEDDING
options.disp = 0; options.isreal = 1; options.issym = 1; 
[Y,eigenvals] = eigs(M,d+1,0,options);
Y = Y(:,1:d)'*sqrt(N);   % bottom evect is [1,1,1,1...] with eval 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- ISOMAP Algorithm
% Written by Tenenbaum, de Silva, and Langford (2000). 
function [Y, R] = Isomap(D, K, options); 
%%%%% Step 0: Initialization and Parameters %%%%%
N = size(D,1); 
INF =  1000*max(max(D))*N;  %% effectively infinite distance
dims = options.dims; 
comp = 1;  % Assume one component.
Y.coords = cell(length(dims),1); 
R = zeros(1,length(dims)); 
%%%%% Step 1: Construct neighborhood graph %%%%%
[tmp, ind] = sort(D); 
for i=1:N
    D(i,ind((2+K):end,i)) = INF; 
end;
D = min(D,D');    %% Make sure distance matrix is symmetric
%%%%% Step 2: Compute shortest paths %%%%%
for k=1:N
     D = min(D,repmat(D(:,k),[1 N])+repmat(D(k,:),[N 1])); 
end
%%%%% Remove outliers from graph %%%%%
n_connect = sum(~(D==INF));        %% number of points each point connects to
[tmp, firsts] = min(D==INF);       %% first point each point connects to
[comps, I, J] = unique(firsts);    %% represent each connected component once
size_comps = n_connect(comps);     %% size of each connected component
[tmp, comp_order] = sort(size_comps);  %% sort connected components by size
comps = comps(comp_order(end:-1:1));    
size_comps = size_comps(comp_order(end:-1:1)); 
n_comps = length(comps);               %% number of connected components
if (comp>n_comps)                
     comp=1;                              %% default: use largest component
end
Y.index = find(firsts==comps(comp)); 
D = D(Y.index, Y.index); 
N = length(Y.index); 
%%%%% Step 3: Construct low-dimensional embeddings (Classical MDS) %%%%%
opt.disp = 0; 
[vec, val] = eigs(-.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2)), max(dims), 'LR', opt); 
h = real(diag(val)); 
[foo,sorth] = sort(h);  sorth = sorth(end:-1:1); 
val = real(diag(val(sorth,sorth))); 
vec = vec(:,sorth); 
D = reshape(D,N^2,1); 
for di = 1:length(dims)
     if (dims(di)<=N)
         Y.coords{di} = real(vec(:,1:dims(di)).*(ones(N,1)*sqrt(val(1:dims(di)))'))'; 
         r2 = 1-corrcoef(reshape(real(L2_distance(Y.coords{di}, Y.coords{di},0)),N^2,1),D).^2; 
         R(di) = r2(2,1); 
     end
end
clear D; 
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
% --- diffusionKernel function
% Written by R. Coifman & S. Lafon.
function [Y] = diffusionKernel (X,sigmaK,alpha,d)
D = L2_distance(X',X',1);
K = exp(-(D/sigmaK).^2);
p = sum(K);
p = p(:);
K1 = K./((p*p').^alpha);
v = sqrt(sum(K1));
v = v(:);
A = K1./(v*v');
if sigmaK >= 0.5
    thre = 1e-7;  
    M = max(max(A));
    A = sparse(A.*double(A>thre*M));
    [U,S,V] = svds(A,d+1);   %Sparse version.
    U = U./(U(:,1)*ones(1,d+1));
else
    [U,S,V] = svd(A,0);   %Full version.
    U = U./(U(:,1)*ones(1,size(U,1)));
end;
Y = U(:,2:d+1);
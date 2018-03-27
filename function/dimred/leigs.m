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
% --- leigs function for Laplacian eigenmap.
% Written by Belkin & Niyogi, 2002.
function [E,V] = leigs(DATA, TYPE, PARAM, NE) 
n = size(DATA,1);
A = sparse(n,n);
step = 100;  
for i1=1:step:n    
    i2 = i1+step-1;
    if (i2> n) 
      i2=n;
    end;
    XX= DATA(i1:i2,:);  
    dt = L2_distance(XX',DATA',0);
    [Z,I] = sort ( dt,2);
    for i=i1:i2
      for j=2:PARAM+1
	        A(i,I(i-i1+1,j))= Z(i-i1+1,j); 
	        A(I(i-i1+1,j),i)= Z(i-i1+1,j); 
      end;    
    end;
end;
W = A;
[A_i, A_j, A_v] = find(A);  % disassemble the sparse matrix
for i = 1: size(A_i)  
    W(A_i(i), A_j(i)) = 1;
end;
D = sum(W(:,:),2);   
L = spdiags(D,0,speye(size(W,1)))-W;
opts.tol = 1e-9;
opts.issym=1; 
opts.disp = 0; 
[E,V] = eigs(L,NE,'sm',opts);
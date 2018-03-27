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
% --- function mdsFast for Multi-Dimensional Scaling
% Written by Michael D. Lee.
% Lee recommends metric=2, iterations=50, learnrate=0.05.
function [points]=mdsFast(d,dim)
[n, check] = size(d);
iterations = 30;
lr=0.05;   % learnrate
r=2;   % metric
% normalise distances to lie between 0 and 1
reshift=min(min(d));
d=d-reshift;
rescale=max(max(d));
d=d/rescale;
% calculate the variance of  the distance matrix
dbar=(sum(sum(d))-trace(d))/n/(n-1);
temp=(d-dbar*ones(n)).^2;
vard=.5*(sum(sum(temp))-trace(temp));
% initialize variables
its=0;
p=rand(n,dim)*.01-.005;
dh=zeros(n);
rinv=1/r;  %PT
kk=1:dim;  %PT 
% main loop for the given number of iterations
while (its<iterations)
   its=its+1;
   % randomly permute the objects to determine the order
   % in which they are pinned for this iteration
   pinning_order=randperm(n);
   for i=1:n
      m=pinning_order(i);
      % having pinned an object, move all of the other on each dimension
      % according to the learning rule   
      
      %PT: Vectorized the procedure, gives factor 6 speed up for n=100 dim=2
      indx=[1:m-1 m+1:n];                                                       
      pmat=repmat(p(m,:),[n 1])-p;                                              
      dhdum=sum(abs(pmat).^r,2).^rinv;
      dh(m,indx)=dhdum(indx)';
      dh(indx,m)=dhdum(indx);
      dhmat=lr*repmat((dhdum(indx)-d(m,indx)').*(dhdum(indx).^(1-r)),[1 dim]);
      p(indx,kk)=p(indx,kk)+dhmat.*abs(pmat(indx,kk)).^(r-1).*sign(pmat(indx,kk));
                    %plus sign in learning rule is due the sign of pmat
   end;
end;
points = p*rescale+reshift;
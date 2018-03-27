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
% --- L2_distance function
% Written by Roland Bunschoten, University of Amsterdam, 1999
function d = L2_distance(a,b,df)
if (size(a,1) == 1)
  a = [a; zeros(1,size(a,2))]; 
  b = [b; zeros(1,size(b,2))]; 
end
aa=sum(a.*a); bb=sum(b.*b); ab=a'*b; 
d = sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
d = real(d); 
if (df==1)
  d = d.*(1-eye(size(d)));
end
addpath('..') 
addpath('../rundata')
addpath('../fdaM')
addpath('../function')
addpath('../function/dimred')

load('for_dimred.mat') %need to include data, d, type

if ~exist('k')
  k = 8;
end
if ~exist('s')
  s = 10;
end
if ~exist('a')
  a = 1;
end

if strcmp(type, 'pca')
    tmpY = pca(data, d);
elseif strcmp(type, 'ltsa')
    [T, NI] = LTSA(data', d, k);
    tmpY = T';
elseif strcmp(type, 'll3')
    tmpY = lle(data', k, d)';
elseif strcmp(type, 'hlle')
    [Y, mse] = HLLE(data', k, d);
    tmpY = Y';
elseif strcmp(type, 'laplacian')
    [E, V] = leigs(data, 'nn', k, d+1);
    tmpY = E(:,1:d);
elseif strcmp(type, 'diffusemap')
    tmpY = diffusionKernel(data, s, a, d);
end

M0 = tmpY;
save('../rundata/from_dimred', 'M0');
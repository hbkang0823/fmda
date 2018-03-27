function  ss = inprod_FEM_basis(basis1, basis2, nrand)
%  INPROD_BASIS  Computes matrix of inner products of FEM bases.
%  NRAND:   number of uniform random points per triangle for approximation
%  last modified 1 Mar 2016
%  Hyun Bin Kang

    if ~(isa_basis(basis1) && isa_basis(basis2))
        error ('The two first arguments are not basis objects.');
    end

    params1 = getbasispars(basis1);
    params2 = getbasispars(basis2);

    if sum(sum(params1.p ~= params2.p)) + sum(sum(params1.t ~= params2.t)) ~= 0
        error('mesh is different');
    end
    
    params = params1;
    
    nodeStruct.order     = params.order;        %
	nodeStruct.nodes     = params.nodes;        %
	nodeStruct.nodeindex = params.nodeindex;    %
	nodeStruct.J         = params.J;            %
	nodeStruct.metric    = params.metric;       %
    
	ss = mass(nodeStruct);
        
% %     params = params1;
% %     ntri = size(params.t, 1);
% %     nbasis = getnbasis(basis1);
% %     tri_basismat1 = zeros(ntri, getnbasis(basis1));
% %     tri_basismat2 = zeros(ntri, getnbasis(basis2));
% %     trivols = zeros(ntri,1);
% % 
% %     if nargin < 3
% %         nrand = 100;
% %     end
% %     tmp_pts = zeros(nrand,2);
% %     
% % %    r1 = linspace(0,1,nrand);
% % %    r2 = linspace(0,1,nrand);
% % 
% %     for itri = 1:ntri
% %         x1 = params.p(params.t(itri,1),1); y1 = params.p(params.t(itri,1),2);
% %         x2 = params.p(params.t(itri,2),1); y2 = params.p(params.t(itri,2),2);
% %         x3 = params.p(params.t(itri,3),1); y3 = params.p(params.t(itri,3),2);
% % 
% %         r1 = rand(nrand,1);
% %         r2 = rand(nrand,1);
% %         tmp_pts(:,1) = (1-sqrt(r1))*x1 + sqrt(r1).*(1-r2)*x2 + sqrt(r1).*r2*x3;
% %         tmp_pts(:,2) = (1-sqrt(r1))*y1 + sqrt(r1).*(1-r2)*y2 + sqrt(r1).*r2*y3;
% % 
% %         tmp_triptavg1 = eval_FEM_basis(tmp_pts(:,1), tmp_pts(:,2), basis1);
% %         tmp_triptavg1 = (1/nrand) * sum(tmp_triptavg1);
% %         tri_basismat1(itri,:) = tmp_triptavg1;
% % 
% %         tmp_triptavg2 = eval_FEM_basis(tmp_pts(:,1), tmp_pts(:,2), basis2);
% %         tmp_triptavg2 = (1/nrand) * sum(tmp_triptavg2);
% %         tri_basismat2(itri,:) = tmp_triptavg2;
% %         
% %         trivols(itri) = 1/2 * abs(x1*y2 + x2*y3 + x3*y1 - x2*y1 - x3*y2 - x1*y3);
% %     end
% % 
% %     ss = tri_basismat1' * diag(trivols) * tri_basismat2;
    
end
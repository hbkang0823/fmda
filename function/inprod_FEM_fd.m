function  inprodmat = inprod_FEM_fd(fdobj1, fdobj2)
%  INPROD   Computes matrix of inner products of FEM fd objects.
%  last modified 11 Feb 2016
%  Hyun Bin Kang

%  FDOBJ1

    if isa_fd(fdobj1)
        coef1  = getcoef(fdobj1);
        coefd1 = size(coef1);
        if length(coefd1) > 2
            error('Functional data object must be univariate');
        end
        nrep1     = coefd1(2);
        basisobj1 = getbasis(fdobj1);
        type1     = getbasistype(basisobj1);
        params1   = getbasispars(basisobj1);
    else
        error('FDOBJ1 is not a functional data object.');
    end

    %  FDOBJ2

    if isa_fd(fdobj2)
        coef2  = getcoef(fdobj2);
        coefd2 = size(coef2);
        if length(coefd2) > 2
            error('Functional data object must be univariate');
        end
        nrep2     = coefd2(2);
        basisobj2 = getbasis(fdobj2);
        type2     = getbasistype(basisobj2);
        params2   = getbasispars(basisobj2);
    end

    if sum(sum(params1.p ~= params2.p)) + sum(sum(params1.t ~= params2.t)) ~= 0
        error('mesh is different');
    end

    Jmat = inprod_FEM_basis(basisobj1, basisobj2);
    inprodmat = coef1'*Jmat*coef2;
end
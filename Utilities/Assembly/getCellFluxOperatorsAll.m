function op = getCellFluxOperatorsAll(G, varargin)
% Setup the operators P and S that are necessary to compute norm of cell flux
%
% The operator P : facetbl (face-valued) to cellvecttbl (cell-valued vector)
% maps face-valued integrated fluxes  to a reconstruction of the vector flux in each cell
%
% The operator S :  cellvecttbl (cell-valued vector) to celltbl (cell-values scalar)
% simply sums up the component
%   
% Let u be a face-values integrated flux (in facetbl)
%   
% To obtain an approximation of the norm of corresponding flux
%
%      v = P*u;
%      v = v.^2;
%      v = S*v;
%      v = sqrt(v);
    %G=cartGrid([3,3,3])
    %G =computeGeometry(G);
    Neigh = G.faces.neighbors(G.cells.faces(:,1),:);
    cellNo=rldecode([1:G.cells.num]',diff(G.cells.facePos)');    
    %% 
    nc     = G.cells.num;
    nhf= numel(G.cells.faces(:,1));
    nf = G.faces.num;
    dims   = G.griddim;
    m  = diff(G.cells.facePos);
    n      = G.griddim*ones(nc, 1);
    [I, J] = blockDiagIndex(n, m);
    N      = G.faces.normals(G.cells.faces(:,1), :);
    sign   = 2*(cellNo == Neigh(:,1)) - 1;
    sN     = bsxfun(@times, N, sign);
    N      = sparse(I, J, reshape(sN', [], 1), dims*nc, nhf)';
    
    %% 
    opt = struct('invertBlocks', 'mex');
    opt = merge_options(opt, varargin{:});
    bi = blockInverter(opt);
    NTN = N'*N;
    NTNinv = bi(NTN, n);
    
    %% 

    f2hf = sparse((1 : nhf)', G.cells.faces(:,1), sign, nhf, nf);
    P = NTNinv*N'*f2hf;
    S = sparse(rldecode(1:nc, dims), [1 : dims*nc]', 1, nc, dims*nc);
    op = struct('S', S, 'P', P);
    
end
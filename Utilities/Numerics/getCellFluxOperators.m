function op = getCellFluxOperators(G)
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
    
    tbls = setupSimpleTables(G);

    celltbl = tbls.celltbl;
    facetbl = tbls.facetbl;
    cellfacetbl = tbls.cellfacetbl;

    N = G.faces.neighbors;
    intInx = all(N ~= 0, 2); % same definition as in setupOperatorsTPFA
    intfacetbl.faces = find(intInx);
    intfacetbl = IndexArray(intfacetbl);

    cellintfacetbl = crossIndexArray(cellfacetbl, intfacetbl, {'faces'});

    cn = cellintfacetbl.get('cells');
    cf = cellintfacetbl.get('faces');
    sgn = 2*double(G.faces.neighbors(cf, 1) == cn) - 1;

    vecttbl.vect = (1 : G.griddim)';
    vecttbl = IndexArray(vecttbl);

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = vecttbl;
    gen.tbl2 = vecttbl;
    gen.replacefds1 = {{'vect', 'vect1'}};
    gen.replacefds2 = {{'vect', 'vect2'}};
    gen.mergefds = {};

    vect12tbl = gen.eval();

    intfacevecttbl = crossIndexArray(intfacetbl, vecttbl, {}, 'optpureproduct', true);
    facevecttbl = crossIndexArray(facetbl, vecttbl, {}, 'optpureproduct', true);
    cellvecttbl = crossIndexArray(celltbl, vecttbl, {}, 'optpureproduct', true);
    cellvect12tbl = crossIndexArray(celltbl, vect12tbl, {}, 'optpureproduct', true);
    cellintfacevecttbl = crossIndexArray(cellintfacetbl, vecttbl, {}, 'optpureproduct', true);
    cellfacevecttbl = crossIndexArray(cellfacetbl, vecttbl, {}, 'optpureproduct', true);

    N = G.faces.normals;
    N = reshape(N', [], 1); % N is in facevecttbl

    map = TensorMap();
    map.fromTbl = facevecttbl;
    map.toTbl = intfacevecttbl;
    map.mergefds = {'faces', 'vect'};
    map = map.setup();

    N = map.eval(N); % N is in intfacevecttbl

    prod = TensorProd();
    prod.tbl1 = cellintfacetbl;
    prod.tbl2 = intfacevecttbl;
    prod.tbl3 = cellintfacevecttbl;
    prod.mergefds = {'faces'};
    prod = prod.setup();

    N = prod.eval(sgn, N); % N is in cellintfacevecttbl

    % We compute NtN

    prod = TensorProd();
    prod.tbl1 = cellintfacevecttbl;
    prod.tbl2 = cellintfacevecttbl;
    prod.tbl3 = cellvect12tbl;
    prod.replacefds1 = {{'vect', 'vect1'}};
    prod.replacefds2 = {{'vect', 'vect2'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'faces'};
    prod = prod.setup();

    NtN = prod.eval(N, N); % NtN is in cellvect12tbl

    %% We setup the tensor so that we can compute block inverse

    prod = TensorProd();
    prod.tbl1 = cellvect12tbl;
    prod.tbl2 = cellvecttbl;
    prod.tbl3 = cellvecttbl;
    prod.replacefds1 = {{'vect1', 'vect'}};
    prod.replacefds2 = {{'vect', 'vect2'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'vect2'};

    A = SparseTensor();
    A = A.setFromTensorProd(NtN, prod);
    A = A.getMatrix();

    % fetch block diagonal inverter
    opt.invertBlocks = 'matlab';
    bi = blockInverter(opt);
    sz = vecttbl.num*ones(celltbl.num, 1);
    invNtN = bi(A, sz);

    map = TensorMap();
    map.fromTbl = cellvecttbl;
    map.toTbl = cellvect12tbl;
    map.replaceFromTblfds = {{'vect', 'vect1'}};
    map.mergefds = {'cells', 'vect1'};
    ind1 = map.getDispatchInd();

    map = TensorMap();
    map.fromTbl = cellvecttbl;
    map.toTbl = cellvect12tbl;
    map.replaceFromTblfds = {{'vect', 'vect2'}};
    map.mergefds = {'cells', 'vect2'};
    ind2 = map.getDispatchInd();

    ind = sub2ind([cellvecttbl.num, cellvecttbl.num], ind1, ind2);

    invNtN = invNtN(ind); % invNtN is in cellvect12tbl

    % We set up mapping from face fluxes to signed fluxes on internal faces

    prod = TensorProd();
    prod.tbl1 = cellintfacevecttbl;
    prod.tbl2 = cellintfacetbl;
    prod.tbl3 = cellintfacevecttbl;
    prod.mergefds = {'cells', 'faces'};
    prod = prod.setup();

    sN = prod.eval(N, sgn); % sN is in cellintfacevecttbl

    % We compute P = (invNtN)*sN

    prod = TensorProd();
    prod.tbl1 = cellvect12tbl;
    prod.tbl2 = cellintfacevecttbl;
    prod.tbl3 = cellintfacevecttbl;
    prod.replacefds1 = {{'vect1', 'vect'}};
    prod.replacefds2 = {{'vect', 'vect2'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'vect2'};
    prod = prod.setup;

    P = prod.eval(invNtN, sN); % P is in cellintfacevectbl

    % Setup matrix for vector flux reconstruction, from facetbl to cellvecttbl.
    
    prod = TensorProd();
    prod.tbl1 = cellintfacevecttbl;
    prod.tbl2 = intfacetbl;
    prod.tbl3 = cellvecttbl;
    prod.reducefds = {'faces'};

    P_T = SparseTensor();
    P_T = P_T.setFromTensorProd(P, prod);

    P = P_T.getMatrix();

    % Setup matrix that compute the sum over the dimension
    map = TensorMap();
    map.fromTbl = cellvecttbl;
    map.toTbl = celltbl;
    map.mergefds = {'cells'};

    S_T = SparseTensor();
    S_T = S_T.setFromTensorMap(map);

    S = S_T.getMatrix();

    op.P = P;
    op.S = S;
    
end

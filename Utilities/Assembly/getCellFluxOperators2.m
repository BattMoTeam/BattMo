function op = getCellFluxOperators2(G)
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
    
    doOptimized = false;
    
    tbls = setupSimpleTables(G);

    celltbl = tbls.celltbl;
    facetbl = tbls.facetbl;
    cellfacetbl = tbls.cellfacetbl;
    
    if doOptimized
        map = TensorMap();
        map.fromTbl = celltbl;
        map.toTbl = cellfacetbl;
        map.mergefds = {'cells'};
        cell_from_cellface = map.getDispatchInd();
        map.fromTbl = facetbl;
        map.toTbl = cellfacetbl;
        map.mergefds = {'faces'};
        face_from_cellface = map.getDispatchInd();
    end

    vecttbl.vect = (1 : G.griddim)';
    vecttbl = IndexArray(vecttbl);

    cn = cellfacetbl.get('cells');
    cf = cellfacetbl.get('faces');
    sgn = 2*double(G.faces.neighbors(cf, 1) == cn) - 1;

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = vecttbl;
    gen.tbl2 = vecttbl;
    gen.replacefds1 = {{'vect', 'vect1'}};
    gen.replacefds2 = {{'vect', 'vect2'}};
    gen.mergefds = {};

    vect12tbl = gen.eval();

    if doOptimized
        opts = {'optpureproduct', true, 'virtual', true};
    else
        opts = {'optpureproduct', true};
    end
    facevecttbl       = crossIndexArray(facetbl    , vecttbl  , {}, opts{:});
    facevecttbl       = crossIndexArray(facetbl       , vecttbl  , {}, opts{:});
    cellvecttbl       = crossIndexArray(celltbl       , vecttbl  , {}, opts{:});
    cellvect12tbl     = crossIndexArray(celltbl       , vect12tbl, {}, opts{:});
    cellfacevecttbl   = crossIndexArray(cellfacetbl, vecttbl  , {}, opts{:});
    cellfacevect12tbl = crossIndexArray(cellfacetbl, vect12tbl, {}, opts{:});
    cellfacevecttbl   = crossIndexArray(cellfacetbl   , vecttbl  , {}, opts{:});

    if doOptimized
        % some shorthands
        d_num   = vecttbl.num;
        if_num  = intfacetbl.num;
        icf_num = cellintfacetbl.num;
        c_num   = celltbl.num;
        f_num   = facetbl.num;
        cf_num  = cellfacetbl.num;
        
        face_from_intface = intfacetbl.get('faces');
        
        map = TensorMap();
        map.fromTbl = cellfacetbl;
        map.toTbl = cellintfacetbl;
        map.mergefds = {'cells', 'faces'};
        cellface_from_cellintface = map.getDispatchInd();
        
        intface_from_face = zeros(facetbl.num, 1);
        intface_from_face(face_from_intface) = (1 : intfacetbl.num)';
        
        cellintface_from_cellface = zeros(cellfacetbl.num, 1);
        cellintface_from_cellface(cellface_from_cellintface) = (1 : cellintfacetbl.num)';        
    end
    
    N = G.faces.normals;
    N = reshape(N', [], 1); % N is in facevecttbl

    prod = TensorProd();
    prod.tbl1 = cellfacetbl;
    prod.tbl2 = facevecttbl;
    prod.tbl3 = cellfacevecttbl;
    prod.mergefds = {'faces'};
    if doOptimized
        prod.pivottbl = cellintfacevecttbl;
        [r, i] = ind2sub([d_num, icf_num], (1 : cellintfacevecttbl.num)');
        prod.dispind1 = i;
        prod.dispind2 = sub2ind([d_num, f_num], r, intface_from_face(face_from_cellface(cellface_from_cellintface(i))));
        prod.dispind3 = (1 : cellintfacevecttbl.num)';
        prod.issetup = true;
    else
        prod = prod.setup();
    end

    N = prod.eval(sgn, N); % N is in cellfacevecttbl

    % We compute NtN

    prod = TensorProd();
    prod.tbl1 = cellfacevecttbl;
    prod.tbl2 = cellfacevecttbl;
    prod.tbl3 = cellvect12tbl;
    prod.replacefds1 = {{'vect', 'vect1'}};
    prod.replacefds2 = {{'vect', 'vect2'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'faces'};
    
    if doOptimized
        prod.pivottbl = cellintfacevect12tbl;
        [r2, r1, i] = ind2sub([d_num, d_num, icf_num], (1 : cellintfacevect12tbl.num)');
        prod.dispind1 = sub2ind([d_num, icf_num], r1, i);
        prod.dispind2 = sub2ind([d_num, icf_num], r2, i);
        prod.dispind3 = sub2ind([d_num, d_num, c_num], r1, r2, cell_from_cellface(cellface_from_cellintface(i)));
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
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
    if doOptimized
        prod.pivottbl = cellvect12tbl;
        [r2, r1, i] = ind2sub([d_num, d_num, c_num], (1 : cellvect12tbl.num)');
        prod.dispind1 = (1 : cellvect12tbl.num)';
        prod.dispind2 = sub2ind([d_num, c_num], r2, i);
        prod.dispind3 = sub2ind([d_num, c_num], r1, i);
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    A = SparseTensor();
    A = A.setFromTensorProd(NtN, prod);
    A = A.getMatrix();

    % fetch block diagonal inverter
    opt.invertBlocks = 'matlab';
    bi = blockInverter(opt);
    sz = vecttbl.num*ones(celltbl.num, 1);
    invNtN = bi(A, sz);

    if doOptimized
        [r2, r1, i] = ind2sub([d_num, d_num, c_num], (1 : cellvect12tbl.num)');
        ind1 = sub2ind([d_num, c_num], r1, i);
        ind2 = sub2ind([d_num, c_num], r2, i);
    else
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
    end

    ind = sub2ind([cellvecttbl.num, cellvecttbl.num], ind1, ind2);

    invNtN = invNtN(ind); % invNtN is in cellvect12tbl

    % We set up mapping from face fluxes to signed fluxes

    prod = TensorProd();
    prod.tbl1 = cellfacevecttbl;
    prod.tbl2 = cellfacetbl;
    prod.tbl3 = cellfacevecttbl;
    prod.mergefds = {'cells', 'faces'};
    if doOptimized
        prod.pivottbl = cellintfacevecttbl;
        [r, i] = ind2sub([d_num, c_num], (1 : cellintfacevecttbl.num)');
        prod.dispind1 = (1 : cellintfacevecttbl.num)';
        prod.dispind2 = i;
        prod.dispind3 = (1 : cellintfacevecttbl.num)';
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    sN = prod.eval(N, sgn); % sN is in cellfacevecttbl

    % We compute P = (invNtN)*sN

    prod = TensorProd();
    prod.tbl1 = cellvect12tbl;
    prod.tbl2 = cellfacevecttbl;
    prod.tbl3 = cellfacevecttbl;
    prod.replacefds1 = {{'vect1', 'vect'}};
    prod.replacefds2 = {{'vect', 'vect2'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'vect2'};
    
    if doOptimized
        prod.pivottbl = cellintfacevect12tbl;
        [r2, r1, i] = ind2sub([d_num, d_num, icf_num], (1 : cellintfacevect12tbl.num)');
        prod.dispind1 = sub2ind([d_num, d_num, c_num], r2, r1, cell_from_cellface(cellface_from_cellintface(i)));
        prod.dispind2 = sub2ind([d_num, icf_num], r2, i);
        prod.dispind3 = sub2ind([d_num, icf_num], r1, i);
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    P = prod.eval(invNtN, sN); % P is in cellfacevectbl

    % Setup matrix for vector flux reconstruction, from facetbl to cellvecttbl.
    
    prod = TensorProd();
    prod.tbl1 = cellfacevecttbl;
    prod.tbl2 = facetbl;
    prod.tbl3 = cellvecttbl;
    prod.reducefds = {'faces'};
    if doOptimized
        prod.pivottbl = cellintfacevecttbl;
        [r, i] = ind2sub([d_num, icf_num], (1 : cellintfacevecttbl.num)');
        prod.dispind1 = (1 : cellintfacevecttbl.num)';
        prod.dispind2 = intface_from_face(face_from_cellface(cellface_from_cellintface(i)));
        prod.dispind3 = sub2ind([d_num, c_num], r, cell_from_cellface(cellface_from_cellintface(i)));
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    P_T = SparseTensor();
    P_T = P_T.setFromTensorProd(P, prod);

    P = P_T.getMatrix();

    % Setup matrix that compute the sum over the dimension
    map = TensorMap();
    map.fromTbl = cellvecttbl;
    map.toTbl = celltbl;
    map.mergefds = {'cells'};
    if doOptimized
        map.pivottbl = cellvecttbl;
        [r, i] = ind2sub([d_num, c_num], (1 : cellvecttbl.num)');
        map.dispind1 = (1 : cellvecttbl.num)';
        map.dispind2 = i;
        map.issetup = true;
    else
        map = map.setup();
    end
    
    S_T = SparseTensor();
    S_T = S_T.setFromTensorMap(map);

    S = S_T.getMatrix();

    op.P = P;
    op.S = S;
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}

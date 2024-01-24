function operators = localSetupOperators(G, varargin)
    opts = struct('assembleCellFluxOperator', false);
    opts = merge_options(opts, varargin{:});

    nc = G.cells.num;
    cst = ones(nc, 1);
    rock = struct('perm', cst, 'poro', cst);
    operators = setupOperatorsTPFA(G, rock);

    hT = computeTrans(G, rock);

    tbls = setupSimpleTables(G);
    cellfacetbl = tbls.cellfacetbl;
    celltbl     = tbls.celltbl;
    facetbl     = tbls.facetbl;
    intfacetbl  = tbls.intfacetbl;

    map = TensorMap();
    map.fromTbl = celltbl;
    map.toTbl = cellfacetbl;
    map.mergefds = {'cells'};
    map = map.setup();

    M = SparseTensor();
    M = M.setFromTensorMap(map);
    M = M.getMatrix;

    map = TensorMap();
    map.fromTbl = cellfacetbl;
    map.toTbl = intfacetbl;
    map.mergefds = {'faces'};
    map = map.setup();

    P = SparseTensor();
    P = P.setFromTensorMap(map);
    P = P.getMatrix;

    T = operators.T;

    operators.harmFace    = @(cellvalue) getHarmFace(cellvalue, P, hT, T, M);
    operators.harmFaceBC  = @(cvalue, faces) getFaceHarmBC(G, cvalue, faces);
    operators.transFaceBC = @(faces) getTransFaceBC(G, faces);

    %% setup the sign for *external* faces
    cells  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    faces  = G.cells.faces(:, 1);

    extfaces = find(any(G.faces.neighbors == 0, 2));
    extcells = sum(G.faces.neighbors(extfaces, :), 2);

    extsgn = 2*(extcells == G.faces.neighbors(extfaces, 1)) - 1;

    sgn = nan(G.faces.num, 1);
    sgn(extfaces) = extsgn;

    operators.sgn = sgn;

    %% setup cell flux reconstruction operator
    if opts.assembleCellFluxOperator
        %operators.cellFluxOp = getCellFluxOperators2(G);
        operators.cellFluxOp = getCellFluxOperatorsAll(G);
    end
end

function facevalue = getHarmFace(cellvalue, P, hT, T, M)

    if numelValue(cellvalue) > 1
        facevalue = 1./(P*(1./(hT.*(M*cellvalue))));
    else
        facevalue = cellvalue.*T;
    end

end
function [T, cells] = getFaceHarmBC(G, cvalue, faces)
    [t, cells] = getTransFaceBC(G, faces);
    T = t.*cvalue(cells);
end

function [t, cells] = getTransFaceBC(G, faces)
    cells = sum(G.faces.neighbors(faces, :), 2);
    cn = sqrt(sum((G.faces.centroids(faces, :) - G.cells.centroids(cells, :)).^2, 2));
    t = G.faces.areas(faces)./cn;
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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

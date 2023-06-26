function operators = localSetupOperators(G, varargin)
    opts = struct('assembleCellFluxOperator', false);
    opts = merge_options(opts, varargin{:});

    nc = G.cells.num;
    cst = ones(nc, 1);
    rock = struct('perm', cst, 'poro', cst);
    operators = setupOperatorsTPFA(G, rock);

    hT = computeTrans(G, rock);

    if isa(G.cells.c_length_factors, 'cell')

        weight = G.cells.c_length_factors{1} + 0 * G.cells.c_length_factors{2};
        operators.pv = operators.pv .* weight;
        G.cells.volumes = G.cells.volumes .* weight;
        hT = hT ./ weight;

        assert(size(operators.pv.jac{1}, 2) == 2);
        assert(size(G.cells.volumes.jac{1}, 2) == 2);
        assert(size(hT.jac{1}, 2) == 2);
        
    elseif isa(G.cells.c_length_factors, 'struct')

        % Original
        pv0 = operators.pv;
        vol0 = G.cells.volumes;
        hT0 = hT;

        % Initialize
        ad = G.cells.c_length_factors.NegativeElectrode.length_factor + G.cells.c_length_factors.PositiveElectrode.length_factor;
        operators.pv = pv0 + 0 * ad;
        G.cells.volumes = vol0 + 0 * ad;
        hT = hT0 + 0 * ad;
        
        % Weight
        eldes = {'NegativeElectrode', 'PositiveElectrode'};
        for k = 1:2
            elde = eldes{k};
            k2 = rem(k, 2) + 1;
            other_elde = eldes{k2};
            
            cells  = G.cells.c_length_factors.(elde).cells;
            cf     = G.cells.cf_length_factors.(elde).cf;
            weight = G.cells.c_length_factors.(elde).length_factor + 0 * G.cells.c_length_factors.(other_elde).length_factor;

            operators.pv(cells) = pv0(cells) * weight;
            G.cells.volumes(cells) = vol0(cells) * weight;
            hT(cf) = hT0(cf) ./ weight;
        end

        assert(size(operators.pv.jac{1}, 2) == 2);
        assert(size(G.cells.volumes.jac{1}, 2) == 2);
        assert(size(hT.jac{1}, 2) == 2);

    else
        keyboard;
    end

    tbls = setupSimpleTables(G);
    cellfacetbl = tbls.cellfacetbl;
    celltbl     = tbls.celltbl;
    %facetbl     = tbls.facetbl;
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

    operators.harmFace = @(cellvalue) 1./(P*(1./(hT.*(M*cellvalue))));

    operators.harmFaceBC = @(cvalue, faces) getFaceHarmBC(G, cvalue, faces);

    %% setup the sign for *external* faces
    % cells  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    % faces  = G.cells.faces(:, 1);

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

function [T, cells] = getFaceHarmBC(G, cvalue, faces)
    cells = sum(G.faces.neighbors(faces, :), 2);
    cn = sqrt(sum((G.faces.centroids(faces, :) - G.cells.centroids(cells, :)).^2, 2));
    t = G.faces.areas(faces)./cn;

    if isa(G.cells.c_length_factors, 'cell')

        weight = G.cells.c_length_factors{1} + 0 * G.cells.c_length_factors{2};
        t = t ./ weight;

        assert(size(t.jac{1}, 2) == 2);
        
    elseif isa(G.cells.c_length_factors, 'struct')

        ad = G.cells.c_length_factors{1} + G.cells.c_length_factors{2};
        t0 = t + 0 * ad;
        
        eldes = {'NegativeElectrode', 'PositiveElectrode'};
        for k = 1:2
            elde = eldes{k};
            other_elde = eldes{k2};
            cf     = G.cells.cf_length_factors.(elde).cf;
            weight = G.cells.c_length_factors.(elde).length_factor + 0 * G.cells.c_length_factors.(other_elde).length_factor;
            t(cf) = t0(cf) ./ weight;
        end
            
        assert(size(t.jac{1}, 2) == 2);

    else
        keyboard;
    end
        
    T = t.*cvalue(cells);
end


%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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

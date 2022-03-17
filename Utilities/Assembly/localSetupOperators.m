function operators = localSetupOperators(G, varargin)
    opts = struct('assembleCellFluxOperator', false);
    opts = merge_options(opts, varargin{:});
    
    nc = G.cells.num;
    cst = ones(nc, 1);
    rock = struct('perm', cst, 'poro', cst);
    operators = setupOperatorsTPFA(G, rock);
    operators.allDiv = getAllDiv(G);
    %operators.harmFace = getFaceHarmMean(G);
    operators.harmFace =@(cellvalue) getTrans(G).*(1./operators.faceAvg(1./cellvalue));
    operators.harmFaceBC = @(cvalue, faces) getFaceHarmBC(G, cvalue, faces);
    
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

function [T, cells] = getFaceHarmBC(G, cvalue, faces)
    cells = sum(G.faces.neighbors(faces, :), 2);
    cn = sqrt(sum((G.faces.centroids(faces, :) - G.cells.centroids(cells, :)).^2, 2));
    t = G.faces.areas(faces)./cn;
    T = t.*cvalue(cells);
end

function hm = getFaceHarmMean(G)
    internal = all(G.faces.neighbors>0, 2);
    N = G.faces.neighbors(internal, :);
    ni = sum(internal);
    cd = sqrt(sum((G.cells.centroids(N(:, 1), :) - G.cells.centroids(N(:, 2), :)).^2, 2)); % NB
    t = G.faces.areas(internal)./cd;
    A = sparse([[1:ni]'; [1:ni]'], N, 1, ni, G.cells.num);
    hm = @(cellvalue) 2.*t./(A*(1./cellvalue));
end

function t = getTrans(G)
    internal = all(G.faces.neighbors>0, 2);
    N = G.faces.neighbors(internal, :);
    cd = sqrt(sum((G.cells.centroids(N(:, 1), :) - G.cells.centroids(N(:, 2), :)).^2, 2)); % NB
    t = G.faces.areas(internal)./cd;
end


function tp = getTwoPointOperator(G)
% Mappings from cells to its faces
    cells = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    faces = G.cells.faces(:, 1);
    % Vector from cell to face centroid
    C = G.faces.centroids(faces, :) - G.cells.centroids(cells, :);
    % Oriented normals
    sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1;
    N = bsxfun(@times, sgn, G.faces.normals(faces, :));
    % Make function
    cn = sum(C.*N, 2)./sum(C.*C, 2);
    tp = @(lambda) cn.*lambda(cells);
end

function ha = getHarmonicAvgOperator(G)
% Harmonic averaging operator
    faces = G.cells.faces(:, 1);
    M = sparse(faces, 1:numel(faces), 1, G.faces.num, numel(faces));
    ha = @(T) 1./(M*(1./T));
end

function allDiv = getAllDiv(G)
    nc = G.cells.num;
    nf = G.faces.num;
    Nall = G.faces.neighbors;
    internal = all(Nall>0, 2);
    ifn = find(internal);
    efn = find(~internal);
    inf = numel(ifn);
    N = Nall(internal, :);
    Nb = Nall(~internal, :);
    signb = 2*(Nb(:, 1)>0) - 1;
    Nb = sum(Nb, 2);
    C = sparse([ifn; ifn], N, ones(inf, 1) * [1,- 1], nf, nc);
    C = C + sparse(efn, Nb, signb, nf, nc);
    allDiv = @(x) (C'*x)./G.cells.volumes;
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

function tbls = setupSimpleTables(G)
    
    nc = G.cells.num;
    nf = G.faces.num;
    nn = G.nodes.num;
    
    celltbl.cells = (1 : nc)';
    celltbl = IndexArray(celltbl);

    facetbl.faces = (1 : nf)';
    facetbl = IndexArray(facetbl);

    nodetbl.nodes = (1 : nn)';
    nodetbl = IndexArray(nodetbl);
    
    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl = IndexArray(cellfacetbl);
    
    facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos));
    facenodetbl.nodes = G.faces.nodes(:, 1);
    facenodetbl = IndexArray(facenodetbl);
    
    tbls = struct('celltbl'    , celltbl    , ...
                  'facetbl'    , facetbl    , ...
                  'nodetbl'    , nodetbl    , ...
                  'cellfacetbl', cellfacetbl, ...
                  'facenodetbl', facenodetbl);
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

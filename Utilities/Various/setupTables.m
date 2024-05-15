function tbls = setupTables(G, varargin)

    opt.includetbls = {};
    opt = merge_options(opt, varargin{:});
    includetbls = opt.includetbls;

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


    if ismember('cellnodetbl', includetbls)

        cellfacenodetbl = crossIndexArray(cellfacetbl, facenodetbl, {'faces'});
        cellnodetbl = projIndexArray(cellfacenodetbl, {'cells', 'nodes'});

        tbls.cellnodetbl = cellnodetbl;
        
    end
    
    if any(ismember({'intfacetbl', 'cellintfacetbl'}, includetbls))

        intfaces = all(G.faces.neighbors> 0, 2);
        intfacetbl.faces = find(intfaces);
        intfacetbl = IndexArray(intfacetbl);

        cellintfacetbl = crossIndexArray(intfacetbl, cellfacetbl, {'faces'});

        tbls.intfacetbl     = intfacetbl;
        tbls.cellintfacetbl = cellintfacetbl;

    end

    if ismember('extfacetbl', includetbls)

        extfaces = any(G.faces.neighbors == 0, 2);
        extfacetbl.faces = find(extfaces);
        extfacetbl = IndexArray(extfacetbl);

        tbls.extfacetbl = extfacetbl;

    end

    if ismember('vectbl', includetbls)

        vectbl.vec = (1 : G.griddim)';
        vectbl = IndexArray(vectbl);

        facevectbl           = crossIndexArray(facetbl          , vectbl       , {}, 'optpureproduct', true);
        nodevectbl           = crossIndexArray(nodetbl          , vectbl       , {}, 'optpureproduct', true);
        cellvectbl           = crossIndexArray(celltbl          , vectbl       , {}, 'optpureproduct', true);
        cellfacevectbl       = crossIndexArray(cellfacetbl      , vectbl       , {}, 'optpureproduct', true);
        facenodevectbl       = crossIndexArray(facenodetbl      , vectbl       , {}, 'optpureproduct', true);

        tbls.vectbl         = vectbl;
        tbls.facevectbl     = facevectbl;
        tbls.nodevectbl     = nodevectbl;
        tbls.cellvectbl     = cellvectbl;
        tbls.cellfacevectbl = cellfacevectbl;
        tbls.facenodevectbl = facenodevectbl;

    end

    if ismember('facenode12tbl', includetbls)

        p = G.faces.nodePos;
        f = G.faces.nodes;

        next = (2 : size(f, 1) + 1) .';
        next(p(2 : end) - 1) = p(1 : end-1);
        nextf = f(next);

        facenode12tbl = replacefield(facenodetbl, {{'nodes', 'nodes1'}});
        facenode12tbl = facenode12tbl.addInd('nodes2', nextf);

        tbls.facenode12tbl = facenode12tbl;

    end

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

function G = genSubGrid(globG, cells)
    nc = globG.cells.num;
    rcells = (1 : nc)';
    rcells(cells) = [];
    [G, cellmap, facemap, nodemap] = removeCells(globG, rcells);
    G = computeGeometry(G);
    mappings = struct('parentGrid', globG     , ...
                      'cellmap'   , cellmap, ...
                      'facemap'   , facemap, ...
                      'nodemap'   , nodemap);
    G.mappings = mappings;
end

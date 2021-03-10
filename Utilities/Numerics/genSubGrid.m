function G = genSubGrid(G0, cells)
    nc = G0.cells.num;
    rcells = (1 : nc)';
    rcells(cells) = [];
    [G, cellmap, facemap, nodemap] = removeCells(G0, rcells);
    G = computeGeometry(G);
    mappings = struct('parentGrid', G0     , ...
                      'cellmap'   , cellmap, ...
                      'facemap'   , facemap, ...
                      'nodemap'   , nodemap);
    G.mappings = mappings;
end

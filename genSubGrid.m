function G = genSubGrid(G0, cells)
    nc = 1 : G0.cells.num;
    rcells = (1 : nc)';
    rcells(cells) = [];
    [G, cellmap, facemap, nodemap] = removeCells(G0, cells);
    mappings = struct('parentGrid', G0     , ...
                      'cellmap'   , cellmap, ...
                      'facemap'   , facemap, ...
                      'nodemap'   , nodemap);
    G.mappings = mappings;
end

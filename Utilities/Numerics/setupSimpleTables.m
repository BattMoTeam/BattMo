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
    
    tbls = struct('celltbl'         , celltbl         , ...
                  'facetbl'         , facetbl         , ...
                  'nodetbl'         , nodetbl         , ...
                  'cellfacetbl'     , cellfacetbl);
end

function coupTerm = setupCcneNeCoupTerm(model)

    ne = model.ne;
    ccne = model.ccne;

    Gne = ne.G;
    Gccne = ccne.G;

    G = Gne.mappings.parentGrid;

    netbls = setupSimpleTables(Gne);
    ccnetbls = setupSimpleTables(Gccne);
    tbls = setupSimpleTables(G);            
    
    necelltbl = netbls.celltbl;
    necelltbl = necelltbl.addInd('globcells', Gne.mappings.cellmap);
    nefacetbl = netbls.facetbl;
    nefacetbl = nefacetbl.addInd('globfaces', Gne.mappings.facemap);

    necellfacetbl = netbls.cellfacetbl;
    necellfacetbl = crossIndexArray(necellfacetbl, necelltbl, {'cells'});
    necellfacetbl = crossIndexArray(necellfacetbl, nefacetbl, {'faces'});
    
    ccnecelltbl = ccnetbls.celltbl;
    ccnecelltbl = ccnecelltbl.addInd('globcells', Gccne.mappings.cellmap);
    ccnefacetbl = ccnetbls.facetbl;
    ccnefacetbl = ccnefacetbl.addInd('globfaces', Gccne.mappings.facemap);
    
    ccnecellfacetbl = ccnetbls.cellfacetbl;
    ccnecellfacetbl = crossIndexArray(ccnecellfacetbl, ccnecelltbl, {'cells'});
    ccnecellfacetbl = crossIndexArray(ccnecellfacetbl, ccnefacetbl, {'faces'});
    
    gen = CrossIndexArrayGenerator();
    gen.tbl1 = necellfacetbl;
    gen.tbl2 = ccnecellfacetbl;
    gen.replacefds1 = {{'cells', 'necells'}, {'faces', 'nefaces'}, {'globcells', 'neglobcells'}};
    gen.replacefds2 = {{'cells', 'ccnecells'}, {'faces', 'ccnefaces'}, {'globcells', 'ccneglobcells'}};
    gen.mergefds = {'globfaces'};
    
    cell12facetbl = gen.eval();

    ccnefaces = cell12facetbl.get('ccnefaces');
    nefaces = cell12facetbl.get('nefaces');
    ccnecells = cell12facetbl.get('ccnecells');
    necells = cell12facetbl.get('necells');            
    
    compnames = {'ccne', 'ne'};
    coupTerm = couplingTerm('ccne-ne', compnames);
    coupTerm.couplingfaces =  [ccnefaces, nefaces];
    coupTerm.couplingcells = [ccnecells, necells];

end
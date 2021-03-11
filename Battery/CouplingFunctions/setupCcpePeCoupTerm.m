function coupTerm = setupCcpePeCoupTerm(model)

    pe = model.pe;
    ccpe = model.ccpe;

    Gpe = pe.G;
    Gccpe = ccpe.G;

    G = Gpe.mappings.parentGrid;

    petbls = setupSimpleTables(Gpe);
    ccpetbls = setupSimpleTables(Gccpe);
    tbls = setupSimpleTables(G);            
    
    pecelltbl = petbls.celltbl;
    pecelltbl = pecelltbl.addInd('globcells', Gpe.mappings.cellmap);
    pefacetbl = petbls.facetbl;
    pefacetbl = pefacetbl.addInd('globfaces', Gpe.mappings.facemap);

    pecellfacetbl = petbls.cellfacetbl;
    pecellfacetbl = crossIndexArray(pecellfacetbl, pecelltbl, {'cells'});
    pecellfacetbl = crossIndexArray(pecellfacetbl, pefacetbl, {'faces'});
    
    ccpecelltbl = ccpetbls.celltbl;
    ccpecelltbl = ccpecelltbl.addInd('globcells', Gccpe.mappings.cellmap);
    ccpefacetbl = ccpetbls.facetbl;
    ccpefacetbl = ccpefacetbl.addInd('globfaces', Gccpe.mappings.facemap);
    
    ccpecellfacetbl = ccpetbls.cellfacetbl;
    ccpecellfacetbl = crossIndexArray(ccpecellfacetbl, ccpecelltbl, {'cells'});
    ccpecellfacetbl = crossIndexArray(ccpecellfacetbl, ccpefacetbl, {'faces'});
    
    gen = CrossIndexArrayGenerator();
    gen.tbl1 = pecellfacetbl;
    gen.tbl2 = ccpecellfacetbl;
    gen.replacefds1 = {{'cells', 'pecells'}, {'faces', 'pefaces'}, {'globcells', 'peglobcells'}};
    gen.replacefds2 = {{'cells', 'ccpecells'}, {'faces', 'ccpefaces'}, {'globcells', 'ccpeglobcells'}};
    gen.mergefds = {'globfaces'};
    
    cell12facetbl = gen.eval();

    ccpefaces = cell12facetbl.get('ccpefaces');
    pefaces = cell12facetbl.get('pefaces');
    ccpecells = cell12facetbl.get('ccpecells');
    pecells = cell12facetbl.get('pecells');            
    
    compnames = {'ccpe', 'pe'};
    coupTerm = couplingTerm('ccpe-pe', compnames);
    coupTerm.couplingfaces =  [ccpefaces, pefaces];
    coupTerm.couplingcells = [ccpecells, pecells];
    
end

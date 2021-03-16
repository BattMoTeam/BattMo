classdef BatteryInputParams2D < BatteryInputParams
    
    methods
        
        function params = setupVariousParams(params)
            
            params.SOC  = 0.5;
            params.T    = 298.15;
            params.J    = 0.1;
            params.Ucut = 2;
        
        end

        function params = setupSubModels(params)
            
            fac = 1;
            
            sepnx  = 30*fac;
            nenx   = 30*fac;
            penx   = 30*fac;
            ccnenx = 20*fac;
            ccpenx = 20*fac;

            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            ny = 10*fac;

            xlength = 1e-6*[10; 100; 50; 80; 10];
            ylength = 1e-2;

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = ylength/ny;
            y = rldecode(y, ny);
            y = [0; cumsum(y)];

            G = tensorGrid(x, y);
            G = computeGeometry(G);
            params.G = G;

            %% setup elyte
            nx = sum(nxs);

            istart = ccnenx + 1;
            ni = nenx + sepnx + penx;
            cells = pickTensorCells(istart, ni, nx, ny);
            params.elyte = orgLiPF6('elyte', G, cells);

            %% setup ne
            istart = ccnenx + 1;
            cells = pickTensorCells(istart, nenx, nx, ny);
            params.ne = graphiteElectrode('ne', G, cells);

            %% setup pe
            istart = ccnenx + nenx + sepnx + 1;
            cells = pickTensorCells(istart, penx, nx, ny);
            params.pe = nmc111Electrode('pe', G, cells);

            %% setup ccne
            istart = 1;
            cells = pickTensorCells(istart, ccnenx, nx, ny);
            params.ccne = currentCollector('ccne', G, cells);

            %% setup ccpe
            istart = ccnenx + nenx + sepnx + penx + 1;
            cells = pickTensorCells(istart, ccpenx, nx, ny);
            params.ccpe = currentCollector('ccpe', G, cells);

            %% setup sep
            istart = ccnenx + nenx + 1;
            cells = pickTensorCells(istart, sepnx, nx, ny);
            params.sep = celgard2500('sep', G, cells);
        
        end
        
        function coupTerm = setupCcneBcCoupTerm(params)

            ccne = params.ccne;
            G = ccne.G;

            % We pick up the faces at the top of Cccne
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf - myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccne'};
            coupTerm = couplingTerm('bc-ccne', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupCcneNeCoupTerm(params)

            ne = params.ne;
            ccne = params.ccne;

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
        
        function coupTerm = setupCcpeBcCoupTerm(params)

            ccpe = params.ccpe;
            G = ccpe.G;

            % We pick up the faces at the top of Cccpe
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf - myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccpe'};
            coupTerm = couplingTerm('bc-ccpe', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupCcpePeCoupTerm(params)

            pe = params.pe;
            ccpe = params.ccpe;

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

        function coupTerm = setupNeElyteCoupTerm(params)
            
            ne = params.ne;
            elyte = params.elyte;
            
            Gne = ne.G;
            Gelyte = elyte.G;
            
            % parent Grid
            G = Gne.mappings.parentGrid;
            
            % All the cells from ne are coupled with elyte
            cells1 = (1 : Gne.cells.num)';
            pcells = Gne.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(Gelyte.mappings.cellmap) = (1 : Gelyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'ne', 'elyte'};
            coupTerm = couplingTerm('ne-elyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
        end
        
        function coupTerm = setupPeElyteCoupTerm(params)
            
            pe = params.pe;
            elyte = params.elyte;
            
            Gpe = pe.G;
            Gelyte = elyte.G;
            
            % parent Grid
            G = Gpe.mappings.parentGrid;
            
            % All the cells from pe are coupled with elyte
            cells1 = (1 : Gpe.cells.num)';
            pcells = Gpe.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(Gelyte.mappings.cellmap) = (1 : Gelyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'pe', 'elyte'};
            coupTerm = couplingTerm('pe-elyte', compnames);
            coupTerm.couplingcells = [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling between faces
            
        end

    end
    
end

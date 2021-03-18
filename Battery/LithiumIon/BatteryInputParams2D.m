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
            NegativeElectrodenx   = 30*fac;
            PositiveElectrodenx   = 30*fac;
            NegativeCurrentCollectornx = 20*fac;
            PositiveCurrentCollectornx = 20*fac;

            nxs = [NegativeCurrentCollectornx; NegativeElectrodenx; sepnx; PositiveElectrodenx; PositiveCurrentCollectornx];
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

            %% setup Electrolyte
            nx = sum(nxs);

            istart = NegativeCurrentCollectornx + 1;
            ni = NegativeElectrodenx + sepnx + PositiveElectrodenx;
            cells = pickTensorCells(istart, ni, nx, ny);
            params.Electrolyte = orgLiPF6('Electrolyte', G, cells);

            %% setup NegativeElectrode
            istart = NegativeCurrentCollectornx + 1;
            cells = pickTensorCells(istart, NegativeElectrodenx, nx, ny);
            params.NegativeElectrode = GraphiteElectrode('NegativeElectrode', G, cells);

            %% setup PositiveElectrode
            istart = NegativeCurrentCollectornx + NegativeElectrodenx + sepnx + 1;
            cells = pickTensorCells(istart, PositiveElectrodenx, nx, ny);
            params.PositiveElectrode = NMC111Electrode('PositiveElectrode', G, cells);

            %% setup NegativeCurrentCollector
            istart = 1;
            cells = pickTensorCells(istart, NegativeCurrentCollectornx, nx, ny);
            params.NegativeCurrentCollector = CurrentCollector('NegativeCurrentCollector', G, cells);

            %% setup PositiveCurrentCollector
            istart = NegativeCurrentCollectornx + NegativeElectrodenx + sepnx + PositiveElectrodenx + 1;
            cells = pickTensorCells(istart, PositiveCurrentCollectornx, nx, ny);
            params.PositiveCurrentCollector = CurrentCollector('PositiveCurrentCollector', G, cells);

            %% setup sep
            istart = NegativeCurrentCollectornx + NegativeElectrodenx + 1;
            cells = pickTensorCells(istart, sepnx, nx, ny);
            params.sep = celgard2500('sep', G, cells);
        
        end
        
        function coupTerm = setupNegativeCurrentCollectorBcCoupTerm(params)

            NegativeCurrentCollector = params.NegativeCurrentCollector;
            G = NegativeCurrentCollector.G;

            % We pick up the faces at the top of CNegativeCurrentCollector
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf - myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'NegativeCurrentCollector'};
            coupTerm = couplingTerm('bc-NegativeCurrentCollector', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupNegativeCurrentCollectorNegativeElectrodeCoupTerm(params)

            NegativeElectrode = params.NegativeElectrode;
            NegativeCurrentCollector = params.NegativeCurrentCollector;

            GNegativeElectrode = NegativeElectrode.G;
            GNegativeCurrentCollector = NegativeCurrentCollector.G;

            G = GNegativeElectrode.mappings.parentGrid;

            NegativeElectrodetbls = setupSimpleTables(GNegativeElectrode);
            NegativeCurrentCollectortbls = setupSimpleTables(GNegativeCurrentCollector);
            tbls = setupSimpleTables(G);            
            
            NegativeElectrodecelltbl = NegativeElectrodetbls.celltbl;
            NegativeElectrodecelltbl = NegativeElectrodecelltbl.addInd('globcells', GNegativeElectrode.mappings.cellmap);
            NegativeElectrodefacetbl = NegativeElectrodetbls.facetbl;
            NegativeElectrodefacetbl = NegativeElectrodefacetbl.addInd('globfaces', GNegativeElectrode.mappings.facemap);

            NegativeElectrodecellfacetbl = NegativeElectrodetbls.cellfacetbl;
            NegativeElectrodecellfacetbl = crossIndexArray(NegativeElectrodecellfacetbl, NegativeElectrodecelltbl, {'cells'});
            NegativeElectrodecellfacetbl = crossIndexArray(NegativeElectrodecellfacetbl, NegativeElectrodefacetbl, {'faces'});
            
            NegativeCurrentCollectorcelltbl = NegativeCurrentCollectortbls.celltbl;
            NegativeCurrentCollectorcelltbl = NegativeCurrentCollectorcelltbl.addInd('globcells', GNegativeCurrentCollector.mappings.cellmap);
            NegativeCurrentCollectorfacetbl = NegativeCurrentCollectortbls.facetbl;
            NegativeCurrentCollectorfacetbl = NegativeCurrentCollectorfacetbl.addInd('globfaces', GNegativeCurrentCollector.mappings.facemap);
            
            NegativeCurrentCollectorcellfacetbl = NegativeCurrentCollectortbls.cellfacetbl;
            NegativeCurrentCollectorcellfacetbl = crossIndexArray(NegativeCurrentCollectorcellfacetbl, NegativeCurrentCollectorcelltbl, {'cells'});
            NegativeCurrentCollectorcellfacetbl = crossIndexArray(NegativeCurrentCollectorcellfacetbl, NegativeCurrentCollectorfacetbl, {'faces'});
            
            gen = CrossIndexArrayGenerator();
            gen.tbl1 = NegativeElectrodecellfacetbl;
            gen.tbl2 = NegativeCurrentCollectorcellfacetbl;
            gen.replacefds1 = {{'cells', 'NegativeElectrodecells'}, {'faces', 'NegativeElectrodefaces'}, {'globcells', 'NegativeElectrodeglobcells'}};
            gen.replacefds2 = {{'cells', 'NegativeCurrentCollectorcells'}, {'faces', 'NegativeCurrentCollectorfaces'}, {'globcells', 'NegativeCurrentCollectorglobcells'}};
            gen.mergefds = {'globfaces'};
            
            cell12facetbl = gen.eval();

            NegativeCurrentCollectorfaces = cell12facetbl.get('NegativeCurrentCollectorfaces');
            NegativeElectrodefaces = cell12facetbl.get('NegativeElectrodefaces');
            NegativeCurrentCollectorcells = cell12facetbl.get('NegativeCurrentCollectorcells');
            NegativeElectrodecells = cell12facetbl.get('NegativeElectrodecells');            
            
            compnames = {'NegativeCurrentCollector', 'NegativeElectrode'};
            coupTerm = couplingTerm('NegativeCurrentCollector-NegativeElectrode', compnames);
            coupTerm.couplingfaces =  [NegativeCurrentCollectorfaces, NegativeElectrodefaces];
            coupTerm.couplingcells = [NegativeCurrentCollectorcells, NegativeElectrodecells];

        end
        
        function coupTerm = setupPositiveCurrentCollectorBcCoupTerm(params)

            PositiveCurrentCollector = params.PositiveCurrentCollector;
            G = PositiveCurrentCollector.G;

            % We pick up the faces at the top of CPositiveCurrentCollector
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf - myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'PositiveCurrentCollector'};
            coupTerm = couplingTerm('bc-PositiveCurrentCollector', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupPositiveCurrentCollectorPositiveElectrodeCoupTerm(params)

            PositiveElectrode = params.PositiveElectrode;
            PositiveCurrentCollector = params.PositiveCurrentCollector;

            GPositiveElectrode = PositiveElectrode.G;
            GPositiveCurrentCollector = PositiveCurrentCollector.G;

            G = GPositiveElectrode.mappings.parentGrid;

            PositiveElectrodetbls = setupSimpleTables(GPositiveElectrode);
            PositiveCurrentCollectortbls = setupSimpleTables(GPositiveCurrentCollector);
            tbls = setupSimpleTables(G);            
            
            PositiveElectrodecelltbl = PositiveElectrodetbls.celltbl;
            PositiveElectrodecelltbl = PositiveElectrodecelltbl.addInd('globcells', GPositiveElectrode.mappings.cellmap);
            PositiveElectrodefacetbl = PositiveElectrodetbls.facetbl;
            PositiveElectrodefacetbl = PositiveElectrodefacetbl.addInd('globfaces', GPositiveElectrode.mappings.facemap);

            PositiveElectrodecellfacetbl = PositiveElectrodetbls.cellfacetbl;
            PositiveElectrodecellfacetbl = crossIndexArray(PositiveElectrodecellfacetbl, PositiveElectrodecelltbl, {'cells'});
            PositiveElectrodecellfacetbl = crossIndexArray(PositiveElectrodecellfacetbl, PositiveElectrodefacetbl, {'faces'});
            
            PositiveCurrentCollectorcelltbl = PositiveCurrentCollectortbls.celltbl;
            PositiveCurrentCollectorcelltbl = PositiveCurrentCollectorcelltbl.addInd('globcells', GPositiveCurrentCollector.mappings.cellmap);
            PositiveCurrentCollectorfacetbl = PositiveCurrentCollectortbls.facetbl;
            PositiveCurrentCollectorfacetbl = PositiveCurrentCollectorfacetbl.addInd('globfaces', GPositiveCurrentCollector.mappings.facemap);
            
            PositiveCurrentCollectorcellfacetbl = PositiveCurrentCollectortbls.cellfacetbl;
            PositiveCurrentCollectorcellfacetbl = crossIndexArray(PositiveCurrentCollectorcellfacetbl, PositiveCurrentCollectorcelltbl, {'cells'});
            PositiveCurrentCollectorcellfacetbl = crossIndexArray(PositiveCurrentCollectorcellfacetbl, PositiveCurrentCollectorfacetbl, {'faces'});
            
            gen = CrossIndexArrayGenerator();
            gen.tbl1 = PositiveElectrodecellfacetbl;
            gen.tbl2 = PositiveCurrentCollectorcellfacetbl;
            gen.replacefds1 = {{'cells', 'PositiveElectrodecells'}, {'faces', 'PositiveElectrodefaces'}, {'globcells', 'PositiveElectrodeglobcells'}};
            gen.replacefds2 = {{'cells', 'PositiveCurrentCollectorcells'}, {'faces', 'PositiveCurrentCollectorfaces'}, {'globcells', 'PositiveCurrentCollectorglobcells'}};
            gen.mergefds = {'globfaces'};
            
            cell12facetbl = gen.eval();

            PositiveCurrentCollectorfaces = cell12facetbl.get('PositiveCurrentCollectorfaces');
            PositiveElectrodefaces = cell12facetbl.get('PositiveElectrodefaces');
            PositiveCurrentCollectorcells = cell12facetbl.get('PositiveCurrentCollectorcells');
            PositiveElectrodecells = cell12facetbl.get('PositiveElectrodecells');            
            
            compnames = {'PositiveCurrentCollector', 'PositiveElectrode'};
            coupTerm = couplingTerm('PositiveCurrentCollector-PositiveElectrode', compnames);
            coupTerm.couplingfaces =  [PositiveCurrentCollectorfaces, PositiveElectrodefaces];
            coupTerm.couplingcells = [PositiveCurrentCollectorcells, PositiveElectrodecells];
            
        end

        function coupTerm = setupNegativeElectrodeElectrolyteCoupTerm(params)
            
            NegativeElectrode = params.NegativeElectrode;
            Electrolyte = params.Electrolyte;
            
            GNegativeElectrode = NegativeElectrode.G;
            GElectrolyte = Electrolyte.G;
            
            % parent Grid
            G = GNegativeElectrode.mappings.parentGrid;
            
            % All the cells from NegativeElectrode are coupled with Electrolyte
            cells1 = (1 : GNegativeElectrode.cells.num)';
            pcells = GNegativeElectrode.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(GElectrolyte.mappings.cellmap) = (1 : GElectrolyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'NegativeElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('NegativeElectrode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
        end
        
        function coupTerm = setupPositiveElectrodeElectrolyteCoupTerm(params)
            
            PositiveElectrode = params.PositiveElectrode;
            Electrolyte = params.Electrolyte;
            
            GPositiveElectrode = PositiveElectrode.G;
            GElectrolyte = Electrolyte.G;
            
            % parent Grid
            G = GPositiveElectrode.mappings.parentGrid;
            
            % All the cells from PositiveElectrode are coupled with Electrolyte
            cells1 = (1 : GPositiveElectrode.cells.num)';
            pcells = GPositiveElectrode.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(GElectrolyte.mappings.cellmap) = (1 : GElectrolyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'PositiveElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('PositiveElectrode-Electrolyte', compnames);
            coupTerm.couplingcells = [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling between faces
            
        end

    end
    
end

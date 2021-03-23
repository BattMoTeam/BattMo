classdef BatteryInputParams2D < BatteryInputParams
    
    methods
        
        function params = setupVariousParams(params)
            
            params.SOC  = 0.5;
            params.T    = 298.15;
            params.J    = 0.1;
            params.Ucut = 2;
        
        end

        function params = setupSubModels(params)
        % Abbreviations used in this function:
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
        % ncc   : NegaticeCurrentCollector
        % pcc   : PositiveCurrentCollector
            
            fac = 1;
            
            sepnx  = 30*fac;
            ne_nx  = 30*fac;
            pe_nx  = 30*fac;
            ncc_nx = 20*fac;
            pcc_nx = 20*fac;

            nxs = [ncc_nx; ne_nx; sepnx; pe_nx; pcc_nx];
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

            istart = ncc_nx + 1;
            ni = ne_nx + sepnx + pe_nx;
            cells = pickTensorCells(istart, ni, nx, ny);
            params.Electrolyte = orgLiPF6(G, cells);

            %% setup NegativeElectrode
            istart = ncc_nx + 1;
            cells = pickTensorCells(istart, ne_nx, nx, ny);
            params.NegativeElectrode = GraphiteElectrode(G, cells);

            %% setup PositiveElectrode
            istart = ncc_nx + ne_nx + sepnx + 1;
            cells = pickTensorCells(istart, pe_nx, nx, ny);
            params.PositiveElectrode = NMC111Electrode(G, cells);

            %% setup NegativeCurrentCollector
            istart = 1;
            cells = pickTensorCells(istart, ncc_nx, nx, ny);
            params.NegativeCurrentCollector = CurrentCollector(G, cells);

            %% setup PositiveCurrentCollector
            istart = ncc_nx + ne_nx + sepnx + pe_nx + 1;
            cells = pickTensorCells(istart, pcc_nx, nx, ny);
            params.PositiveCurrentCollector = CurrentCollector(G, cells);

            %% setup sep
            istart = ncc_nx + ne_nx + 1;
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
        % Abbreviations used in this function:
        % ne    : NegativeElectrode
        % ncc   : NegativeCurrentCollector
            
            ne = params.NegativeElectrode;
            ncc = params.NegativeCurrentCollector;

            G_ne = ne.G;
            G_ncc = ncc.G;
            
            G = G_ne.mappings.parentGrid;

            tbls_ne  = setupSimpleTables(G_ne);
            tbls_ncc = setupSimpleTables(G_ncc);
            tbls     = setupSimpleTables(G);
            
            celltbl_ne = tbls_ne.celltbl;
            celltbl_ne = celltbl_ne.addInd('globcells', G_ne.mappings.cellmap);
            facetbl_ne = tbls_ne.facetbl;
            facetbl_ne = facetbl_ne.addInd('globfaces', G_ne.mappings.facemap);

            cellfacetbl_ne = tbls_ne.cellfacetbl;
            cellfacetbl_ne = crossIndexArray(cellfacetbl_ne, celltbl_ne, {'cells'});
            cellfacetbl_ne = crossIndexArray(cellfacetbl_ne, facetbl_ne, {'faces'});
            
            celltbl_ncc = tbls_ncc.celltbl;
            celltbl_ncc = celltbl_ncc.addInd('globcells', G_ncc.mappings.cellmap);
            facetbl_ncc = tbls_ncc.facetbl;
            facetbl_ncc = facetbl_ncc.addInd('globfaces', G_ncc.mappings.facemap);
            
            cellfacetbl_ncc = tbls_ncc.cellfacetbl;
            cellfacetbl_ncc = crossIndexArray(cellfacetbl_ncc, celltbl_ncc, {'cells'});
            cellfacetbl_ncc = crossIndexArray(cellfacetbl_ncc, facetbl_ncc, {'faces'});
            
            gen = CrossIndexArrayGenerator();
            gen.tbl1 = cellfacetbl_ne;
            gen.tbl2 = cellfacetbl_ncc;
            gen.replacefds1 = {{'cells', 'cells_ne'}, {'faces', 'faces_ne'}};
            gen.replacefds2 = {{'cells', 'cells_ncc'}, {'faces', 'faces_ncc'}};
            gen.mergefds = {'globfaces'};
            
            cell12facetbl = gen.eval();

            faces_ncc = cell12facetbl.get('faces_ncc');
            faces_ne  = cell12facetbl.get('faces_ne');
            cells_ncc = cell12facetbl.get('cells_ncc');
            cells_ne  = cell12facetbl.get('cells_ne');
            
            compnames = {'NegativeCurrentCollector', 'NegativeElectrode'};
            coupTerm = couplingTerm('NegativeCurrentCollector-NegativeElectrode', compnames);
            coupTerm.couplingfaces =  [faces_ncc, faces_ne];
            coupTerm.couplingcells = [cells_ncc, cells_ne];

        end
        
        function coupTerm = setupPositiveCurrentCollectorBcCoupTerm(params)
        % Abbreviations used in this function:
        % pcc   : PositiveCurrentCollector
            
            pcc = params.PositiveCurrentCollector;
            G = pcc.G;

            % We pick up the faces at the top of PositiveCurrentCollector
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
        % Abbreviations used in this function:
        % pe    : NegativeElectrode
        % pcc   : NegativeCurrentCollector
            
            pe = params.PositiveElectrode;
            pcc = params.PositiveCurrentCollector;

            G_pe = pe.G;
            G_pcc = pcc.G;

            G = G_pe.mappings.parentGrid;

            tbls_pe  = setupSimpleTables(G_pe);
            tbls_pcc = setupSimpleTables(G_pcc);
            tbls     = setupSimpleTables(G);            
            
            celltbl_pe = tbls_pe.celltbl;
            celltbl_pe = celltbl_pe.addInd('globcells', G_pe.mappings.cellmap);
            facetbl_pe = tbls_pe.facetbl;
            facetbl_pe = facetbl_pe.addInd('globfaces', G_pe.mappings.facemap);

            cellfacetbl_pe = tbls_pe.cellfacetbl;
            cellfacetbl_pe = crossIndexArray(cellfacetbl_pe, celltbl_pe, {'cells'});
            cellfacetbl_pe = crossIndexArray(cellfacetbl_pe, facetbl_pe, {'faces'});
            
            celltbl_pcc = tbls_pcc.celltbl;
            celltbl_pcc = celltbl_pcc.addInd('globcells', G_pcc.mappings.cellmap);
            facetbl_pcc = tbls_pcc.facetbl;
            facetbl_pcc = facetbl_pcc.addInd('globfaces', G_pcc.mappings.facemap);
            
            cellfacetbl_pcc = tbls_pcc.cellfacetbl;
            cellfacetbl_pcc = crossIndexArray(cellfacetbl_pcc, celltbl_pcc, {'cells'});
            cellfacetbl_pcc = crossIndexArray(cellfacetbl_pcc, facetbl_pcc, {'faces'});
            
            gen = CrossIndexArrayGenerator();
            gen.tbl1 = cellfacetbl_pe;
            gen.tbl2 = cellfacetbl_pcc;
            gen.replacefds1 = {{'cells', 'cells_pe'}, {'faces', 'faces_pe'}};
            gen.replacefds2 = {{'cells', 'cells_pcc'}, {'faces', 'faces_pcc'}};
            gen.mergefds = {'globfaces'};
            
            cell12facetbl = gen.eval();

            faces_pcc = cell12facetbl.get('faces_pcc');
            faces_pe  = cell12facetbl.get('faces_pe');
            cells_pcc = cell12facetbl.get('cells_pcc');
            cells_pe  = cell12facetbl.get('cells_pe');
            
            compnames = {'PositiveCurrentCollector', 'PositiveElectrode'};
            coupTerm = couplingTerm('PositiveCurrentCollector-PositiveElectrode', compnames);
            coupTerm.couplingfaces =  [faces_pcc, faces_pe];
            coupTerm.couplingcells = [cells_pcc, cells_pe];
            
        end

        function coupTerm = setupNegativeElectrodeElectrolyteCoupTerm(params)
        % Abbreviations used in this function:
        % elyte : Electrolyte
        % ne    : NegativeElectrode

            ne = params.NegativeElectrode;
            elyte = params.Electrolyte;
            
            G_ne = ne.G;
            G_elyte = elyte.G;
            
            % parent Grid
            G = G_ne.mappings.parentGrid;
            
            % All the cells from NegativeElectrode are coupled with Electrolyte
            cells1 = (1 : G_ne.cells.num)';
            pcells = G_ne.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'NegativeElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('NegativeElectrode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
        end
        
        function coupTerm = setupPositiveElectrodeElectrolyteCoupTerm(params)
        % Abbreviations used in this function:
        % elyte : Electrolyte
        % pe    : NegativeElectrode
            
            pe = params.PositiveElectrode;
            elyte = params.Electrolyte;
            
            G_pe = pe.G;
            G_elyte = elyte.G;
            
            % parent Grid
            G = G_pe.mappings.parentGrid;
            
            % All the cells from PositiveElectrode are coupled with Electrolyte
            cells1 = (1 : G_pe.cells.num)';
            pcells = G_pe.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'PositiveElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('PositiveElectrode-Electrolyte', compnames);
            coupTerm.couplingcells = [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling between faces
            
        end

    end
    
end

classdef ElectrodeInputParams2D < ElectrodeInputParams
    
    
    methods
        function paramobj = ElectrodeInputParams(params)
        % params struct should contain valid fields for ComponentInputParams
        %
        % and valid fields for the methods (see implementation of those methods)
        %
        % - setupElectrodeActiveComponent
        % - setupCurrentCollector
        % - setupCurrentCollectorBcCoupTerm
        % - setupCurrentCollectorElectrodeActiveComponentCoupTerm
            
            paramobj = paramobj@ComponentInputParams(params);
            paramobj.eac = paramobj.setupElectrodeActiveComponent(params);
            paramobj.cc  = paramobj.setupCurrentCollector(params);
            paramobj.couplingTerms = paramobj.setupCouplingTerms(params);
            
        end

        function coupTerm = setupCurrentCollectorBcCoupTerm(paramobj, params)
        % Abbreviations used in this function:
        % cc   : CurrentCollector
            
            cc_paramobj = paramobj.cc;
            G = cc_paramobj.G;

            % We pick up the faces at the top of PositiveCurrentCollector
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf - myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'CurrentCollector'};
            coupTerm = couplingTerm('bc-CurrentCollector', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;
            
        end
        
        function coupTerm = setupCurrentCollectorElectrodeActiveComponentCoupTerm(paramobj, params)
        % Abbreviations used in this function:
        % eac : ElectrodeActiveComponent
        % cc  : CurrentCollector
            
            eac_paramobj = paramobj.eac;
            cc_paramobj = paramobj.cc;

            G_eac = eac_paramobj.G;
            G_cc = cc_paramobj.G;
            
            G = G_eac.mappings.parentGrid;

            tbls_eac = setupSimpleTables(G_eac);
            tbls_cc  = setupSimpleTables(G_cc);
            tbls     = setupSimpleTables(G);
            
            celltbl_eac = tbls_eac.celltbl;
            celltbl_eac = celltbl_eac.addInd('globcells', G_eac.mappings.cellmap);
            facetbl_eac = tbls_eac.facetbl;
            facetbl_eac = facetbl_eac.addInd('globfaces', G_eac.mappings.facemap);

            cellfacetbl_eac = tbls_eac.cellfacetbl;
            cellfacetbl_eac = crossIndexArray(cellfacetbl_eac, celltbl_eac, {'cells'});
            cellfacetbl_eac = crossIndexArray(cellfacetbl_eac, facetbl_eac, {'faces'});
            
            celltbl_cc = tbls_cc.celltbl;
            celltbl_cc = celltbl_cc.addInd('globcells', G_cc.mappings.cellmap);
            facetbl_cc = tbls_cc.facetbl;
            facetbl_cc = facetbl_cc.addInd('globfaces', G_cc.mappings.facemap);
            
            cellfacetbl_cc = tbls_cc.cellfacetbl;
            cellfacetbl_cc = crossIndexArray(cellfacetbl_cc, celltbl_cc, {'cells'});
            cellfacetbl_cc = crossIndexArray(cellfacetbl_cc, facetbl_cc, {'faces'});
            
            gen = CrossIndexArrayGenerator();
            gen.tbl1 = cellfacetbl_eac;
            gen.tbl2 = cellfacetbl_cc;
            gen.replacefds1 = {{'cells', 'cells_eac'}, {'faces', 'faces_eac'}};
            gen.replacefds2 = {{'cells', 'cells_cc'}, {'faces', 'faces_cc'}};
            gen.mergefds = {'globfaces'};
            
            cell12facetbl = gen.eval();

            faces_cc  = cell12facetbl.get('faces_cc');
            faces_eac = cell12facetbl.get('faces_eac');
            cells_cc  = cell12facetbl.get('cells_cc');
            cells_eac = cell12facetbl.get('cells_eac');
            
            compnames = {'CurrentCollector', 'ElectrodeActiveComponent'}
            coupTerm = couplingTerm('CurrentCollector-ElectrodeActiveComponent', compnames);
            coupTerm.couplingfaces =  [faces_cc, faces_eac];
            coupTerm.couplingcells = [cells_cc, cells_eac];

        end

    end
    
end

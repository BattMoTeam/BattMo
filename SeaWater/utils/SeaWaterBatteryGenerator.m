classdef SeaWaterBatteryGenerator < BatteryGenerator

    methods

        function gen = SeaWaterBatteryGenerator()
            gen = gen@BatteryGenerator();
        end
        
        function paramobj = setupBatteryInputParams(gen, paramobj, params)
            [paramobj, gen] = gen.setupGrid(paramobj, params);
            paramobj.Electrolyte = gen.setupElectrolyte(paramobj.Electrolyte, params);
            paramobj = gen.setupElectrodes(paramobj, params);
            paramobj = gen.setupElectrodeElectrolyteCoupTerm(paramobj);
        end
        
        function paramobj = setupElectrolyte(gen, paramobj, params)
            paramobj = gen.setupElectrolyteGrid(paramobj, params);
        end
        

        function paramobj = setupElectrodes(gen, paramobj, params)
           
            ctde = 'Cathode';
            ande = 'Anode';
            % setup Cathode
            paramobj.(ctde) = gen.setupElectrodeGrid(paramobj.(ctde), params.(ctde));
            paramobj.(ctde) = gen.setupElectrodeBcCoupTerm(paramobj.(ctde), params.(ctde));
            % setup Anode 
            paramobj.(ande) = gen.setupElectrodeGrid(paramobj.(ande), params.(ande));
            paramobj.(ande) = gen.setupElectrodeBcCoupTerm(paramobj.(ande), params.(ande));
            
        end
        
        function paramobj = setupElectrodeBcCoupTerm(gen, paramobj, params)
            
            % default setup
            compname = 'Electrode';
            compnames = {compname};
            coupname = sprintf('Exterior-%s', compname);
            coupTerm = couplingTerm(coupname, compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;
            
            paramobj.externalCouplingTerm = coupTerm;
            
        end        

                

        function paramobj = setupElectrodeElectrolyteCoupTerm(gen, paramobj, params)
        % paramobj is instance of SeaWaterBatteryInputParams
        % setup paramobj.couplingTerms
            ctde  = 'Cathode';
            ande  = 'Anode';
            elyte = 'Electrolyte';
            
            couplingTerms = {};
            
            G_ctde = paramobj.(ctde).G;
            G_elyte = paramobj.(elyte).G;
            
            % parent Grid
            G = G_ctde.mappings.parentGrid;
            
            % All the cells from Cathode are coupled with Electrolyte
            cells1 = (1 : G_ctde.cells.num)';
            pcells = G_ctde.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'Cathode', 'Electrolyte'};
            coupTerm = couplingTerm('Cathode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
            couplingTerms{end + 1} = coupTerm;
            
            G_ande = paramobj.(ande).G;
            G_elyte = paramobj.(elyte).G;
            
            % parent Grid
            G = G_ande.mappings.parentGrid;
            
            % All the cells from Anode are coupled with Electrolyte
            cells1 = (1 : G_ande.cells.num)';
            pcells = G_ande.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'Anode', 'Electrolyte'};
            coupTerm = couplingTerm('Anode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
            couplingTerms{end + 1} = coupTerm;
            
            paramobj.couplingTerms = couplingTerms;

        end

        function paramobj = setupCurrentCollectorElectrodeActiveComponentCoupTerm(gen, paramobj, params)
        % paramobj is instance of ElectrodeInputParams
        % setup paramobj.couplingTerm
            eac = 'ElectrodeActiveComponent';
            cc  = 'CurrentCollector';
            
            compnames = {'CurrentCollector', 'ElectrodeActiveComponent'};
            coupTerm = couplingTerm('CurrentCollector-ElectrodeActiveComponent', compnames);
            
            G_eac = paramobj.(eac).G;
            G_cc = paramobj.(cc).G;
            
            cctbl.faces = (1 : G_cc.faces.num)';
            cctbl.globfaces = G_cc.mappings.facemap;
            cctbl = IndexArray(cctbl);
            
            
            eactbl.faces = (1 : G_eac.faces.num)';
            eactbl.globfaces = G_eac.mappings.facemap;
            eactbl = IndexArray(eactbl);

            gen = CrossIndexArrayGenerator();
            gen.tbl1 = cctbl;
            gen.tbl2 = eactbl;
            gen.replacefds1 = {{'faces', 'faces1'}};
            gen.replacefds2 = {{'faces', 'faces2'}};
            gen.mergefds = {'globfaces'};
            tbl = gen.eval();
            
            cc_coupfaces = tbl.get('faces1');
            eac_coupfaces = tbl.get('faces2');
            
            cc_coupcells = sum(G_cc.faces.neighbors(cc_coupfaces, :), 2);
            eac_coupcells = sum(G_eac.faces.neighbors(eac_coupfaces, :), 2);
            
            coupTerm.couplingfaces = [cc_coupfaces, eac_coupfaces];
            coupTerm.couplingcells = [cc_coupcells, eac_coupcells];
            
            paramobj.couplingTerm = coupTerm;
            
        end

        function paramobj = setupCurrentCollectorBcCoupTerm(gen, paramobj, params)
        % paramobj is instance of CurrentCollectorInputParams
        % setup paramobj.couplingTerm
            
            % default setup
            compnames = {'CurrentCollector'};
            coupTerm = couplingTerm('bc-CurrentCollector', compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;
            
            paramobj.couplingTerm = coupTerm;
            
        end

        
    end
    
    
end



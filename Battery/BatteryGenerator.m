classdef BatteryGenerator
% Object that add grids and coupling terms to a paramobj instance of BatteryInputParams (through method 
% updateBatteryInputParams)
%
% This class goes through the whole setup and is meant to be used as a base class.
%
% This class (or rather an subclass of it) can also be used to setup values of the paramobj instance that is sent to it
% (just overload the updateBatteryInputParams method by adding the desired setup)
%
% Example : BatteryGenerator1D.m, BatteryGenerator2D.m

    properties
        % Global grid 
        % It is stored here because is shared by many setup functions
        % It is setup by function setupBatteryGrid
        G
    end

    methods
        
        function paramobj = updateBatteryInputParams(gen, paramobj)
            error('virtual function - should call setupBatteryInputParams with some argument for params');
        end
        
        function paramobj = setupBatteryInputParams(gen, paramobj, params)
        % main function : add grid and coupling to paramobj structure
            [paramobj, gen] = gen.setupGrid(paramobj, params);
            paramobj.elyte = gen.setupElectrolyte(paramobj.elyte, params);
            paramobj = gen.setupElectrodes(paramobj, params);
            paramobj = gen.setupThermalModel(paramobj, params);
            paramobj = gen.setupElectrodeElectrolyteCoupTerm(paramobj);
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
        % setup paramobj.G and update gen.G 
            error('virtual function');
        end
        
        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
            paramobj.thermal.G = gen.G;
            coupTerm = couplingTerm('ThermalConvectiveCooling', {'ThermalModel'});
            coupTerm.couplingcells = params.couplingcells;
            coupTerm.couplingfaces = params.couplingfaces;
            paramobj.thermal.couplingTerm = coupTerm;
        end
        
        function paramobj = setupElectrolyte(gen, paramobj, params)
        % paramobj is instance of ElectrolyteInputParams
            paramobj = gen.setupElectrolyteGrid(paramobj, params);
            paramobj.sep = gen.setupSeparatorGrid(paramobj.sep, params.sep);
        end

        
        function paramobj = setupElectrolyteGrid(gen, paramobj, params)
        % paramobj is instance of ElectrolyteInputParams
        % setup paramobj.G

            % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
            
        end
        
        function paramobj = setupSeparatorGrid(gen, paramobj, params)
        % paramobj is instance of SeparatorInputParams
        % setup paramobj.G
            
           % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
            
        end

        function paramobj = setupElectrodes(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
            
            % setup Negative Electrode 
            paramobj.ne = gen.setupElectrode(paramobj.ne, params.ne);
            % setup Positive Electrode
            paramobj.pe = gen.setupElectrode(paramobj.pe, params.pe);
            
        end
                
        function paramobj = setupElectrode(gen, paramobj, params)
        % paramobj is instance of ElectrodeInputParams
            
            % setup Electrode grid
            paramobj = gen.setupElectrodeGrid(paramobj, params);
            % setup Electrode active component (eac)
            paramobj.eac = gen.setupElectrodeActiveComponentGrid(paramobj.eac, params.eac);
            % setup current collector (cc)
            paramobj.cc = gen.setupCurrentCollector(paramobj.cc, params.cc);
            % setup coupling term between eac and cc
            paramobj = gen.setupCurrentCollectorElectrodeActiveComponentCoupTerm(paramobj, params);
            
        end
        
        function paramobj = setupElectrodeGrid(gen, paramobj, params)
        % paramobj is instance of ElectrodeInputParams
        % setup paramobj.G
            
            % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
        end

        function paramobj = setupElectrodeActiveComponentGrid(gen, paramobj, params)
        % paramobj is instance of ElectrodeActiveComponentInputParams
        % setup paramobj.G
            
            % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
        end
        
        function paramobj = setupCurrentCollector(gen, paramobj, params)
        % paramobj is instance of CurrentCollectorInputParams
            paramobj = gen.setupCurrentCollectorGrid(paramobj, params);
            paramobj = gen.setupCurrentCollectorBcCoupTerm(paramobj, params);
        end

        function paramobj = setupCurrentCollectorGrid(gen, paramobj, params)
        % paramobj is instance of CurrentCollectorInputParams
        % setup paramobj.G

           % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
        end       

        function paramobj = setupElectrodeElectrolyteCoupTerm(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
        % setup paramobj.couplingTerms
            
            couplingTerms = {};
            
            G_ne = paramobj.ne.eac.G;
            G_elyte = paramobj.elyte.G;
            
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
            
            couplingTerms{end + 1} = coupTerm;
            
            G_pe = paramobj.pe.eac.G;
            G_elyte = paramobj.elyte.G;
            
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
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
            couplingTerms{end + 1} = coupTerm;
            
            paramobj.couplingTerms = couplingTerms;

        end

        function paramobj = setupCurrentCollectorElectrodeActiveComponentCoupTerm(gen, paramobj, params)
        % paramobj is instance of ElectrodeInputParams
        % setup paramobj.couplingTerm
            
            compnames = {'CurrentCollector', 'ElectrodeActiveComponent'};
            coupTerm = couplingTerm('CurrentCollector-ElectrodeActiveComponent', compnames);
            
            G_eac = paramobj.eac.G;
            G_cc = paramobj.cc.G;
            
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



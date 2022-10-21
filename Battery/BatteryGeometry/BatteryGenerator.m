classdef BatteryGenerator
% Base class that add grids and coupling terms to a paramobj instance of BatteryInputParams (through method 
% updateBatteryInputParams)
%
% This class goes through the whole grid setup and is meant to be used as a base class.
%
% Example : BatteryGenerator1D.m, BatteryGenerator2D.m, BatteryGenerator3D

    properties
        % Global grid 
        % It is stored here because is shared by many setup functions
        % It is setup by function setupBatteryGrid
        G
        
    end

    methods
        
        function [paramobj, gen] = updateBatteryInputParams(gen, paramobj, params)
        % this function is the main class function as it returns an updated paramobj object with grid structure
            error('virtual function');
        end
        
        function [paramobj, gen] = setupBatteryInputParams(gen, paramobj, params)
        % main function : add grid and coupling to paramobj structure
            [paramobj, gen] = gen.setupGrid(paramobj, params);
            % We check if params contains some Electrolyte field 
            if isfield(params, 'Electrolyte')
                params_electrolyte = params.Electrolyte;
            else
                params_electrolyte = [];
            end
            paramobj.Electrolyte = gen.setupElectrolyte(paramobj.Electrolyte, params_electrolyte);
            paramobj = gen.setupElectrodes(paramobj, params);
            if gen.use_thermal
                params_thermal = [];
                if isfield(params, 'ThermalModel')
                    params_thermal = params.ThermalModel;
                end
                paramobj = gen.setupThermalModel(paramobj, params_thermal);
            end
            paramobj = gen.setupElectrodeElectrolyteCoupTerm(paramobj);
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
        % setup paramobj.G and update gen.G 
            error('virtual function');
        end
        
        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
            paramobj.ThermalModel.G = gen.G;
            coupTerm = couplingTerm('ThermalConvectiveCooling', {'ThermalModel'});
            coupTerm.couplingcells = params.couplingcells;
            coupTerm.couplingfaces = params.couplingfaces;
            paramobj.ThermalModel.couplingTerm = coupTerm;
        end
        
        function paramobj = setupElectrolyte(gen, paramobj, params)
        % paramobj is instance of ElectrolyteInputParams
            paramobj = gen.setupElectrolyteGrid(paramobj, params);
            paramobj.Separator = gen.setupSeparatorGrid(paramobj.Separator, params.Separator);
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
            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            eldes = {ne, pe};

            % We add the electrode type to the params structure (This information can be used by derived classes and may
            % then simplify setup)
            params.(ne).electrodeType = ne;
            params.(pe).electrodeType = pe;
            
            % setup Negative Electrode
            paramobj.(ne) = gen.setupElectrode(paramobj.(ne), params.(ne));
            % setup Positive Electrode
            paramobj.(pe) = gen.setupElectrode(paramobj.(pe), params.(pe));

            
        end
                
        function paramobj = setupElectrode(gen, paramobj, params)
        % paramobj is instance of ElectrodeInputParams
        
            % shorthands 
            am = 'ActiveMaterial';
            cc = 'CurrentCollector';
            
            % setup Electrode grid
            paramobj = gen.setupElectrodeGrid(paramobj, params);
            % setup Electrode active component (am)
            paramobj.(am) = gen.setupActiveMaterialGrid(paramobj.(am), params.(am));
            % setup current collector (cc)
            if paramobj.include_current_collector
                % We add the electrode type to the params structure. (This information can be used by derived classes
                % and may then simplify setup)
                params.(cc).electrodeType = params.electrodeType;
                paramobj.(cc) = gen.setupCurrentCollector(paramobj.(cc), params.(cc));
                % setup coupling term between am and cc
                paramobj = gen.setupCurrentCollectorActiveMaterialCoupTerm(paramobj, params);
            else
                paramobj.(am) = gen.setupActiveMaterialBcCoupTerm(paramobj.(am), params.(am));
            end
            
        end
        
        function paramobj = setupElectrodeGrid(gen, paramobj, params)
        % paramobj is instance of ElectrodeInputParams
        % setup paramobj.G
            
            % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
        end

        function paramobj = setupActiveMaterialGrid(gen, paramobj, params)
        % paramobj is instance of ActiveMaterialInputParams
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
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            am   = 'ActiveMaterial';
            
            couplingTerms = {};
            
            G_ne = paramobj.(ne).(am).G;
            G_elyte = paramobj.(elyte).G;
            
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
            
            G_pe = paramobj.(pe).(am).G;
            G_elyte = paramobj.(elyte).G;
            
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

        function paramobj = setupCurrentCollectorActiveMaterialCoupTerm(gen, paramobj, params)
        % paramobj is instance of ElectrodeInputParams
        % setup paramobj.couplingTerm
            am = 'ActiveMaterial';
            cc  = 'CurrentCollector';
            
            compnames = {'CurrentCollector', 'ActiveMaterial'};
            coupTerm = couplingTerm('CurrentCollector-ActiveMaterial', compnames);
            
            G_am = paramobj.(am).G;
            G_cc = paramobj.(cc).G;
            
            cctbl.faces = (1 : G_cc.faces.num)';
            cctbl.globfaces = G_cc.mappings.facemap;
            cctbl = IndexArray(cctbl);
            
            
            eactbl.faces = (1 : G_am.faces.num)';
            eactbl.globfaces = G_am.mappings.facemap;
            eactbl = IndexArray(eactbl);

            gen = CrossIndexArrayGenerator();
            gen.tbl1 = cctbl;
            gen.tbl2 = eactbl;
            gen.replacefds1 = {{'faces', 'faces1'}};
            gen.replacefds2 = {{'faces', 'faces2'}};
            gen.mergefds = {'globfaces'};
            tbl = gen.eval();
            
            cc_coupfaces = tbl.get('faces1');
            am_coupfaces = tbl.get('faces2');
            
            cc_coupcells = sum(G_cc.faces.neighbors(cc_coupfaces, :), 2);
            am_coupcells = sum(G_am.faces.neighbors(am_coupfaces, :), 2);
            
            coupTerm.couplingfaces = [cc_coupfaces, am_coupfaces];
            coupTerm.couplingcells = [cc_coupcells, am_coupcells];
            
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
            
            paramobj.externalCouplingTerm = coupTerm;
            
        end

        function paramobj = setupActiveMaterialBcCoupTerm(gen, paramobj, params)
            
            % default setup
            compname = 'ActiveMaterial';
            compnames = {compname};
            coupname = sprintf('Exterior-%s', compname);
            coupTerm = couplingTerm(coupname, compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;
            
            paramobj.externalCouplingTerm = coupTerm;
            
        end        

    end
    
    
end





%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}

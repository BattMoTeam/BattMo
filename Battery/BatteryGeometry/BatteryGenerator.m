classdef BatteryGenerator
% Base class that add grids and coupling terms to a paramobj instance of BatteryInputParams (through method 
% updateBatteryInputParams)
%
% This class goes through the whole grid setup and is meant to be used as a base class.
%
% Example : :class:`BatteryGenerator1D <BatteryGeometry.BatteryGenerator1D>`, :class:`BatteryGenerator2D <BatteryGeometry.BatteryGenerator2D>`, :class:`BatteryGenerator3D <BatteryGeometry.BatteryGenerator3D>`

    properties
        
        % Parent grid, instance of Grid
        parentGrid
        
    end

    methods
        
        function [paramobj, gen] = updateBatteryInputParams(gen, paramobj, params)
        %
        % This function is the main class function as it returns an updated :code:`paramobj` object with grid structure
        %
            
            error('virtual function');
            
        end
        
        function [paramobj, gen] = setupBatteryInputParams(gen, paramobj, params)
        % This is the main function which should be called by the derived class at the end of the :meth:`updateBatteryInputParams` method
        % The function set up the grid and the coupling terms and add those in the :code:`paramobj` structure which is
        % an instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`
        % 
            
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
        % setup :code:`paramobj.G` and update  :attr:`G`
        % Here, :code:`paramobj` is an instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`
            
            error('virtual function');
            
        end
        
        function paramobj = setupThermalModel(gen, paramobj, params)
        % Method that setups the grid and the coupling for the thermal model
            
            G = genSubGrid(gen.parentGrid, (1 : gen.parentGrid.topology.cells.num)');
            G = setupCellFluxOperators(G);
            coupTerm = couplingTerm('ThermalConvectiveCooling', {'ThermalModel'});
            coupTerm.couplingcells = params.couplingcells;
            coupTerm.couplingfaces = params.couplingfaces;
            paramobj.ThermalModel.couplingTerm = coupTerm;
            
        end
        
        function paramobj = setupElectrolyte(gen, paramobj, params)
        % Method that setups the grid and the coupling for the electrolyte model
        % Here, :code:`paramobj` is instance of :class:`ElectrolyteInputParams <Electrochemistry.ElectrolyteInputParams>`
            
            paramobj = gen.setupElectrolyteGrid(paramobj, params);
            paramobj.Separator = gen.setupSeparatorGrid(paramobj.Separator, params.Separator);
        end

        
        function paramobj = setupElectrolyteGrid(gen, paramobj, params)
        % Setup the grid for the electrolyte
        % Here, :code:`paramobj` is instance of :class:`ElectrolyteInputParams <Electrochemistry.ElectrolyteInputParams>`
            
            % Default setup
            paramobj.G = genSubGrid(gen.parentGrid, params.cellind);
            
        end
        
        function paramobj = setupSeparatorGrid(gen, paramobj, params)
        % Setup the grid for the separator
        % Here, :code:`paramobj` is instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`
            
           % Default setup
            paramobj.G = genSubGrid(gen.parentGrid, params.cellind);
            
        end

        function paramobj = setupElectrodes(gen, paramobj, params)
        % Method that setups the grid and the coupling terms for both electrodes

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
        % Method that setups the grid and the coupling terms for an electrode
        % Here, :code:`paramobj` is instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`
            
            % shorthands 
            am = 'ActiveMaterial';
            cc = 'CurrentCollector';
            
            % setup Electrode grid
            paramobj = gen.setupElectrodeGrid(paramobj, params);
            % setup Electrode active component (am)
            params.(am).electrode_case = paramobj.electrode_case;
            paramobj.(am) = gen.setupActiveMaterialGrid(paramobj.(am), params.(am));
            % setup current collector (cc)
            if paramobj.include_current_collectors
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
        % Setup the grid for an electrode
        % Here, :code:`paramobj` is instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`
            
            paramobj.G = genSubGrid(gen.parentGrid, params.cellind);
            
        end

        function paramobj = setupActiveMaterialGrid(gen, paramobj, params)
        % Setup the grid for the active material
        % Here, :code:`paramobj` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
            
            paramobj.G = genSubGrid(gen.parentGrid, params.cellind);

            switch params.electrode_case
              case 'default'
                
                paramobj.G = genSubGrid(gen.parentGrid, params.cellind);
                
              case 'composite'
                
                G = genSubGrid(gen.parentGrid, params.cellind);
                paramobj.G          = G;
                paramobj.FirstMaterial.G = G;
                paramobj.SecondMaterial.G  = G;
                
              otherwise
                error('electrode_case not recognized')
            end

        end
        
        function paramobj = setupCurrentCollector(gen, paramobj, params)
        % Method that setups the grid and coupling terms for the current collectors
        % Here, :code:`paramobj` is instance of :class:`CurrentCollectorInputParams <Electrochemistry.CurrentCollectorInputParams>`

            paramobj = gen.setupCurrentCollectorGrid(paramobj, params);
            paramobj = gen.setupCurrentCollectorBcCoupTerm(paramobj, params);
            
        end

        function paramobj = setupCurrentCollectorGrid(gen, paramobj, params)
        % Setup a grid for a current collector
        % Here, :code:`paramobj` is instance of :class:`CurrentCollectorInputParams <Electrochemistry.CurrentCollectorInputParams>`

            paramobj.G = genSubGrid(gen.parentGrid, params.cellind);
            
        end       

        function paramobj = setupElectrodeElectrolyteCoupTerm(gen, paramobj, params)
        % Setup the coupling terms between the electrode and electrolyte
        % Here, :code:`paramobj` is instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`
            
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            am   = 'ActiveMaterial';
            
            couplingTerms = {};
            
            G_ne    = paramobj.(ne).(am).G;
            G_elyte = paramobj.(elyte).G;
            
            % parent Grid
            G = G_ne.parentGrid;
            
            % All the cells from NegativeElectrode are coupled with Electrolyte
            cells1 = (1 : G_ne.getNumberOfCells())';
            pcells = G_ne.mappings.cellmap(cells1);
            
            mapping = zeros(G.getNumberOfCells(), 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.getNumberOfCells())';
            cells2 = mapping(pcells);
            
            compnames = {'NegativeElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('NegativeElectrode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
            couplingTerms{end + 1} = coupTerm;
            
            G_pe = paramobj.(pe).(am).G;
            G_elyte = paramobj.(elyte).G;
            
            % All the cells from PositiveElectrode are coupled with Electrolyte
            cells1 = (1 : G_pe.getNumberOfCells())';
            pcells = G_pe.mappings.cellmap(cells1);
            
            mapping = zeros(G.getNumberOfCells(), 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.getNumberOfCells())';
            cells2 = mapping(pcells);
            
            compnames = {'PositiveElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('PositiveElectrode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
            couplingTerms{end + 1} = coupTerm;
            
            paramobj.couplingTerms = couplingTerms;

        end

        function paramobj = setupCurrentCollectorActiveMaterialCoupTerm(gen, paramobj, params)
        % Setup the coupling term between the current collector and the active material
        % Here, :code:`paramobj` is instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`

            am = 'ActiveMaterial';
            cc  = 'CurrentCollector';
            
            compnames = {'CurrentCollector', 'ActiveMaterial'};
            coupTerm = couplingTerm('CurrentCollector-ActiveMaterial', compnames);
            
            G_am = paramobj.(am).G;
            G_cc = paramobj.(cc).G;
            
            cctbl.faces = (1 : G_cc.topology.faces.num)';
            cctbl.globfaces = G_cc.mappings.facemap;
            cctbl = IndexArray(cctbl);
            
            eactbl.faces = (1 : G_am.topology.faces.num)';
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
            
            cc_coupcells = sum(G_cc.topology.faces.neighbors(cc_coupfaces, :), 2);
            am_coupcells = sum(G_am.topology.faces.neighbors(am_coupfaces, :), 2);
            
            coupTerm.couplingfaces = [cc_coupfaces, am_coupfaces];
            coupTerm.couplingcells = [cc_coupcells, am_coupcells];
            
            paramobj.couplingTerm = coupTerm;
            
        end

        function paramobj = setupCurrentCollectorBcCoupTerm(gen, paramobj, params)
        % Setup the boundary locations for the current collector
            
            % default setup
            compnames = {'CurrentCollector'};
            coupTerm = couplingTerm('bc-CurrentCollector', compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;
            
            paramobj.externalCouplingTerm = coupTerm;
            
        end

        function paramobj = setupActiveMaterialBcCoupTerm(gen, paramobj, params)
        % Setup the boundary locations for the active material (used in the absence of current collectors)
            
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
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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

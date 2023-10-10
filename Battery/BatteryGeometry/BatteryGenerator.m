classdef BatteryGenerator
% Base class that add grids and coupling terms to a paramobj instance of BatteryInputParams (through method 
% updateBatteryInputParams)
%
% This class goes through the whole grid setup and is meant to be used as a base class.
%
% Example : :class:`BatteryGenerator1D <BatteryGeometry.BatteryGenerator1D>`, :class:`BatteryGenerator2D <BatteryGeometry.BatteryGenerator2D>`, :class:`BatteryGenerator3D <BatteryGeometry.BatteryGenerator3D>`

    properties
        % Global grid 
        % It is stored here because is shared by many setup functions. 
        % It is constructed by the methode  :meth:`setupGrid`
        G
        
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

            params = pickField(params, 'Electrolyte');
            paramobj.Electrolyte = gen.setupElectrolyte(paramobj.Electrolyte, params);
            
            params = pickField(params, 'Separator');
            paramobj.Separator = gen.setupSeparator(paramobj.Separator, params);
            
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
            
            paramobj.ThermalModel.G = gen.G;
            coupTerm = couplingTerm('ThermalConvectiveCooling', {'ThermalModel'});
            coupTerm.couplingcells = params.couplingcells;
            coupTerm.couplingfaces = params.couplingfaces;
            paramobj.ThermalModel.couplingTerm = coupTerm;
        end
        
        function paramobj = setupElectrolyte(gen, paramobj, params)
        % Method that setups the grid and the coupling for the electrolyte model
        % Here, :code:`paramobj` is instance of :class:`ElectrolyteInputParams <Electrochemistry.ElectrolyteInputParams>`
            paramobj = gen.setupElectrolyteGrid(paramobj, params);
        end

        function paramobj = setupSeparator(gen, paramobj, params)
        % Method that setups the grid for the separator model
        % Here, :code:`paramobj` is instance of :class:`SeparatorInputParams <Electrochemistry.SeparatorInputParams>`            
            paramobj = gen.setupSeparatorGrid(paramobj, params);
        end

        
        function paramobj = setupElectrolyteGrid(gen, paramobj, params)
        % Setup the grid for the electrolyte
        % Here, :code:`paramobj` is instance of :class:`ElectrolyteInputParams <Electrochemistry.ElectrolyteInputParams>`
            
            % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
            
        end
        
        function paramobj = setupSeparatorGrid(gen, paramobj, params)
        % Setup the grid for the separator
        % Here, :code:`paramobj` is instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`
            
           % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
            
        end

        function paramobj = setupElectrodes(gen, paramobj, params)
        % Method that setups the grid and the coupling terms for both electrodes

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            eldes = {ne, pe};

            % We add the electrode type to the params structure (This information can be used by derived classes and may
            % then simplify setup)
            params.(ne).electrode_type = ne;
            params.(pe).electrode_type = pe;
            
            % setup Negative Electrode
            paramobj.(ne) = gen.setupElectrode(paramobj.(ne), params.(ne));
            % setup Positive Electrode
            paramobj.(pe) = gen.setupElectrode(paramobj.(pe), params.(pe));

            
        end
                
        function paramobj = setupElectrode(gen, paramobj, params)
        % Method that setups the grid and the coupling terms for an electrode
        % Here, :code:`paramobj` is instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`
            
            % shorthands 
            co = 'Coating';
            cc = 'CurrentCollector';
            
            % setup Electrode grid
            paramobj = gen.setupElectrodeGrid(paramobj, params);
            % setup Electrode coating component (co)
            paramobj.(co) = gen.setupCoatingGrid(paramobj.(co), params.(co));
            % setup current collector (cc)
            if paramobj.include_current_collectors
                % We add the electrode type to the params structure. (This information can be used by derived classes
                % and may then simplify setup)
                params.(cc).electrode_type = params.electrode_type;
                paramobj.(cc) = gen.setupCurrentCollector(paramobj.(cc), params.(cc));
                % setup coupling term between am and cc
                paramobj = gen.setupCurrentCollectorCoatingCoupTerm(paramobj, params);
            else
                paramobj.(co) = gen.setupCoatingBcCoupTerm(paramobj.(co), params.(co));
            end
            
        end
        
        function paramobj = setupElectrodeGrid(gen, paramobj, params)
        % Setup the grid for an electrode
        % Here, :code:`paramobj` is instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`
            
            paramobj.G = genSubGrid(gen.G, params.cellind);
        end

        function paramobj = setupCoatingGrid(gen, paramobj, params)
        % Setup the grid for the active material
        % Here, :code:`paramobj` is instance of :class:`CoatingInputParams <Electrochemistry.CoatingInputParams>`
            
            paramobj.G = genSubGrid(gen.G, params.cellind);

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

            paramobj.G = genSubGrid(gen.G, params.cellind);
            
        end       

        function paramobj = setupElectrodeElectrolyteCoupTerm(gen, paramobj, params)
        % Setup the coupling terms between the electrode and electrolyte
        % Here, :code:`paramobj` is instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`
            
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            co    = 'Coating';
            
            couplingTerms = {};
            
            G_ne = paramobj.(ne).(co).G;
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
            
            G_pe = paramobj.(pe).(co).G;
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

        function paramobj = setupCurrentCollectorCoatingCoupTerm(gen, paramobj, params)
        % Setup the coupling term between the current collector and the active material
        % Here, :code:`paramobj` is instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`

            co = 'Coating';
            cc = 'CurrentCollector';
            
            compnames = {'CurrentCollector', 'Coating'};
            coupTerm = couplingTerm('CurrentCollector-Coating', compnames);
            
            G_am = paramobj.(co).G;
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
        % Setup the boundary locations for the current collector
            
            % default setup
            compnames = {'CurrentCollector'};
            coupTerm = couplingTerm('bc-CurrentCollector', compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;
            
            paramobj.externalCouplingTerm = coupTerm;
            
        end

        function paramobj = setupCoatingBcCoupTerm(gen, paramobj, params)
        % Setup the boundary locations for the active material (used in the absence of current collectors)
            
            % default setup
            compname = 'Coating';
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

classdef BatteryGenerator
% Base class that add grids and coupling terms to a inputparams instance of BatteryInputParams (through method
% updateBatteryInputParams)
%
% This class goes through the whole grid setup and is meant to be used as a base class.
%
% Example : :class:`BatteryGeneratorP2D <BatteryGeometry.BatteryGeneratorP2D>`, :class:`BatteryGeneratorP3D <BatteryGeometry.BatteryGeneratorP3D>`, :class:`BatteryGeneratorP4D <BatteryGeometry.BatteryGeneratorP4D>`

    properties

        % Global grid
        % It is stored here because is shared by many setup functions.
        parentGrid

    end

    methods

        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams, params)
        %
        % This function is the main class function as it returns an updated :code:`inputparams` object with grid structure
        %

            error('virtual function');

        end


        function [inputparams, gen] = setupBatteryInputParams(gen, inputparams, params)
        % This is the main function which should be called by the derived class at the end of the :meth:`updateBatteryInputParams` method
        % The function set up the grid and the coupling terms and add those in the :code:`inputparams` structure which is
        % an instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`
        %

            [inputparams, gen] = gen.setupGrid(inputparams, params);

            params = pickField(params, 'Electrolyte');
            inputparams.Electrolyte = gen.setupElectrolyte(inputparams.Electrolyte, params);

            params = pickField(params, 'Separator');
            inputparams.Separator = gen.setupSeparator(inputparams.Separator, params);

            inputparams = gen.setupElectrodes(inputparams, params);

            if gen.use_thermal
                params_thermal = [];
                if isfield(params, 'ThermalModel')
                    params_thermal = params.ThermalModel;
                end
                inputparams = gen.setupThermalModel(inputparams, params_thermal);
            end

            inputparams = gen.setupElectrodeElectrolyteCoupTerm(inputparams);

        end


        function [inputparams, gen] = setupGrid(gen, inputparams, params)
        % setup :code:`inputparams.G` and update  :attr:`G`
        % Here, :code:`inputparams` is an instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`

            error('virtual function');

        end


        function inputparams = setupThermalModel(gen, inputparams, params)
        % Method that setups the grid and the coupling for the thermal model

            inputparams.ThermalModel.G = gen.parentGrid;
            coupTerm = couplingTerm('ThermalConvectiveCooling', {'ThermalModel'});
            coupTerm.couplingcells = params.couplingcells;
            coupTerm.couplingfaces = params.couplingfaces;
            inputparams.ThermalModel.couplingTerm = coupTerm;

        end


        function inputparams = setupElectrolyte(gen, inputparams, params)
        % Method that setups the grid and the coupling for the electrolyte model
        % Here, :code:`inputparams` is instance of :class:`ElectrolyteInputParams <Electrochemistry.ElectrolyteInputParams>`

            inputparams = gen.setupElectrolyteGrid(inputparams, params);

        end


        function inputparams = setupSeparator(gen, inputparams, params)
        % Method that setups the grid for the separator model
        % Here, :code:`inputparams` is instance of :class:`SeparatorInputParams <Electrochemistry.SeparatorInputParams>`

            inputparams = gen.setupSeparatorGrid(inputparams, params);

        end


        function inputparams = setupElectrolyteGrid(gen, inputparams, params)
        % Setup the grid for the electrolyte
        % Here, :code:`inputparams` is instance of :class:`ElectrolyteInputParams <Electrochemistry.ElectrolyteInputParams>`

            inputparams.G = genSubGrid(gen.parentGrid, params.cellind);

        end


        function inputparams = setupSeparatorGrid(gen, inputparams, params)
        % Setup the grid for the separator
        % Here, :code:`inputparams` is instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`

            inputparams.G = genSubGrid(gen.parentGrid, params.cellind);

        end


        function inputparams = setupElectrodes(gen, inputparams, params)
        % Method that setups the grid and the coupling terms for both electrodes

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            % We add the electrode type to the params structure (This information can be used by derived classes and may
            % then simplify setup)
            params.(ne).electrode_type = ne;
            params.(pe).electrode_type = pe;

            % setup Negative Electrode
            inputparams.(ne) = gen.setupElectrode(inputparams.(ne), params.(ne));
            % setup Positive Electrode
            inputparams.(pe) = gen.setupElectrode(inputparams.(pe), params.(pe));

        end


        function inputparams = setupElectrode(gen, inputparams, params)
        % Method that setups the grid and the coupling terms for an electrode
        % Here, :code:`inputparams` is instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`

            % shorthands
            co = 'Coating';
            cc = 'CurrentCollector';

            % setup Electrode grid
            inputparams = gen.setupElectrodeGrid(inputparams, params);

            % setup Electrode coating component (co)
            inputparams.(co) = gen.setupCoatingGrid(inputparams.(co), params.(co));

            % setup current collector (cc)
            if inputparams.include_current_collectors
                % We add the electrode type to the params
                % structure. (This information can be used by derived
                % classes and may then simplify setup)
                params.(cc).electrode_type = params.electrode_type;
                inputparams.(cc) = gen.setupCurrentCollector(inputparams.(cc), params.(cc));
                % setup coupling term between am and cc
                inputparams = gen.setupCurrentCollectorCoatingCoupTerm(inputparams, params);
            else
                inputparams.(co) = gen.setupCoatingBcCoupTerm(inputparams.(co), params.(co));
            end

        end


        function inputparams = setupElectrodeGrid(gen, inputparams, params)
        % Setup the grid for an electrode
        % Here, :code:`inputparams` is instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`

            inputparams.G = genSubGrid(gen.parentGrid, params.cellind);

        end


        function inputparams = setupCoatingGrid(gen, inputparams, params)
        % Setup the grid for the active material
        % Here, :code:`inputparams` is instance of :class:`CoatingInputParams <Electrochemistry.CoatingInputParams>`

            inputparams.G = genSubGrid(gen.parentGrid, params.cellind);

        end


        function inputparams = setupCurrentCollector(gen, inputparams, params)
        % Method that setups the grid and coupling terms for the current collectors
        % Here, :code:`inputparams` is instance of :class:`CurrentCollectorInputParams <Electrochemistry.CurrentCollectorInputParams>`

            inputparams = gen.setupCurrentCollectorGrid(inputparams, params);
            inputparams = gen.setupCurrentCollectorBcCoupTerm(inputparams, params);

        end


        function inputparams = setupCurrentCollectorGrid(gen, inputparams, params)
        % Setup a grid for a current collector
        % Here, :code:`inputparams` is instance of :class:`CurrentCollectorInputParams <Electrochemistry.CurrentCollectorInputParams>`

            inputparams.G = genSubGrid(gen.parentGrid, params.cellind);

        end


        function inputparams = setupElectrodeElectrolyteCoupTerm(gen, inputparams, ~)
        % Setup the coupling terms between the electrode and electrolyte
        % Here, :code:`inputparams` is instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            co    = 'Coating';

            couplingTerms = {};

            G_ne = inputparams.(ne).(co).G;
            G_elyte = inputparams.(elyte).G;

            % All the cells from NegativeElectrode are coupled with Electrolyte
            cells1 = (1 : G_ne.getNumberOfCells())';
            cells2 = G_elyte.mappings.invcellmap(G_ne.mappings.cellmap(cells1));

            compnames = {'NegativeElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('NegativeElectrode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty

            couplingTerms{end + 1} = coupTerm;

            G_pe = inputparams.(pe).(co).G;
            G_elyte = inputparams.(elyte).G;

            % All the cells from PositiveElectrode are coupled with Electrolyte
            cells1 = (1 : G_pe.getNumberOfCells())';
            cells2 = G_elyte.mappings.invcellmap(G_pe.mappings.cellmap(cells1));

            compnames = {'PositiveElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('PositiveElectrode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty

            couplingTerms{end + 1} = coupTerm;

            inputparams.couplingTerms = couplingTerms;

        end


        function inputparams = setupCurrentCollectorCoatingCoupTerm(gen, inputparams, ~)
        % Setup the coupling term between the current collector and the active material
        % Here, :code:`inputparams` is instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`

            co = 'Coating';
            cc = 'CurrentCollector';

            compnames = {'CurrentCollector', 'Coating'};
            coupTerm = couplingTerm('CurrentCollector-Coating', compnames);

            G_am = inputparams.(co).G;
            G_cc = inputparams.(cc).G;

            cctbl.faces = (1 : G_cc.getNumberOfFaces())';
            cctbl.globfaces = G_cc.mappings.facemap;
            cctbl = IndexArray(cctbl);

            eactbl.faces = (1 : G_am.getNumberOfFaces())';
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

            inputparams.couplingTerm = coupTerm;

        end


        function inputparams = setupCurrentCollectorBcCoupTerm(gen, inputparams, params)
        % Setup the boundary locations for the current collector

            compnames = {'CurrentCollector'};
            coupTerm = couplingTerm('bc-CurrentCollector', compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;

            inputparams.externalCouplingTerm = coupTerm;

        end


        function inputparams = setupCoatingBcCoupTerm(gen, inputparams, params)
        % Setup the boundary locations for the active material (used in the absence of current collectors)

            compname = 'Coating';
            compnames = {compname};
            coupname = sprintf('Exterior-%s', compname);
            coupTerm = couplingTerm(coupname, compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;

            inputparams.externalCouplingTerm = coupTerm;

        end

    end


end





%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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

classdef CoinCellBatteryGenerator < BatteryGenerator

    properties

        %
        % Design params for the components
        %
        % The format is a `matlab table structure <https://se.mathworks.com/help/matlab/tables.html>`_ where the rows are labeled with
        %
        % - NegativeCurrentCollector
        % - NegativeCoating
        % - Separator
        % - PositiveCoating
        % - PositiveCurrentCollector
        %
        % and the corresponding columns
        %
        % - thickness
        % - diameter
        % - numCellLayers : Discretization number in the layer
        %
        % see example :code:`runCR`
        compDims

        numRadial % Discretization number in the radial direction
        numAngular % Discretization number if the angular direction

        tag
        tagdict

        negativeExtCurrentFaces
        positiveExtCurrentFaces

        externalHeatTransferCoefficient = 1e3;

        include_current_collectors
        use_thermal

    end

    methods

        function gen = CoinCellBatteryGenerator()
            gen = gen@BatteryGenerator();
        end

        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams, params)

            gen.include_current_collectors = inputparams.include_current_collectors;
            gen.use_thermal = inputparams.use_thermal;

            fdnames = {'compDims', ...
                       'numRadial', ...
                       'numAngular'};
            gen = dispatchParams(gen, params, fdnames);

            inputparams = gen.setupBatteryInputParams(inputparams, []);

        end

        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            gen = coinCellGrid(gen);
            inputparams.G = gen.G;

        end

        function inputparams = setupElectrolyte(gen, inputparams, params)


            inds = cell2mat(values(gen.tagdict, {'PositiveCoating', 'Separator', 'NegativeCoating', 'Electrolyte'}));
            cellind = ismember(gen.tag, inds);
            params.cellind = find(cellind);

            inputparams = setupElectrolyte@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupSeparator(gen, inputparams, params)

            inds = gen.tagdict('Separator');
            cellind = ismember(gen.tag, inds);

            params.cellind = find(cellind);

            inputparams = setupSeparator@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupElectrodes(gen, inputparams, params)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';
            co = 'Coating';

            cellind = ismember(gen.tag, gen.tagdict('NegativeCoating'));
            params.(ne).(co).cellind = find(cellind);
            cellind = ismember(gen.tag, gen.tagdict('NegativeCurrentCollector'));
            params.(ne).(cc).cellind = find(cellind);
            params.(ne).(cc).extfaces = gen.negativeExtCurrentFaces;
            params.(ne).cellind = [params.(ne).(co).cellind; params.(ne).(cc).cellind];

            params.(ne).bcfaces = pi;

            cellind = ismember(gen.tag, gen.tagdict('PositiveCoating'));
            params.(pe).(co).cellind = find(cellind);
            cellind = ismember(gen.tag, gen.tagdict('PositiveCurrentCollector'));
            params.(pe).(cc).cellind = find(cellind);
            params.(pe).(cc).extfaces = gen.positiveExtCurrentFaces;
            params.(pe).cellind = [params.(pe).(co).cellind; params.(pe).(cc).cellind];

            params.(ne).bcfaces = exp(1);

            inputparams = setupElectrodes@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupCurrentCollectorBcCoupTerm(gen, inputparams, params)

            G = inputparams.G; % grid of the current collector
            extfaces = params.extfaces;

            globG = G.mappings.parentGrid;
            facemap = G.mappings.facemap;
            invfacemap = zeros(globG.faces.num, 1);
            invfacemap(facemap) = (1 : G.faces.num)';

            clear params
            params.bcfaces = invfacemap(extfaces);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

            inputparams = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupThermalModel(gen, inputparams, ~)

            [couplingfaces, couplingcells] = boundaryFaces(gen.G);

            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);

            coef = gen.externalHeatTransferCoefficient;
            inputparams.ThermalModel.externalHeatTransferCoefficient = coef;

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

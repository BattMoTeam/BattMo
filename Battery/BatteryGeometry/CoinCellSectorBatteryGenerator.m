classdef CoinCellSectorBatteryGenerator < BatteryGenerator

    properties

        % Design params
        thickness
        diameter
        angle
        offset
        numCellLayers
        nR

        tag     % cell-valued vector giving component number (indexing is given by tagdict)
        tagdict % dictionary giving the component number

        positiveExtCurrentFaces
        negativeExtCurrentFaces

        thermalCoolingFaces
        thermalExchangeFaces
        thermalExchangeFacesTag

        use_thermal

    end

    methods

        function gen = CoinCellSectorBatteryGenerator()
            gen = gen@BatteryGenerator();
        end

        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams, params)

            gen.thickness     = params.thickness;
            gen.diameter      = params.diameter;
            gen.angle         = params.angle;
            %gen.offset       = params.offset;
            gen.numCellLayers = params.numCellLayers;
            gen.nR            = params.nR;

            gen.use_thermal = inputparams.use_thermal;

            [inputparams, gen] = gen.setupBatteryInputParams(inputparams, []);

        end

        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            gen = coinCellSectorGrid(gen);
            inputparams.G = gen.G;

        end

        function inputparams = setupElectrolyte(gen, inputparams, params)

            inds = [gen.tagdict('PositiveActiveMaterial');
                    gen.tagdict('ElectrolyteSeparator');
                    gen.tagdict('NegativeActiveMaterial')];
            cellind = ismember(gen.tag, inds);
            params.cellind = find(cellind);
            inds = gen.tagdict('ElectrolyteSeparator');
            cellind = ismember(gen.tag, inds);
            params.Separator.cellind = find(cellind);

            inputparams = setupElectrolyte@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupElectrodes(gen, inputparams, params)

            % shorthands
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            am = 'ActiveMaterial';

            cellind = ismember(gen.tag, gen.tagdict('NegativeActiveMaterial'));
            params.(ne).(am).cellind = find(cellind);
            cellind = ismember(gen.tag, gen.tagdict('NegativeCurrentCollector'));
            params.(ne).(cc).cellind = find(cellind);
            params.(ne).(cc).extfaces = gen.negativeExtCurrentFaces;
            params.(ne).cellind = [params.(ne).(am).cellind; params.(ne).(cc).cellind];

            cellind = ismember(gen.tag, gen.tagdict('PositiveActiveMaterial'));
            params.(pe).(am).cellind = find(cellind);
            cellind = ismember(gen.tag, gen.tagdict('PositiveCurrentCollector'));
            params.(pe).(cc).cellind = find(cellind);
            params.(pe).(cc).extfaces = gen.positiveExtCurrentFaces;
            params.(pe).cellind = [params.(pe).(am).cellind; params.(pe).(cc).cellind];

            inputparams = setupElectrodes@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupCurrentCollectorBcCoupTerm(gen, inputparams, params)
        % inputparams is instance of CurrentCollectorInputParams

            G = inputparams.G; % grid of the current collector
            extfaces = params.extfaces;

            clear params
            globG = G.mappings.parentGrid;
            facemap = G.mappings.facemap;
            invfacemap = zeros(globG.faces.num, 1);
            invfacemap(facemap) = (1 : G.faces.num)';

            params.bcfaces = invfacemap(extfaces);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

            inputparams = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupThermalModel(gen, inputparams, ~)
        % inputparams is instance of BatteryInputParams

            G = gen.G;

            couplingfaces = gen.thermalCoolingFaces;
            couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);
            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);


            couplingfaces = gen.thermalExchangeFaces;
            couplingtags  = gen.thermalExchangeFacesTag;
            couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);

            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);

            thermal = 'ThermalModel'; % shorcut

            if isempty(inputparams.(thermal).externalHeatTransferCoefficientTopFaces) || ...
                    isempty(inputparams.(thermal).externalHeatTransferCoefficientSideFaces)
                inputparams.(thermal).externalHeatTransferCoefficient = ...
                    inputparams.(thermal).externalHeatTransferCoefficient*ones(numel(couplingfaces), 1);
            else
                externalHeatTransferCoefficient = nan(numel(couplingfaces), 1);
                externalHeatTransferCoefficient(couplingtags == 1) = inputparams.(thermal).externalHeatTransferCoefficientTopFaces;
                externalHeatTransferCoefficient(couplingtags == 2) = inputparams.(thermal).externalHeatTransferCoefficientSideFaces;
                inputparams.(thermal).externalHeatTransferCoefficient = externalHeatTransferCoefficient;
            end
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

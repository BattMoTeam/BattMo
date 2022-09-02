classdef CoinCellBatteryGenerator < BatteryGenerator

    properties

        % Design params
        compdims
        meshSize
        offset = [0, 0] % Negative and positive electrode offset half this distance from center
        
        % Sector model specific design params
        use_sector
        angle

        tag     % cell-valued vector giving component number (indexing is given by tagdict)
        tagdict % dictionary giving the component number

        negativeExtCurrentFaces
        positiveExtCurrentFaces

        % Heat parameters
        externalHeatTransferCoefficientTab = 1e3;
        externalHeatTransferCoefficient = 1e3;

        % For sector model
        thermalCoolingFaces
        thermalExchangeFaces = [];
        thermalExchangeFacesTag = [];

        use_thermal
    end

    methods

        function gen = CoinCellBatteryGenerator()
            gen = gen@BatteryGenerator();
        end

        function [paramobj, gen] = updateBatteryInputParams(gen, paramobj, params)

            gen.use_thermal = paramobj.use_thermal;

            gen.compdims = params.compdims;
            gen.meshSize = params.meshSize;
            gen.use_sector = params.use_sector;
            gen.angle = params.angle;
            %gen.offset = params.offset;
            assert(all(gen.offset == 0), 'Coin cells with non-zero offset is not yet implemented')
            
            [paramobj, gen] = gen.setupBatteryInputParams(paramobj, []);

        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)

            gen = coinCellGrid(gen);
            paramobj.G = gen.G;

        end

        function paramobj = setupElectrolyte(gen, paramobj, params)

            tagdict = gen.tagdict;
            tag = gen.tag;

            inds = cell2mat(values(tagdict, {'PositiveActiveMaterial', 'ElectrolyteSeparator', 'NegativeActiveMaterial', 'Electrolyte'}));
            cellind = ismember(tag, inds);
            params.cellind = find(cellind);

            inds = tagdict('ElectrolyteSeparator');
            cellind = ismember(tag, inds);
            params.Separator.cellind = find(cellind);

            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);

        end

        function paramobj = setupElectrodes(gen, paramobj, params)

            % shorthands
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            am = 'ActiveMaterial';

            tagdict = gen.tagdict;
            tag = gen.tag;

            cellind = ismember(tag, tagdict('NegativeActiveMaterial'));
            params.(ne).(am).cellind = find(cellind);
            cellind = ismember(tag, tagdict('NegativeCurrentCollector'));
            params.(ne).(cc).cellind = find(cellind);
            params.(ne).(cc).extfaces = gen.negativeExtCurrentFaces;
            params.(ne).cellind = [params.(ne).(am).cellind; params.(ne).(cc).cellind];

            cellind = ismember(tag, tagdict('PositiveActiveMaterial'));
            params.(pe).(am).cellind = find(cellind);
            cellind = ismember(tag, tagdict('PositiveCurrentCollector'));
            params.(pe).(cc).cellind = find(cellind);
            params.(pe).(cc).extfaces = gen.positiveExtCurrentFaces;
            params.(pe).cellind = [params.(pe).(am).cellind; params.(pe).(cc).cellind];

            paramobj = setupElectrodes@BatteryGenerator(gen, paramobj, params);

        end


        function paramobj = setupCurrentCollectorBcCoupTerm(gen, paramobj, params)
        % paramobj is instance of CurrentCollectorInputParams

            G = paramobj.G; % grid of the current collector
            extfaces = params.extfaces;

            globG = G.mappings.parentGrid;
            facemap = G.mappings.facemap;
            invfacemap = zeros(globG.faces.num, 1);
            invfacemap(facemap) = (1 : G.faces.num)';

            clear params
            params.bcfaces = invfacemap(extfaces);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

            paramobj = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, paramobj, params);

        end

        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams

        % Cooling on external faces

            [couplingfaces, couplingcells] = boundaryFaces(gen.G);

            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);

            coef = gen.externalHeatTransferCoefficient;
            paramobj.ThermalModel.externalHeatTransferCoefficient = coef;

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

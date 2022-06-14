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

        function [paramobj, gen] = updateBatteryInputParams(gen, paramobj, params)

            gen.thickness = params.thickness;
            gen.diameter  = params.diameter;
            gen.angle     = params.angle;
            %gen.offset    = params.offset;
            gen.numCellLayers    = params.numCellLayers;
            gen.nR        = params.nR;
            
            gen.use_thermal = paramobj.use_thermal;

            [paramobj, gen] = gen.setupBatteryInputParams(paramobj, []);

        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)

            gen = coinCellSectorGrid(gen);

            paramobj.G = gen.G;

        end

        function paramobj = setupElectrolyte(gen, paramobj, params)

            tagdict = gen.tagdict;
            tag = gen.tag;

            inds = [tagdict('PositiveActiveMaterial'); tagdict('ElectrolyteSeparator'); tagdict('NegativeActiveMaterial')];
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

            clear params
            globG = G.mappings.parentGrid;
            facemap = G.mappings.facemap;
            invfacemap = zeros(globG.faces.num, 1);
            invfacemap(facemap) = (1 : G.faces.num)';

            params.bcfaces = invfacemap(extfaces);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

            paramobj = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, paramobj, params);

        end

        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams

            G = gen.G;

            couplingfaces = gen.thermalCoolingFaces;
            couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);
            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);

            
            couplingfaces = gen.thermalExchangeFaces;
            couplingtags  = gen.thermalExchangeFacesTag;
            couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);

            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);

            thermal = 'ThermalModel'; % shorcut

            if isempty(paramobj.(thermal).externalHeatTransferCoefficientTopFaces) | ...
                    isempty(paramobj.(thermal).externalHeatTransferCoefficientSideFaces)
                paramobj.(thermal).externalHeatTransferCoefficient = ...
                    paramobj.(thermal).externalHeatTransferCoefficient*ones(numel(couplingfaces), 1);
            else
                externalHeatTransferCoefficient = nan(numel(couplingfaces), 1);
                externalHeatTransferCoefficient(couplingtags == 1) = paramobj.(thermal).externalHeatTransferCoefficientTopFaces;
                externalHeatTransferCoefficient(couplingtags == 2) = paramobj.(thermal).externalHeatTransferCoefficientSideFaces;
                paramobj.(thermal).externalHeatTransferCoefficient = externalHeatTransferCoefficient;
            end
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

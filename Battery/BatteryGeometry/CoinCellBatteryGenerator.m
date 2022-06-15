classdef CoinCellBatteryGenerator < BatteryGenerator

    properties

        % Design params
        thickness
        diameter
        offset
        numCellLayers
        meshSize
        
        tag     % cell-valued vector giving component number (indexing is given by tagdict)
        tagdict % dictionary giving the component number

        negativeExtCurrentFaces
        positiveExtCurrentFaces

        % Heat parameters
        externalHeatTransferCoefficientTab = 1e3;
        externalHeatTransferCoefficient = 1e3;

        use_thermal
    end

    methods

        function gen = CoinCellBatteryGenerator()
            gen = gen@BatteryGenerator();
        end

        function [paramobj, gen] = updateBatteryInputParams(gen, paramobj, params)

            paramobj.include_current_collectors = true;
            gen.use_thermal = paramobj.use_thermal;

            gen.thickness     = params.thickness;
            gen.diameter      = params.diameter;
            gen.offset        = params.offset;
            gen.numCellLayers = params.numCellLayers;
            gen.meshSize      = params.meshSize;
            
            [paramobj, gen] = gen.setupBatteryInputParams(paramobj, []);

        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)

            gen = coinCellGrid(gen);

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

            % G = gen.G;
            % figure, hold on
            % plotGrid(G, 'facecolor', 'none');
            % plotGrid(G, params.cellind, 'facecolor', 'b', 'facealpha', 0.2);
            % plotGrid(G, params.Separator.cellind, 'facecolor', 'r');
            % view(3)
            % keyboard;
            
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

            % G = gen.G;
            % figure, hold on
            % plotGrid(G, 'facecolor', 'none')
            % plotGrid(G, params.(ne).(am).cellind, 'facecolor', 'b')
            % plotGrid(G, params.(ne).(cc).cellind, 'facecolor', 'r')
            % view(3)

            % G = gen.G;
            % figure, hold on
            % plotGrid(G, 'facecolor', 'none')
            % plotFaces(G, params.(ne).(cc).extfaces)
            % view(3)
            % keyboard
            
            cellind = ismember(tag, tagdict('PositiveActiveMaterial'));
            params.(pe).(am).cellind = find(cellind);
            cellind = ismember(tag, tagdict('PositiveCurrentCollector'));
            params.(pe).(cc).cellind = find(cellind);
            params.(pe).(cc).extfaces = gen.positiveExtCurrentFaces;
            params.(pe).cellind = [params.(pe).(am).cellind; params.(pe).(cc).cellind];

            % G = gen.G;
            % figure, hold on
            % plotGrid(G, 'facecolor', 'none')
            % plotGrid(G, params.(pe).(am).cellind, 'facecolor', 'b')
            % plotGrid(G, params.(pe).(cc).cellind, 'facecolor', 'r')
            % view(3)
            
            % G = gen.G;
            % figure, hold on
            % plotGrid(G, 'facecolor', 'none')
            % plotFaces(G, params.(pe).(cc).extfaces)
            % view(3)
            % keyboard

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

            figure, hold on
            plotGrid(G, 'facecolor', 'none');
            plotFaces(G, params.bcfaces),view(3)
            figure, hold on
            plotGrid(G, 'facecolor', 'none');
            plotGrid(G, params.bccells),view(3)
            keyboard;
            
            paramobj = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, paramobj, params);

        end

        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams

            % G = gen.G;

            % couplingfaces = gen.thermalCoolingFaces;
            % couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);
            % params = struct('couplingfaces', couplingfaces, ...
            %                 'couplingcells', couplingcells);
            % paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);

            
            % couplingfaces = gen.thermalExchangeFaces;
            % couplingtags  = gen.thermalExchangeFacesTag;
            % couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);

            % params = struct('couplingfaces', couplingfaces, ...
            %                 'couplingcells', couplingcells);
            % paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);

            % thermal = 'ThermalModel'; % shorcut

            % if isempty(paramobj.(thermal).externalHeatTransferCoefficientTopFaces) | ...
            %         isempty(paramobj.(thermal).externalHeatTransferCoefficientSideFaces)
            %     paramobj.(thermal).externalHeatTransferCoefficient = ...
            %         paramobj.(thermal).externalHeatTransferCoefficient*ones(numel(couplingfaces), 1);
            % else
            %     externalHeatTransferCoefficient = nan(numel(couplingfaces), 1);
            %     externalHeatTransferCoefficient(couplingtags == 1) = paramobj.(thermal).externalHeatTransferCoefficientTopFaces;
            %     externalHeatTransferCoefficient(couplingtags == 2) = paramobj.(thermal).externalHeatTransferCoefficientSideFaces;
            %     paramobj.(thermal).externalHeatTransferCoefficient = externalHeatTransferCoefficient;
        % end

        % Cooling on external faces

            [couplingfaces, couplingcells] = boundaryFaces(gen.G);
            
            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);



            coef = gen.externalHeatTransferCoefficient;%*ones(bcfacetbl.num, 1);
            %coef(ind) = gen.externalHeatTransferCoefficientTab;
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

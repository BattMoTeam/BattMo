classdef CoinCellBatteryGenerator < BatteryGenerator

    properties

        %
        % Design params for the components
        %
        % The format is a `matlab table structure <https://se.mathworks.com/help/matlab/tables.html>`_ where the rows are labeled with
        %
        % - NegativeCurrentCollector
        % - NegativeActiveMaterial
        % - ElectrolyteSeparator
        % - PositiveActiveMaterial
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

        function [paramobj, gen] = updateBatteryInputParams(gen, paramobj, params)

            gen.include_current_collectors = paramobj.include_current_collectors;
            gen.use_thermal = paramobj.use_thermal;
            
            fdnames = {'compDims', ...
                       'numRadial', ...
                       'numAngular'};
            gen = dispatchParams(gen, params, fdnames);
            
            paramobj = gen.setupBatteryInputParams(paramobj, []);
            
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

            params.(ne).bcfaces = pi;
            
            cellind = ismember(tag, tagdict('PositiveActiveMaterial'));
            params.(pe).(am).cellind = find(cellind);
            cellind = ismember(tag, tagdict('PositiveCurrentCollector'));
            params.(pe).(cc).cellind = find(cellind);
            params.(pe).(cc).extfaces = gen.positiveExtCurrentFaces;
            params.(pe).cellind = [params.(pe).(am).cellind; params.(pe).(cc).cellind];

            params.(ne).bcfaces = exp(1);
            
            paramobj = setupElectrodes@BatteryGenerator(gen, paramobj, params);

        end


        function paramobj = setupCurrentCollectorBcCoupTerm(gen, paramobj, params)

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

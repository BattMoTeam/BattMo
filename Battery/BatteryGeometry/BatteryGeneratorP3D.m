classdef BatteryGeneratorP3D < BatteryGenerator
% Setup 2D grid

    properties


        % Vector of components lengths in x direction
        %
        % - x(1) : length of negative current collector (default: 10 micrometer)
        % - x(2) : length of negative active material (default: 100 micrometer)
        % - x(3) : length of separator (default: 50 micrometer)
        % - x(4) : length of positive active material (default: 80 micrometer)
        % - x(5) : length of positive current collector (default: 10 micrometer)
        %
        xlength = 1e-6*[10; 100; 50; 80; 10];

        ylength = 1e-2; % length in y direction (default: 1 cm)

        ccnenx = 10; % discretization number for negative current collector (default: 10)
        nenx   = 30; % discretization number for negative active material  (default: 30)
        sepnx  = 30; % discretization number for separator (default: 30)
        penx   = 30; % discretization number for positive active material  (default: 30)
        ccpenx = 10; % discretization number for positive current collector (default: 10)

        ny = 10; % discretization number in y direction (default: 10)

        include_current_collectors % boolean: true if grid for current collectors should be included

        use_thermal % boolean true if grid for thermal model should be setup.

        externalHeatTransferCoefficientTab = 1e3;  % heat transfer coefficient at tab boundary (default: 1e3)
        externalHeatTransferCoefficient = 1e3;     % heat transfer coefficient at boundary (default: 1e3)

    end

    methods

        function gen = BatteryGeneratorP3D()

            gen = gen@BatteryGenerator();

        end

        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams)

            fdnames = {'include_current_collectors', ...
                       'use_thermal'};

            gen = dispatchParams(gen, inputparams, fdnames);

            inputparams = gen.setupBatteryInputParams(inputparams, []);

        end

        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;

            ylength = gen.ylength;

            if gen.include_current_collectors
                ccnenx = gen.ccnenx;
                ccpenx = gen.ccpenx;
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
                xlength = gen.xlength;
            else
                nxs = [nenx; sepnx; penx];
                xlength = gen.xlength(2:4);
            end

            nx = sum(nxs);
            ny = gen.ny;

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = ylength/ny;
            y = rldecode(y, ny);
            y = [0; cumsum(y)];

            G = tensorGrid(x, y);

            parentGrid = Grid(G);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G  = G;
            gen.parentGrid = parentGrid;

        end

        function inputparams = setupElectrolyte(gen, inputparams, params)

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;

            if gen.include_current_collectors
                ccnenx = gen.ccnenx;
                ccpenx = gen.ccpenx;
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            else
                ccnenx = 0;
                ccpenx = 0;
                nxs = [nenx; sepnx; penx];
            end

            ny = gen.ny;
            nx = sum(nxs);

            istart = ccnenx + 1;
            ni = nenx + sepnx + penx;
            params.cellind = pickTensorCells(istart, ni, nx, ny);

            inputparams = setupElectrolyte@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupSeparator(gen, inputparams, params)

            if gen.include_current_collectors
                ccnenx = gen.ccnenx;
                ccpenx = gen.ccpenx;
                nxs = [ccnenx; gen.nenx; gen.sepnx; gen.penx; ccpenx];
            else
                ccnenx = 0;
                ccpenx = 0;
                nxs = [gen.nenx; gen.sepnx; gen.penx];
            end

            nx = sum(nxs);

            istart = ccnenx + gen.nenx + 1;
            ni = gen.sepnx;
            params.cellind = pickTensorCells(istart, ni, nx, gen.ny);

            inputparams = setupSeparator@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupElectrodes(gen, inputparams, params)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';
            co = 'Coating';

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;

            if gen.include_current_collectors
                ccnenx = gen.ccnenx;
                ccpenx = gen.ccpenx;
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            else
                ccnenx = 0;
                ccpenx = 0;
                nxs = [nenx; sepnx; penx];
           end

            nx = sum(nxs);
            ny = gen.ny;

            %% Negative electrode
            istart = 1;
            ni = ccnenx + nenx;
            params.(ne).cellind = pickTensorCells(istart, ni, nx, ny);

            if gen.include_current_collectors
                % Negative electrode - current collector
                istart = 1;
                ni = ccnenx;
                params.(ne).(cc).cellind = pickTensorCells(istart, ni, nx, ny);
            end

            % Negative electrode - electrode active component
            istart = ccnenx + 1;
            ni = nenx;
            params.(ne).(co).cellind = pickTensorCells(istart, ni, nx, ny);

            %% Positive electrode
            istart = ccnenx + nenx + sepnx + 1;
            ni = penx + ccpenx;
            params.(pe).cellind = pickTensorCells(istart, ni, nx, ny);

            if gen.include_current_collectors
                % Positive electrode - current collector
                istart = ccnenx + nenx + sepnx + penx + 1;
                ni = ccpenx;
                params.(pe).(cc).cellind = pickTensorCells(istart, ni, nx, ny);
            end

            % Positive electrode - electrode active component
            istart = ccnenx + nenx + sepnx + 1;
            ni = penx;
            params.(pe).(co).cellind = pickTensorCells(istart, ni, nx, ny);

            if ~gen.include_current_collectors
                % Save electrode type for convenience
                params.(ne).(co).electrode_type = ne;
                params.(pe).(co).electrode_type = pe;
            end

            inputparams = setupElectrodes@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupCurrentCollectorBcCoupTerm(gen, inputparams, params)

            params = gen.findBoundary(inputparams.G, params);
            inputparams = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupCoatingBcCoupTerm(gen, inputparams, params)

            params = gen.findBoundary(inputparams.G, params);
            inputparams = setupCoatingBcCoupTerm@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupThermalModel(gen, inputparams, ~)

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';

            if gen.include_current_collectors
                cc = 'CurrentCollector';
                coupterm_ne = inputparams.(ne).(cc).externalCouplingTerm;
                facemap_ne  = inputparams.(ne).(cc).G.mappings.facemap;
                coupterm_pe = inputparams.(pe).(cc).externalCouplingTerm;
                facemap_pe  = inputparams.(pe).(cc).G.mappings.facemap;
            else
                co = 'Coating';
                coupterm_ne = inputparams.(ne).(co).externalCouplingTerm;
                facemap_ne  = inputparams.(ne).(co).G.mappings.facemap;
                coupterm_pe = inputparams.(pe).(co).externalCouplingTerm;
                facemap_pe  = inputparams.(pe).(co).G.mappings.facemap;
            end

            % The cooling is done on the external faces
            pG = gen.parentGrid;
            extfaces = any(pG.topology.faces.neighbors == 0, 2);
            couplingfaces = find(extfaces);
            couplingcells = sum(pG.topology.faces.neighbors(couplingfaces, :), 2);

            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);

            % We assign the values different values to the tab
            tabtbl.faces = [facemap_ne(coupterm_ne.couplingfaces);
                            facemap_pe(coupterm_pe.couplingfaces)];
            tabtbl = IndexArray(tabtbl);
            bcfacetbl.faces = couplingfaces;
            bcfacetbl = IndexArray(bcfacetbl);

            map = TensorMap();
            map.fromTbl = bcfacetbl;
            map.toTbl = tabtbl;
            map.mergefds = {'faces'};
            ind = map.getDispatchInd();

            coef = gen.externalHeatTransferCoefficient*ones(bcfacetbl.num, 1);
            coef(ind) = gen.externalHeatTransferCoefficientTab;

            inputparams.ThermalModel.externalHeatTransferCoefficient = coef;

        end

    end


    methods (Access = private)

        function params = findBoundary(gen, G, params)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            fc = G.getFaceCentroids();
            xf = fc(:, 1);

            switch params.electrode_type
              case ne
                x0 = min(xf);
              case pe
                x0 = max(xf);
            end

            params.bcfaces = find(abs(xf - x0) < eps*1000);
            params.bccells = sum(G.topology.faces.neighbors(params.bcfaces, :), 2);

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

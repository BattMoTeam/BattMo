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

            fdnames = {'compDims', ...
                       'numRadial', ...
                       'numAngular', ...
                       'include_current_collectors', ...
                       'use_thermal'};
            gen = dispatchParams(gen, params, fdnames);

            inputparams = gen.setupBatteryInputParams(inputparams, []);

        end


        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

        % Components in order from z=0 (top) to z=zmax (bottom) with the
        % surrounding electrolyte last
            compNames = {'NegativeCurrentCollector', ...
                         'NegativeCoating', ...
                         'Separator', ...
                         'PositiveCoating', ...
                         'PositiveCurrentCollector', ...
                         'Electrolyte'};
            assert(numel(intersect(compNames, gen.compDims.Row)) == 5);

            thickness = gen.compDims.thickness;
            numCellLayers = gen.compDims.numCellLayers;
            dz = thickness ./ numCellLayers;

            % Repeat dz for each layer
            dz = rldecode(dz, numCellLayers);

            comptag = 1 : numel(compNames);
            gen.tagdict = containers.Map(compNames, comptag);

            gen.compDims = gen.compDims;

            % Extract basic dimensions
            R = 0.5 * max(gen.compDims.diameter);
            xc0 = [0, 0];

            % Radial points
            hr = R / gen.numRadial;
            diams = unique(gen.compDims.diameter);
            r = 0.5 * diams;
            dr = diff(r);
            x1 = linspace(0, r(1), ceil(r(1)/hr));

            % Make center cell smaller:
            x1 = [0 0.25*x1(2) x1(2:end)];

            % Make sure we have 3 nodes.
            %nr = max(3, ceil(dr/hr));
            nr = ceil(dr/hr);

            for k = 2:numel(r)
                x1 = [x1, linspace(r(k-1), r(k), nr(k-1))];
            end
            x1 = unique(x1)';
            x1(:,2) = zeros(size(x1,1), 1);

            % Curves for transfinite interpolation
            x1(1,:) = [];
            s = (x1(:,1)-x1(1,1))/(x1(end,1)-x1(1,1));
            f1 = x1;
            f2 = f1;
            t = linspace(0, 1, gen.numAngular+1)';
            c = [cos(t*2*pi), sin(t*2*pi)];
            g1 = f1(1,1)*c;
            g2 = f1(end,1)*c;

            G = transfiniteGrid(s, t, f1, f2, g1, g2, 'fill', true);

            % Extrude
            G = makeLayeredGrid(G, dz);
            G = computeGeometry(G);

            %% Tag
            thickness = gen.compDims.thickness;
            thickness_accum = [0; cumsum(thickness)];
            zc = G.cells.centroids(:, 3);
            gen.tag = repmat(gen.tagdict('Electrolyte'), G.cells.num, 1);

            for k = 1:numel(compNames)
                key = compNames{k};
                if ~strcmp(key, 'Electrolyte')
                    t0 = thickness_accum(k);
                    t1 = thickness_accum(k+1);
                    zidx = zc > t0 & zc < t1;

                    rc = vecnorm(G.cells.centroids(:, 1:2) - xc0, 2, 2);
                    ridx = rc < 0.5 * gen.compDims{key, 'diameter'};

                    idx = zidx & ridx;
                    gen.tag(idx) = gen.tagdict(key);
                end
            end
            assert(all(gen.tag > -1));

            %% Electrode faces where current is applied
            minz = min(G.nodes.coords(:, 3)) + 100*eps;
            maxz = max(G.nodes.coords(:, 3)) - 100*eps;
            fmin = G.faces.centroids(:, 3) <= minz;
            fmax = G.faces.centroids(:, 3) >= maxz;
            assert(sum(fmin) == sum(fmax));

            cn = find(gen.tag == gen.tagdict('NegativeCurrentCollector'), 1);
            cp = find(gen.tag == gen.tagdict('PositiveCurrentCollector'), 1);

            if abs(cn - minz) < abs(cp - minz)
                gen.negativeExtCurrentFaces = fmin;
                gen.positiveExtCurrentFaces = fmax;
            else
                gen.negativeExtCurrentFaces = fmax;
                gen.positiveExtCurrentFaces = fmin;
            end

            parentGrid = Grid(G);
            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            gen.parentGrid = parentGrid;

            inputparams.G = G;

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

            % Grid of the current collector
            G = inputparams.G;
            extfaces = params.extfaces;

            facemap = G.mappings.facemap;
            invfacemap = zeros(G.getNumberOfFaces(), 1);
            invfacemap(facemap) = (1 : G.getNumberOfFaces())';

            params.bcfaces = invfacemap(extfaces);
            params.bccells = sum(G.topology.faces.neighbors(params.bcfaces, :), 2);

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

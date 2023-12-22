classdef BatteryGeneratorP4D < BatteryGenerator
% Setup 3D grid with tab

    properties

        %
        % vector of x-lengths
        %
        % - x(1) : x-length of first tab (default: 4cm)
        % - x(2) : x-length between the tabs (default: 2cm)
        % - x(3) : x-length of last tab (default: 4cm)
        %
        xlength = 1e-2*[0.4; 0.2; 0.4];

        %
        % vector of y-lengths
        %
        % - y(1) : y-length of first tab (default : 1mm)
        % - y(2) : y-length between the tabs (default : 2cm)
        % - y(3) : y-length of last tab (default : 1mm)
        %
        ylength = 1e-2*[0.1; 2; 0.1];

        %
        % Vector of components lengths in z direction
        %
        % - z(1) : length of negative current collector (default: 10 micrometer)
        % - z(2) : length of negative active material (default: 100 micrometer)
        % - z(3) : length of separator (default: 50 micrometer)
        % - z(4) : length of positive active material (default: 80 micrometer)
        % - z(5) : length of positive current collector (default: 10 micrometer)
        %
        zlength = 1e-6*[10; 100; 50; 80; 10];

        % Shorthands used below
        % ne    : Negative electrode
        % pe    : Positive electrode
        % am    : Electrode active component
        % cc    : Current collector
        % elyte : Electrolyte


        facz = 1;

        sep_nz   = 3; % discretization number in z-direction for separator (default: 3)
        ne_am_nz = 3; % discretization number in z-direction for positive active material (default: 3)
        pe_am_nz = 3; % discretization number in z-direction for negative active material (default: 3)
        ne_cc_nz = 2; % discretization number in z-direction for negative current collector (default: 3)
        pe_cc_nz = 2; % discretization number in z-direction for positive current collector (default: 3)

        % Discretization resolution in x-direction

        facx = 1;

        int_elyte_nx = 3; % discretization number in x-direction interior region (default: 3)
        ne_cc_nx     = 3; % discretization number in x-direction negative tab region (default: 3)
        pe_cc_nx     = 3; % discretization number in x-direction positive tab region (default: 3)

        % Discretization resolution in y-direction

        facy = 1;

        elyte_ny = 4; % discretization number in y-direction interior region (default: 3)
        ne_cc_ny = 2; % discretization number in y-direction negative tab region (default: 3)
        pe_cc_ny = 2; % discretization number in y-direction positive tab region (default: 3)

        % Utility variables computed once and then shared by methods (should not be set)

        elyte_nz;
        allparams;
        invcellmap;

        externalHeatTransferCoefficientTab = 1e3; % Heat transfer coefficient at the tab
        externalHeatTransferCoefficient = 1e2;    % Heat transfer coefficient for the non-tab regions

        use_thermal

    end

    methods

        function gen = BatteryGeneratorP4D()
            gen = gen@BatteryGenerator();
        end

        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams)

            assert(inputparams.include_current_collectors, 'This geometry must include current collectors');
            gen.use_thermal = inputparams.use_thermal;
            inputparams = gen.setupBatteryInputParams(inputparams, []);

        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)

            % shorthands
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            cc    = 'CurrentCollector';
            elyte = 'Electrolyte';
            sep   = 'Separator';

            gen = gen.applyResolutionFactors();

            nxs = [gen.ne_cc_nx; gen.int_elyte_nx; gen.pe_cc_nx];
            nys = [gen.ne_cc_ny; gen.elyte_ny; gen.pe_cc_ny];
            nzs = [gen.ne_cc_nz; gen.ne_am_nz; gen.sep_nz; gen.ne_am_nz; gen.pe_cc_nz];

            x = gen.xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = gen.ylength./nys;
            y = rldecode(y, nys);
            y = [0; cumsum(y)];

            z = gen.zlength./nzs;
            z = rldecode(z, nzs);
            z = [0; cumsum(z)];

            G = tensorGrid(x, y, z);

            nx = sum(nxs);
            ny = sum(nys);
            nz = sum(nzs);

            dimGrid = [nx; ny; nz];

            gen.elyte_nz = gen.sep_nz + gen.ne_am_nz + gen.pe_am_nz;

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.elyte_nz];
            allparams.(elyte).cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + gen.ne_am_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.sep_nz];
            allparams.(sep).cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            %% setup gen.ne_eac

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.ne_am_nz];
            allparams.(ne).(co).cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            %% setup gen.pe_eac

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + gen.ne_am_nz + gen.sep_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.pe_am_nz];
            allparams.(pe).(co).cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            %% setup gen.ne_cc

            startSubGrid = [1; gen.ne_cc_ny + 1; 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.ne_cc_nz];
            allparams.(ne).(cc).cellind1 = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            % We add the tab

            startSubGrid = [1; 1; 1];
            dimSubGrid   = [gen.ne_cc_nx; gen.ne_cc_ny; gen.ne_cc_nz];
            allparams.(ne).(cc).cellindtab = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            allparams.(ne).(cc).cellind = [allparams.(ne).(cc).cellind1; allparams.(ne).(cc).cellindtab];

            %% setup gen.pe_cc

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + gen.elyte_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.pe_cc_nz];
            allparams.(pe).(cc).cellind1 = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            % We add the tab

            startSubGrid = [gen.ne_cc_nx + gen.int_elyte_nx + 1; gen.ne_cc_ny + gen.elyte_ny + 1; gen.ne_cc_nz + gen.elyte_nz + 1];
            dimSubGrid   = [gen.pe_cc_nx; gen.pe_cc_ny; gen.pe_cc_nz];
            allparams.(pe).(cc).cellindtab = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            allparams.(pe).(cc).cellind = [allparams.(pe).(cc).cellind1; allparams.(pe).(cc).cellindtab];

            cellind = [allparams.(elyte).cellind; allparams.(ne).(co).cellind; allparams.(pe).(co).cellind; allparams.(ne).(cc).cellind; allparams.(pe).(cc).cellind];

            rcellind = setdiff((1 : G.cells.num)', cellind);

            nGlob = G.cells.num;
            [G, cellmap, facemap, nodemap] = removeCells(G, rcellind);
            invcellmap = zeros(nGlob, 1);
            invcellmap(cellmap) = (1 : G.cells.num)';

            G = computeGeometry(G);

            inputparams.G = G;

            gen.invcellmap = invcellmap;
            gen.allparams = allparams;
            gen.G = G;

        end

        function gen = applyResolutionFactors(gen)

            facz = gen.facz;

            gen.sep_nz   = facz*gen.sep_nz;
            gen.ne_am_nz = facz*gen.ne_am_nz;
            gen.pe_am_nz = facz*gen.pe_am_nz;
            gen.ne_cc_nz = facz*gen.ne_cc_nz;
            gen.pe_cc_nz = facz*gen.pe_cc_nz;

            facx = gen.facx;

            gen.int_elyte_nx = facx*gen.int_elyte_nx;
            gen.ne_cc_nx     = facx*gen.ne_cc_nx;
            gen.pe_cc_nx     = facx*gen.pe_cc_nx;

            facy = gen.facy;

            gen.ne_cc_ny = facy*gen.ne_cc_ny;
            gen.pe_cc_ny = facy*gen.pe_cc_ny;
            gen.elyte_ny = facy*gen.elyte_ny;

        end

        function inputparams = setupElectrolyte(gen, inputparams, params)

            params = gen.allparams.Electrolyte;
            imap = gen.invcellmap;
            params.cellind = imap(params.cellind);

            inputparams = setupElectrolyte@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupSeparator(gen, inputparams, params)

            params = gen.allparams.Separator;
            imap = gen.invcellmap;
            params.cellind = imap(params.cellind);

            inputparams = setupSeparator@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupElectrodes(gen, inputparams, params)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';
            co = 'Coating';

            params = gen.allparams;
            imap = gen.invcellmap;

            params.(ne).(co).cellind = imap(params.(ne).(co).cellind);
            params.(ne).(cc).cellind = imap(params.(ne).(cc).cellind);
            params.(ne).(cc).name = 'negative';
            params.(ne).cellind = [params.(ne).(co).cellind; params.(ne).(cc).cellind];

            params.(pe).(co).cellind = imap(params.(pe).(co).cellind);
            params.(pe).(cc).cellind = imap(params.(pe).(cc).cellind);
            params.(pe).(cc).name = 'positive';
            params.(pe).cellind = [params.(pe).(co).cellind; params.(pe).(cc).cellind];

            inputparams = setupElectrodes@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupCurrentCollectorBcCoupTerm(gen, inputparams, params)

            G = inputparams.G;
            yf = G.faces.centroids(:, 2);

            switch params.name
              case 'negative'
                myf = min(yf);
              case 'positive'
                myf = max(yf);
            end

            params.bcfaces = find(abs(yf - myf) < eps*1000);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

            inputparams = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupThermalModel(gen, inputparams, params)

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            cc    = 'CurrentCollector';

            % the cooling is done on the external faces
            G = gen.G;
            extfaces = any(G.faces.neighbors == 0, 2);
            couplingfaces = find(extfaces);
            couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);

            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);

            tabcellinds = [gen.allparams.(pe).(cc).cellindtab; gen.allparams.(ne).(cc).cellindtab];
            tabtbl.cells = gen.invcellmap(tabcellinds);
            tabtbl = IndexArray(tabtbl);

            tbls = setupSimpleTables(G);
            cellfacetbl = tbls.cellfacetbl;

            tabcellfacetbl = crossIndexArray(tabtbl, cellfacetbl, {'cells'});
            tabfacetbl = projIndexArray(tabcellfacetbl, {'faces'});

            bcfacetbl.faces = couplingfaces;
            bcfacetbl = IndexArray(bcfacetbl);

            tabbcfacetbl = crossIndexArray(bcfacetbl, tabfacetbl, {'faces'});

            map = TensorMap();
            map.fromTbl = bcfacetbl;
            map.toTbl = tabbcfacetbl;
            map.mergefds = {'faces'};
            ind = map.getDispatchInd();

            coef = gen.externalHeatTransferCoefficient*ones(bcfacetbl.num, 1);
            coef(ind) = gen.externalHeatTransferCoefficientTab;

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

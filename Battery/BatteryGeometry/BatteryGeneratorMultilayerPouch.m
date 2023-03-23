classdef BatteryGeneratorMultilayerPouch < BatteryGenerator
% Setup 3D grid with tab

    properties

        % Physical dimension

        xlength = 1e-2*[0.4; 0.2; 0.4];
        ylength = 1e-2*[0.1; 2; 0.1];
        unit_cell_thickness = 1e-6*[10; 100; 50; 80; 10];
        n_layers = 5;
        zlength = [];%1e-6*[10; 100; 50; 80; 10];

        % Shorthands used below
        % ne    : Negative electrode
        % pe    : Positive electrode
        % am    : Electrode active component
        % cc    : Current collector
        % elyte : Electrolyte

        % Discretization resolution in z-direction

        facz = 1;

        sep_nz   = 3;
        ne_am_nz = 3;
        pe_am_nz = 3;
        ne_cc_nz = 2;
        pe_cc_nz = 2;

        % Discretization resolution in x-direction

        facx = 1;

        int_elyte_nx = 3;
        ne_cc_nx     = 3;
        pe_cc_nx     = 3;

        % Discretization resolution in y-direction

        facy = 1;

        ne_cc_ny = 2;
        pe_cc_ny = 2;
        elyte_ny = 4;

        % Utility variables computed once and then shared by methods (should not be set)
        elyte_nz;
        allparams;
        invcellmap;

        % Heat parameters
        externalHeatTransferCoefficientTab = 1e3;
        externalHeatTransferCoefficient = 1e3;

        use_thermal

    end

    methods

        function gen = BatteryGeneratorMultilayerPouch()

            gen = gen@BatteryGenerator();

            gen.zlength = gen.unit_cell_thickness;
            for ind = 2:gen.n_layers
                if rem(ind,2) == 0
                    gen.zlength = [gen.zlength; flipud(gen.unit_cell_thickness(1:end-1))];
                else
                    gen.zlength = [gen.zlength; gen.unit_cell_thickness(2:end)];
                end
            end

        end

        function [paramobj, gen] = updateBatteryInputParams(gen, paramobj)

            assert(paramobj.include_current_collectors, 'This geometry must include current collectors');
            gen.use_thermal = paramobj.use_thermal;
            paramobj = gen.setupBatteryInputParams(paramobj, []);

        end

        function [paramobj, gen] = setupGrid(gen, paramobj, ~)

            % shorthands
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            cc    = 'CurrentCollector';
            elyte = 'Electrolyte';
            sep   = 'Separator';

            gen = gen.applyResolutionFactors();

            nxs = [gen.ne_cc_nx; gen.int_elyte_nx; gen.pe_cc_nx];
            nys = [gen.ne_cc_ny; gen.elyte_ny; gen.pe_cc_ny];

            unit_cell_nzs     = [gen.ne_cc_nz; gen.ne_am_nz; gen.sep_nz;  gen.pe_am_nz; gen.pe_cc_nz];
            unit_cell_nzs_tag = {{ne,cc};      {ne,am};      {elyte,sep}; {pe,am};      {pe,cc}};
            nzs = unit_cell_nzs;
            nzs_tag = unit_cell_nzs_tag;

            for ind = 2:gen.n_layers
                if rem(ind,2) == 0
                    nzs = [nzs; flipud(unit_cell_nzs(1:end-1))];
                    nzs_tag = [nzs_tag; flipud(unit_cell_nzs_tag(1:end-1))];
                else
                    nzs = [nzs; unit_cell_nzs(2:end)];
                    nzs_tag = [nzs_tag; unit_cell_nzs_tag(2:end)];
                end
            end

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

            NZ = [0; cumsum(nzs)] + 1;

            % Initialize
            gen.allparams.(elyte).cellind = [];
            for k = 1:numel(nzs)
                gen.allparams = setfield(gen.allparams, nzs_tag{k}{:}, 'cellind', []);
            end

            for k = 2:numel(nzs)-1
                % TODO don't compute I, J
                % if k == 1
                %     % Account for tab
                %     [I, J, K] = ndgrid(1:G.cartDims(1)-1, 1:G.cartDims(2), NZ(k):NZ(k+1)-1);
                % elseif k == numel(nzs)
                %     % Account for tab
                %     [I, J, K] = ndgrid(2:G.cartDims(1), 1:G.cartDims(2), NZ(k):NZ(k+1)-1);
                % else
                %     [I, J, K] = ndgrid(2:G.cartDims(1)-1, 1:G.cartDims(2), NZ(k):NZ(k+1)-1);
                % end

                [I, J, K] = ndgrid(2:G.cartDims(1)-1, 1:G.cartDims(2), NZ(k):NZ(k+1)-1);
                c0 = getfield(gen.allparams, nzs_tag{k}{:}, 'cellind');
                c1 = sub2ind(G.cartDims, I(:), J(:), K(:));
                cellind = [c0; c1];

                % if k == 1
                %     cellind_bottom = c1;
                % elseif k == numel(nzs)
                %     cellind_top = c1;
                % end

                fprintf('%s %s\n', nzs_tag{k}{:})
                gen.allparams = setfield(gen.allparams, nzs_tag{k}{:}, 'cellind', cellind);

                if any(contains(nzs_tag{k}, am)) || any(contains(nzs_tag{k}, sep))
                    gen.allparams.(elyte).cellind = [gen.allparams.(elyte).cellind; cellind];
                end

            end
            gen.allparams.(elyte).cellind = unique(gen.allparams.(elyte).cellind);

            % Bottom tab
            assert(strcmp(nzs_tag{1}{2}, cc));
            [I, J, K] = ndgrid(1:G.cartDims(1)-1, 1:G.cartDims(2), NZ(1):NZ(2)-1);
            c0 = getfield(gen.allparams, nzs_tag{1}{:}, 'cellind');
            cellind_bottom = sub2ind(G.cartDims, I(:), J(:), K(:));
            gen.allparams = setfield(gen.allparams, nzs_tag{1}{:}, 'cellind', [c0; cellind_bottom]);

            [I, J, K] = ndgrid(1, 1:G.cartDims(2), NZ(1):NZ(2)-1);
            cellindtab_bottom = sub2ind(G.cartDims, I(:), J(:), K(:));
            gen.allparams = setfield(gen.allparams, nzs_tag{1}{:}, 'cellindtab', cellindtab_bottom);

            % Top tab
            assert(strcmp(nzs_tag{end}{2}, cc));
            [I, J, K] = ndgrid(2:G.cartDims(1), 1:G.cartDims(2), NZ(end-1):NZ(end)-1);
            c0 = getfield(gen.allparams, nzs_tag{end}{:}, 'cellind');
            cellind_top = sub2ind(G.cartDims, I(:), J(:), K(:));
            gen.allparams = setfield(gen.allparams, nzs_tag{end}{:}, 'cellind', [c0; cellind_top]);

            [I, J, K] = ndgrid(G.cartDims(1), 1:G.cartDims(2), NZ(end-1):NZ(end)-1);
            cellindtab_top = sub2ind(G.cartDims, I(:), J(:), K(:));
            gen.allparams = setfield(gen.allparams, nzs_tag{end}{:}, 'cellindtab', cellindtab_top);

            % figure,plotGrid(G),hold on
            plotGrid(G,gen.allparams.(elyte).cellind,'facecolor','b'),title(elyte)
            view(3)

            figure,plotGrid(G),hold on
            plotGrid(G,gen.allparams.(elyte).(sep).cellind,'facecolor','b'),title(sep)
            view(3)

            figure,plotGrid(G),hold on
            plotGrid(G,gen.allparams.(ne).(cc).cellind,'facecolor','b'),title(ne,cc)
            view(3)

            figure,plotGrid(G),hold on
            plotGrid(G,gen.allparams.(ne).(am).cellind,'facecolor','b'),title(ne,am)
            view(3)

            figure,plotGrid(G),hold on
            plotGrid(G,gen.allparams.(ne).(cc).cellindtab,'facecolor','b'),title(ne,'tab')
            view(3)

            figure,plotGrid(G),hold on
            plotGrid(G,gen.allparams.(pe).(cc).cellind,'facecolor','b'),title(pe,cc)
            view(3)

            figure,plotGrid(G),hold on
            plotGrid(G,gen.allparams.(pe).(am).cellind,'facecolor','b'),title(pe,am)
            view(3)

            figure,plotGrid(G),hold on
            plotGrid(G,gen.allparams.(pe).(cc).cellindtab,'facecolor','b'),title(pe,'tab')
            view(3)

            keyboard;


            cellind = [gen.allparams.(elyte).cellind; gen.allparams.(ne).(am).cellind; gen.allparams.(pe).(am).cellind; gen.allparams.(ne).(cc).cellind; gen.allparams.(pe).(cc).cellind];
            rcellind = setdiff((1 : G.cells.num)', cellind);

            figure,plotGrid(G),hold on, plotGrid(G,rcellind,'facecolor','r'),view(3)%,plotGrid(G,cellind_top,'facecolor', 'b'),view(3),plotGrid(G,cellind_bottom, 'facecolor', 'm'),view(3)
            keyboard;

            nGlob = G.cells.num;
            [G, cellmap] = removeCells(G, rcellind);
            gen.invcellmap = zeros(nGlob, 1);
            gen.invcellmap(cellmap) = (1 : G.cells.num)';

            G = computeGeometry(G);
            paramobj.G = G;
            gen.invcellmap = (1:G.cells.num)'; % TODO check if this is correct
            gen.G = G;

        end

        function gen = applyResolutionFactors(gen)

            gen.sep_nz   = gen.facz*gen.sep_nz;
            gen.ne_am_nz = gen.facz*gen.ne_am_nz;
            gen.pe_am_nz = gen.facz*gen.pe_am_nz;
            gen.ne_cc_nz = gen.facz*gen.ne_cc_nz;
            gen.pe_cc_nz = gen.facz*gen.pe_cc_nz;

            gen.int_elyte_nx = gen.facx*gen.int_elyte_nx;
            gen.ne_cc_nx     = gen.facx*gen.ne_cc_nx;
            gen.pe_cc_nx     = gen.facx*gen.pe_cc_nx;

            gen.ne_cc_ny = gen.facy*gen.ne_cc_ny;
            gen.pe_cc_ny = gen.facy*gen.pe_cc_ny;
            gen.elyte_ny = gen.facy*gen.elyte_ny;

        end

        function paramobj = setupElectrolyte(gen, paramobj, ~)

            params = gen.allparams.Electrolyte;
            imap = gen.invcellmap;
            params.cellind = imap(params.cellind);
            params.Separator.cellind = imap(params.Separator.cellind);

            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);

        end

        function paramobj = setupElectrodes(gen, paramobj, ~)


            % shorthands
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            am = 'ActiveMaterial';

            params = gen.allparams;
            imap = gen.invcellmap;

            params.(ne).(am).cellind = imap(params.(ne).(am).cellind);
            params.(ne).(cc).cellind = imap(params.(ne).(cc).cellind);
            params.(ne).(cc).name = 'negative';
            params.(ne).cellind = [params.(ne).(am).cellind; params.(ne).(cc).cellind];

            params.(pe).(am).cellind = imap(params.(pe).(am).cellind);
            params.(pe).(cc).cellind = imap(params.(pe).(cc).cellind);
            params.(pe).(cc).name = 'positive';
            params.(pe).cellind = [params.(pe).(am).cellind; params.(pe).(cc).cellind];

            paramobj = setupElectrodes@BatteryGenerator(gen, paramobj, params);

        end

        function paramobj = setupCurrentCollectorBcCoupTerm(gen, paramobj, params)

            G = paramobj.G;
            yf = G.faces.centroids(:, 2);

            switch params.name
              case 'negative'
                myf = min(yf);
              case 'positive'
                myf = max(yf);
            end

            params.bcfaces = find(abs(yf - myf) < eps*1000);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

            paramobj = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, paramobj, params);

        end

        function paramobj = setupThermalModel(gen, paramobj, ~)
        % paramobj is instance of BatteryInputParams
        %
        % We recover the external coupling terms for the current collectors

            % shorthands
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
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);

            tabcellinds = [gen.allparams.(pe).(cc).cellindtab; gen.allparams.(ne).(cc).cellindtab];
            tabtbl.cells = tabcellinds;
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

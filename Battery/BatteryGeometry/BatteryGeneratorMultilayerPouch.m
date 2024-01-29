classdef BatteryGeneratorMultilayerPouch < BatteryGenerator
% Setup 3D grid with tab

    properties

        % Physical dimensions (without tabs), default values are given
        pouch_width  = 100*milli*meter;
        pouch_height = 100*milli*meter;

        % For now: Tabs are placed in the center and have the same width, default values are given
        tab_width     = 50*milli*meter;
        ne_tab_height = 40*milli*meter;
        pe_tab_height = 20*milli*meter;

        % If cap_tabs=true, we cap the tabs. This simplifies the numerical simulation without changing the accuracy
        cap_tabs = false;

        % Layer thickness, default values are given
        unit_cell_thickness = 1e-6*[10; 100; 50; 80; 10];

        % Number of layers
        n_layers = 2;

        % Shorthands used below
        % ne    : Negative electrode
        % pe    : Positive electrode
        % am    : Electrode active component
        % cc    : Current collector
        % elyte : Electrolyte

        % Discretization resolution in z-direction

        facz = 1;

        sep_nz   = 3;
        ne_co_nz = 3;
        pe_co_nz = 3;
        ne_cc_nz = 3;
        pe_cc_nz = 3;

        % Discretization resolution in x-direction

        facx = 1;

        elyte_nx = 2;
        tab_nx = 3;

        % Discretization resolution in y-direction

        facy = 1;

        ne_tab_ny = 2;
        pe_tab_ny = 2;
        elyte_ny = 4;

        % Utility variables computed once and then shared by methods (should not be set)
        allparams;
        invcellmap;

        use_thermal

        % helpers
        tbls

    end

    methods

        function gen = BatteryGeneratorMultilayerPouch()

            gen = gen@BatteryGenerator();

        end

        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams)

            assert(inputparams.include_current_collectors, 'This geometry must include current collectors');
            gen.use_thermal = inputparams.use_thermal;
            [inputparams, gen] = gen.setupBatteryInputParams(inputparams, []);

        end

        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            % shorthands
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            cc    = 'CurrentCollector';
            elyte = 'Electrolyte';
            sep   = 'Separator';

            gen = gen.applyResolutionFactors();

            % Setup z
            zlength = gen.unit_cell_thickness;
            for ind = 2 : gen.n_layers
                if rem(ind,2) == 0
                    zlength = [zlength; flipud(gen.unit_cell_thickness(1:end-1))];
                else
                    zlength = [zlength; gen.unit_cell_thickness(2:end)];
                end
            end

            unit_cell_nzs     = [gen.ne_cc_nz; gen.ne_co_nz; gen.sep_nz;  gen.pe_co_nz; gen.pe_cc_nz];
            unit_cell_nzs_tag = {{ne,cc};      {ne,co};      {sep};       {pe,co};      {pe,cc}};
            nzs = unit_cell_nzs;
            nzs_tag = unit_cell_nzs_tag;

            for ind = 2 : gen.n_layers
                if rem(ind,2) == 0
                    nzs = [nzs; flipud(unit_cell_nzs(1:end-1))];
                    nzs_tag = [nzs_tag; flipud(unit_cell_nzs_tag(1:end-1))];
                else
                    nzs = [nzs; unit_cell_nzs(2:end)];
                    nzs_tag = [nzs_tag; unit_cell_nzs_tag(2:end)];
                end
            end

            z = zlength./nzs;
            z = rldecode(z, nzs);
            z = [0; cumsum(z)];

            % Setup widths
            x0 = 0.5*(gen.pouch_width - gen.tab_width);
            dxlength = [x0; gen.pouch_width-gen.tab_width; x0];
            nxs = [gen.elyte_nx; gen.tab_nx; gen.elyte_nx];
            x = dxlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            % Setup heights
            dylength = [gen.ne_tab_height; gen.pouch_height; gen.pe_tab_height];
            nys = [gen.ne_tab_ny; gen.elyte_ny; gen.pe_tab_ny];
            y = dylength./nys;
            y = rldecode(y, nys);
            y = [0; cumsum(y)];

            % Setup grid
            G = tensorGrid(x, y, z);

            % Setup index array for global grid with carthesian indices.
            nc = G.cells.num;
            [indx, indy, indz] = ind2sub([numel(x) - 1, numel(y) - 1, numel(z) - 1], (1 : nc)');
            globcelltbl.globcells = (1 : nc)';
            globcelltbl.indx      = indx;
            globcelltbl.indy      = indy;
            globcelltbl.indz      = indz;
            globcelltbl = IndexArray(globcelltbl);

            % assign dir index to faces
            % dir = 1 : perpendicular to x-axis
            % dir = 2 : perpendicular to y-axis
            % dir = 3 : perpendicular to z-axis

            tbls = setupSimpleTables(G);
            globfacetbl = tbls.facetbl;
            cellfacetbl = tbls.cellfacetbl;
            cellfacetbl = cellfacetbl.addInd('cartInd', G.cells.faces(:, 2));

            dir = zeros(globfacetbl.num, 1);
            cartinds = {{[1; 2], 1};
                        {[3; 4], 2};
                        {[5; 6], 3};
                       };
            for icart = 1 : numel(cartinds)
               cartind = cartinds{icart};
               clear findtbl
               findtbl.cartInd = cartind{1};
               findtbl = IndexArray(findtbl);
               cellfacecartindtbl = crossIndexArray(cellfacetbl, findtbl, {'cartInd'});
               selfacetbl = projIndexArray(cellfacecartindtbl, {'faces'});
               dir(selfacetbl.get('faces')) = cartind{2};
            end
            globfacetbl = globfacetbl.addInd('dir', dir);
            globfacetbl = replacefield(globfacetbl, {{'faces', 'globfaces'}});

            % Integer layers
            NZ = [0; cumsum(nzs)] + 1;

            % Initialize
            for k = 1:numel(nzs)
                gen.allparams = setfield(gen.allparams, nzs_tag{k}{:}, 'cellind', []);
            end
            gen.allparams.(ne).(cc).cellindtab = [];
            gen.allparams.(pe).(cc).cellindtab = [];

            for k = 1 : numel(nzs)

                % TODO don't compute I, J
                % Create interior slabs
                [I, J, K] = ndgrid(1:G.cartDims(1), (gen.ne_tab_ny+1):(G.cartDims(2)-gen.pe_tab_ny), NZ(k):NZ(k+1)-1);

                % Cells in this IJK box
                cbox = sub2ind(G.cartDims, I(:), J(:), K(:));

                % Cells previously marked with this tag
                cprev = getfield(gen.allparams, nzs_tag{k}{:}, 'cellind');

                % All cells, new and old
                cellind = [cbox; cprev];

                % Update
                gen.allparams = setfield(gen.allparams, nzs_tag{k}{:}, 'cellind', cellind);

                % Create tabs for the CCs
                create_tab = false;
                if strcmp(nzs_tag{k}{1}, ne) && strcmp(nzs_tag{k}{2}, cc)
                    % NE tab
                    [I, J, K] = ndgrid((gen.elyte_nx + 1) : (G.cartDims(1) - gen.elyte_nx), 1 : (gen.ne_tab_ny+1), NZ(k) : (NZ(k+1) - 1));
                    create_tab = ~gen.cap_tabs;
                elseif strcmp(nzs_tag{k}{1}, pe) && strcmp(nzs_tag{k}{2}, cc)
                    % PE tab
                    [I, J, K] = ndgrid((gen.elyte_nx+1):(G.cartDims(1)-gen.elyte_nx), (G.cartDims(2)-gen.pe_tab_ny-1):G.cartDims(2), NZ(k):NZ(k+1)-1);
                    create_tab = ~gen.cap_tabs;
                end

                if create_tab
                    cbox = sub2ind(G.cartDims, I(:), J(:), K(:));
                    cprev = getfield(gen.allparams, nzs_tag{k}{:}, 'cellindtab');
                    cellindtab = [cbox; cprev];
                    gen.allparams = setfield(gen.allparams, nzs_tag{k}{:}, 'cellindtab', cellindtab);
                end

            end

            % Electrolyte is the am's of pe and ne, as well as separator
            gen.allparams.(elyte).cellind = [gen.allparams.(ne).(co).cellind;
                                             gen.allparams.(sep).cellind;
                                             gen.allparams.(pe).(co).cellind];

            % CCs are including the tabs
            eldes = {ne, pe};
            for k = 1 : numel(eldes)
                elde = eldes{k};
                gen.allparams.(elde).(cc).cellind = [gen.allparams.(elde).(cc).cellind;
                                                     gen.allparams.(elde).(cc).cellindtab];
            end

            % Remove cells
            cellind = [gen.allparams.(sep).cellind;
                       gen.allparams.(ne).(co).cellind;
                       gen.allparams.(pe).(co).cellind;
                       gen.allparams.(ne).(cc).cellind;
                       gen.allparams.(pe).(cc).cellind];
            rcellind = setdiff((1 : G.cells.num)', cellind);
            nGlob = G.cells.num;
            [G, cellmap, facemap] = removeCells(G, rcellind);

            % Inverse map
            gen.invcellmap = zeros(nGlob, 1);
            gen.invcellmap(cellmap) = (1 : G.cells.num)';

            celltbl.globcells = cellmap;
            celltbl.cells = (1 : G.cells.num)';
            celltbl = IndexArray(celltbl);

            celltbl = crossIndexArray(celltbl, globcelltbl, {'globcells'});

            facetbl.globfaces = facemap;
            facetbl.faces = (1 : G.faces.num)';
            facetbl = IndexArray(facetbl);

            facedirtbl = crossIndexArray(facetbl, globfacetbl, {'globfaces'});
            facedirtbl = projIndexArray(facedirtbl, {'faces', 'dir'});

            % we setup extfacetbl. It contains the exterior faces orthogonal to y-axis.

            extfaces = find(any(G.faces.neighbors == 0, 2));
            extfacetbl.faces = extfaces;
            extfacetbl = IndexArray(extfacetbl);

            extfacedirtbl = crossIndexArray(extfacetbl, facedirtbl, {'faces'});

            seldirtbl.dir = 2;
            seldirtbl = IndexArray(seldirtbl);

            extfacetbl = crossIndexArray(extfacedirtbl, seldirtbl, {'dir'});
            extfacetbl = projIndexArray(extfacetbl, {'faces'});

            alltbls = setupSimpleTables(G);

            cellfacetbl = alltbls.cellfacetbl;

            tbls = struct('celltbl'    , celltbl    , ...
                          'extfacetbl' , extfacetbl , ...
                          'cellfacetbl', cellfacetbl, ...
                          'facedirtbl' , facedirtbl);

            parentGrid = Grid(G);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G = G;

            gen.parentGrid = parentGrid;
            gen.tbls       = tbls;

        end

        function gen = applyResolutionFactors(gen)

            gen.sep_nz   = gen.facz*gen.sep_nz;
            gen.ne_co_nz = gen.facz*gen.ne_co_nz;
            gen.pe_co_nz = gen.facz*gen.pe_co_nz;
            gen.ne_cc_nz = gen.facz*gen.ne_cc_nz;
            gen.pe_cc_nz = gen.facz*gen.pe_cc_nz;

            gen.elyte_nx = gen.facx*gen.elyte_nx;
            gen.tab_nx   = gen.facx*gen.tab_nx;

            gen.ne_tab_ny = gen.facy*gen.ne_tab_ny;
            gen.pe_tab_ny = gen.facy*gen.pe_tab_ny;
            gen.elyte_ny  = gen.facy*gen.elyte_ny;

        end

        function inputparams = setupElectrolyte(gen, inputparams, ~)

            params = gen.allparams.Electrolyte;
            imap = gen.invcellmap;
            params.cellind = imap(params.cellind);
            inputparams = setupElectrolyte@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupSeparator(gen, inputparams, ~)

            params = gen.allparams.Separator;
            imap = gen.invcellmap;
            params.cellind = imap(params.cellind);
            inputparams = setupSeparator@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupElectrodes(gen, inputparams, ~)

            % shorthands
            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';
            co = 'Coating';

            params  = gen.allparams;
            imap    = gen.invcellmap;
            tbls    = gen.tbls;

            celltbl = tbls.celltbl;

            params.(ne).(co).cellind = imap(params.(ne).(co).cellind);
            params.(ne).(cc).cellind = imap(params.(ne).(cc).cellind);
            params.(ne).(cc).name = ne;
            params.(ne).cellind = [params.(ne).(co).cellind; params.(ne).(cc).cellind];

            % setup current collector external coupling when cap_tabs is true.
            if gen.cap_tabs

                tabtbl.indx = gen.elyte_nx + (1 : gen.tab_nx)';
                tabtbl.indy = repmat(gen.ne_tab_ny + 1, numel(tabtbl.indx), 1);
                tabtbl = IndexArray(tabtbl);

                tabtbl = crossIndexArray(celltbl, tabtbl, {'indx', 'indy'});
                tabtbl = projIndexArray(tabtbl, {'cells'});

                cctbl.cells = params.(ne).(cc).cellind;
                cctbl = IndexArray(cctbl);

                tabtbl = crossIndexArray(tabtbl, cctbl, {'cells'});

                params.(ne).(cc).bccells = tabtbl.get('cells'); % indexing in gen.G (not in current collector grid)

            end

            params.(pe).(co).cellind = imap(params.(pe).(co).cellind);
            params.(pe).(cc).cellind = imap(params.(pe).(cc).cellind);
            params.(pe).(cc).name = pe;
            params.(pe).cellind = [params.(pe).(co).cellind; params.(pe).(cc).cellind];

            if gen.cap_tabs

                clear tabtbl;
                tabtbl.indx = gen.elyte_nx + (1 : gen.tab_nx)';
                tabtbl.indy = repmat(gen.ne_tab_ny + gen.elyte_ny, numel(tabtbl.indx), 1);
                tabtbl = IndexArray(tabtbl);

                tabtbl = crossIndexArray(celltbl, tabtbl, {'indx', 'indy'});
                tabtbl = projIndexArray(tabtbl, {'cells'});

                clear cctbl;
                cctbl.cells = params.(pe).(cc).cellind;
                cctbl = IndexArray(cctbl);

                tabtbl = crossIndexArray(tabtbl, cctbl, {'cells'});

                params.(pe).(cc).bccells = tabtbl.get('cells'); % indexing in gen.G (not in current collector grid)

            end

            inputparams = setupElectrodes@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupCurrentCollectorBcCoupTerm(gen, inputparams, params)

            G = inputparams.G;

            if gen.cap_tabs

                tbls = setupSimpleTables(G.mrstFormat());
                celltbl     = tbls.celltbl;
                facetbl     = tbls.facetbl;
                cellfacetbl = tbls.cellfacetbl;

                celltbl = celltbl.addInd('globcells', inputparams.G.mappings.cellmap);
                facetbl = facetbl.addInd('globfaces', inputparams.G.mappings.facemap);

                cellfacetbl = crossIndexArray(cellfacetbl, facetbl, {'faces'});
                cellfacetbl = crossIndexArray(cellfacetbl, celltbl, {'cells'});

                bccelltbl.globcells = params.bccells;
                bccelltbl = IndexArray(bccelltbl);

                tblgen = CrossIndexArrayGenerator();
                tblgen.tbl1 = bccelltbl;
                tblgen.tbl2 = gen.tbls.cellfacetbl;
                tblgen.replacefds2 = {{'cells', 'globcells'}, {'faces', 'globfaces'}};
                tblgen.mergefds = {'globcells'};

                bccellfacetbl = tblgen.eval();

                tblgen = CrossIndexArrayGenerator();
                tblgen.tbl1 = bccellfacetbl;
                tblgen.tbl2 = gen.tbls.extfacetbl;
                tblgen.replacefds2 = {{'faces', 'globfaces'}};
                tblgen.mergefds = {'globfaces'};

                bccellfacetbl = tblgen.eval();

                bccellfacetbl = crossIndexArray(bccellfacetbl, cellfacetbl, {'globcells', 'globfaces'});

                params.bcfaces = bccellfacetbl.get('faces');
                params.bccells = bccellfacetbl.get('cells');

            else

                fc = G.getFaceCentroids();
                yf = fc(:, 2);

                switch params.name
                  case 'NegativeElectrode'
                    myf = min(yf);
                  case 'PositiveElectrode'
                    myf = max(yf);
                end

                params.bcfaces = find(abs(yf - myf) < eps*1000);
                params.bccells = sum(G.topology.faces.neighbors(params.bcfaces, :), 2);

            end

            inputparams = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, inputparams, params);

        end

        function inputparams = setupThermalModel(gen, inputparams, ~)
        % inputparams is instance of BatteryInputParams
        %

            % shorthands
            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';

            if gen.cap_tabs

                tbls = gen.tbls;
                facedirtbl  = tbls.facedirtbl;
                cellfacetbl = tbls.cellfacetbl;

                G = gen.parentGrid;
                extfaces = any(G.topology.faces.neighbors == 0, 2);
                extfacetbl.faces = find(extfaces);
                extfacetbl = IndexArray(extfacetbl);

                seldirtbl.dir = [1; 2];
                seldirtbl = IndexArray(seldirtbl);

                bcfacetbl = crossIndexArray(facedirtbl, seldirtbl, {'dir'});
                bcfacetbl = crossIndexArray(bcfacetbl, extfacetbl, {'faces'});

                bccellfacetbl = crossIndexArray(cellfacetbl, bcfacetbl, {'faces'});

                couplingfaces = bccellfacetbl.get('faces');
                couplingcells = bccellfacetbl.get('cells');

                params = struct('couplingfaces', couplingfaces, ...
                                'couplingcells', couplingcells);
                inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);

            else

                % We use the tab
                eldes = {ne, pe};

                couplingfaces = [];
                couplingcells = [];

                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    cc_couplingfaces = inputparams.(elde).(cc).externalCouplingTerm.couplingfaces;
                    cc_couplingcells = inputparams.(elde).(cc).externalCouplingTerm.couplingcells;

                    mappings = inputparams.G.mappings;
                    couplingfaces = [couplingfaces; mappings.facemap(cc_couplingfaces)];
                    couplingcells = [couplingcells; mappings.cellmap(cc_couplingcells)];
                end

                params = struct('couplingfaces', couplingfaces, ...
                                'couplingcells', couplingcells);
                inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);

            end
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

classdef SpiralBatteryGenerator < BatteryGenerator
% Setup a grid Jelly Roll model

    properties

        nwindings % number of windings in the spiral

        rInner    % Inner Radius correspoding to the empty space in the middle

        %
        % Dictionary of widths for each component. The required key names for the dictionary are
        %
        % - 'Separator'
        % - 'NegativeCoating'
        % - 'NegativeCurrentCollector'
        % - 'PositiveCoating'
        % - 'PositiveCurrentCollector'
        widthDict

        % dicionary with number of cell in radial direction for each component (same keys as in widthDict).
        nrDict

        L            % length of the battery
        nas          % number of cells in the angular direction
        nL           % number of discretization cells in the longitudonal
        refLcoef     % coefficient use in refinement at top/bottom
        angleuniform

        exteriorNegativeElectrodeLayer % boolean, true if we want to have negative electrode on all outer layers.
        
        tag     % cell-valued vector giving component number (indexing is given by tagdict)
        tagdict % dictionary giving the component number

        tabparams % parameters for the tab on the positive current collector

        positiveExtCurrentFaces
        negativeExtCurrentFaces
        thermalExchangeFaces
        thermalExchangeFacesTag

        celltbl
        widthLayer
        nWidthLayer
        heightLayer
        nHeightLayer

        % computed tab width (due to discretization, we cannot enforce the tab widths)

        tabwidths

        % for the tabs (implemented only for aligned tabs now)

        windingnumbers

        % added only if exteriorNegativeElectrodeLayer
        full_G
        full_celltbl

        %
        use_thermal
        
    end


    methods

        function gen = SpiralBatteryGenerator()

            gen = gen@BatteryGenerator();

        end


        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams, params)

            gen.nwindings    = params.nwindings;
            gen.rInner       = params.rInner;
            gen.widthDict    = params.widthDict;
            gen.nrDict       = params.nrDict;
            gen.nas          = params.nas;
            gen.L            = params.L;
            gen.nL           = params.nL;
            gen.refLcoef     = params.refLcoef;
            gen.tabparams    = params.tabparams;
            gen.angleuniform = params.angleuniform;
            
            gen.exteriorNegativeElectrodeLayer = params.exteriorNegativeElectrodeLayer;
            
            gen.use_thermal = inputparams.use_thermal;

            [inputparams, gen] = gen.setupBatteryInputParams(inputparams, []);

        end


        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            output = spiralGrid(gen);

            G = output.G;

            % Transform to differentiable Grid structure.
            parentGrid = Grid(G);
            G          = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            gen.parentGrid = parentGrid;

            % Assign properties
            gen.tag                     = output.tag;
            gen.tagdict                 = output.tagdict;
            gen.positiveExtCurrentFaces = output.positiveExtCurrentFaces;
            gen.negativeExtCurrentFaces = output.negativeExtCurrentFaces;
            gen.thermalExchangeFaces    = output.thermalExchangeFaces;
            gen.thermalExchangeFacesTag = output.thermalExchangeFacesTag;
            gen.widthLayer              = output.widthLayer;
            gen.nWidthLayer             = output.nWidthLayer;
            gen.heightLayer             = output.heightLayer;
            gen.nHeightLayer            = output.nHeightLayer;
            gen.celltbl                 = output.celltbl;
            if gen.exteriorNegativeElectrodeLayer
                gen.full_G       = output.full_G;
                gen.full_celltbl = output.full_celltbl;
            end
         
            inputparams.G = G;

        end


        function UGrids = setupUnRolledGrids(gen, inputparams, varargin)

        % Allow for setting the size of the unrolled grid
            opt = struct('setX', [], ...
                         'setY', [], ...
                         'setZ', []);
            opt = merge_options(opt, varargin{:});

            if gen.exteriorNegativeElectrodeLayer
                G = Grid(gen.full_G);
            else
                G = inputparams.G;
            end

            widthLayer   = gen.widthLayer;
            nWidthLayer  = gen.nWidthLayer;
            heightLayer  = gen.heightLayer;
            nHeightLayer = gen.nHeightLayer;
            nas          = gen.nas;
            nL           = gen.nL;
            nwindings    = gen.nwindings;

            cartG = cartGrid([nas*nwindings, sum(nWidthLayer), nL]); % NB: will be Grid() later

            vecttbl.vect = (1 : cartG.griddim)';
            vecttbl = IndexArray(vecttbl)';

            [indi, indj, indk] = ind2sub([nas*nwindings, sum(nWidthLayer), nL], (1 : cartG.cells.num)');
            cartcelltbl.indi = indi;
            cartcelltbl.indj = indj;
            cartcelltbl.indk = indk;
            cartcelltbl.cells = (1 : cartG.cells.num)';
            cartcelltbl = IndexArray(cartcelltbl);

            cellindjtbl.indj = (1 : sum(nWidthLayer))';
            cellindjtbl = IndexArray(cellindjtbl);

            map = TensorMap();
            map.fromTbl = cellindjtbl;
            map.toTbl = cartcelltbl;
            map.mergefds = {'indj'};
            map = map.setup();

            w = map.eval(widthLayer);

            cellindktbl.indk = (1 : sum(nHeightLayer))';
            cellindktbl = IndexArray(cellindktbl);

            map = TensorMap();
            map.fromTbl = cellindktbl;
            map.toTbl = cartcelltbl;
            map.mergefds = {'indk'};
            map = map.setup();

            h = map.eval(heightLayer);

            vol = G.getVolumes();

            if gen.exteriorNegativeElectrodeLayer
                % the grid is created as if it was a full grid (without removing some positive electrode)
                celltbl = gen.full_celltbl.removeInd({'cells', 'indi', 'indj'});
            else
                celltbl = gen.celltbl.removeInd({'cells', 'indi', 'indj'});
            end
            
            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = cartcelltbl;
            map.replaceFromTblfds = {{'curvindi', 'indi'}, {'curvindj', 'indj'}};
            map.mergefds = {'indi', 'indj', 'indk'};

            cartInd = map.getDispatchInd();

            map = map.setup();

            vol = map.eval(vol);
            l = vol./(h.*w);

            % Here, we assume knowledge on ordering of cartcelltbl;
            l = l(1 : nas*nwindings);

            % we add cartesian indexing to the nodes
            [indi, indj, indk] = ind2sub([nas*nwindings + 1, sum(nWidthLayer) + 1, nL + 1], (1 : cartG.nodes.num)');

            nodetbl.nodes = (1 : cartG.nodes.num)';
            nodetbl.indi = indi;
            nodetbl.indj = indj;
            nodetbl.indk = indk;
            nodetbl = IndexArray(nodetbl);

            xnode = [0; cumsum(l)];
            ynode = [0; cumsum(widthLayer)];
            znode = [0; cumsum(heightLayer)];

            nodeinditbl.indi = (1 : (nas*nwindings + 1))';
            nodeinditbl = IndexArray(nodeinditbl);

            map = TensorMap();
            map.fromTbl = nodeinditbl;
            map.toTbl = nodetbl;
            map.mergefds = {'indi'};
            map = map.setup();

            xnode = map.eval(xnode);

            nodeindjtbl.indj = (1 : (sum(nWidthLayer) + 1))';
            nodeindjtbl = IndexArray(nodeindjtbl);

            map = TensorMap();
            map.fromTbl = nodeindjtbl;
            map.toTbl = nodetbl;
            map.mergefds = {'indj'};
            map = map.setup();

            ynode = map.eval(ynode);

            nodeindktbl.indk = (1 : (nL + 1))';
            nodeindktbl = IndexArray(nodeindktbl);

            map = TensorMap();
            map.fromTbl = nodeindktbl;
            map.toTbl = nodetbl;
            map.mergefds = {'indk'};
            map = map.setup();

            znode = map.eval(znode);

            if ~isempty(opt.setX)
                assert(abs(min(xnode)) < eps);
                assert(abs(max(xnode)) > eps);
                xnode = xnode / max(xnode) * opt.setX;
            end
            if ~isempty(opt.setY)
                error('Not implemented');
            end
            if ~isempty(opt.setZ)
                error('Not implemented');
            end

            cartG.nodes.coords = [xnode, ynode, znode];

            % In the case of exteriorNegativeElectrodeLayer, we need to remove part of the grid

            if gen.exteriorNegativeElectrodeLayer

                celltbl = gen.celltbl.removeInd({'cells', 'indi', 'indj'});               

                celltbl = replacefield(celltbl, {{'curvindi', 'indi'}, {'curvindj', 'indj'}});
                [cartcelltbl, indstruct] = crossIndexArray(celltbl, cartcelltbl, {'indi', 'indj', 'indk'});

                cartcelltbl = cartcelltbl.addInd('roll_cells', indstruct{1}.inds);
                cartcelltbl = sortIndexArray(cartcelltbl, {'roll_cells', 'cells'});
                
                [cartG, cellmap, facemap, nodemap] = extractSubgrid(cartG, cartcelltbl.get('cells'));

                converttbl.newcells = (1 : cartG.cells.num)';
                converttbl.cells = cellmap;
                converttbl = IndexArray(converttbl);

                map = TensorMap();
                map.fromTbl  = cartcelltbl;
                map.toTbl    = converttbl;
                map.mergefds = {'cells'};
                
                cartInd = map.getDispatchInd();
                
            end

            cartG = computeGeometry(cartG);

            %% we setup the mappings between the different grids

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            am    = 'ActiveMaterial';
            cc    = 'CurrentCollector';
            elyte = 'Electrolyte';
            sep   = 'Separator';

            tagdict = gen.tagdict;
            tag     = gen.tag;

            clear celltbl
            celltbl.cartcells = (1 : cartG.cells.num)';
            celltbl.cells = cartInd;
            celltbl = IndexArray(celltbl);

            UGrids.G = cartG;
            UGrids.G.mappings.ind = cartInd;

            % Convert to Grid()
            cartG = Grid(cartG);

            compG = inputparams.(pe).(cc).G;
            comptag = tagdict('PositiveCurrentCollector');
            compCartG = genSubGrid(cartG, find(tag(cartInd) == comptag));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            Gs.CurrentCollector = compCartG;

            compG = inputparams.(pe).(co).G;
            comptag = tagdict('PositiveCoating');
            compCartG = genSubGrid(cartG, find(tag(cartInd) == comptag));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            Gs.Coating = compCartG;

            UGrids.PositiveElectrode = Gs;
            clear Gs;

            compG = inputparams.(ne).(cc).G;
            comptag = tagdict('NegativeCurrentCollector');
            compCartG = genSubGrid(cartG, find(tag(cartInd) == comptag));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            Gs.CurrentCollector = compCartG;

            compG = inputparams.(ne).(co).G;
            comptag = tagdict('NegativeCoating');
            compCartG = genSubGrid(cartG, find(tag(cartInd) == comptag));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            Gs.Coating = compCartG;

            UGrids.NegativeElectrode = Gs;
            clear Gs;

            compG = inputparams.(elyte).G;
            comptags = [tagdict('Separator'); tagdict('NegativeCoating'); tagdict('PositiveCoating')];
            compCartG = genSubGrid(cartG, find(ismember(tag(cartInd), comptags)));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            UGrids.Electrolyte = compCartG;

            UGrids.ThermalModel = UGrids.G;

        end


        function compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl)

            compcelltbl.loccells = (1 : compG.getNumberOfCells())';
            compcelltbl.cells = compG.mappings.cellmap;
            compcelltbl = IndexArray(compcelltbl);
            compcelltbl = crossIndexArray(compcelltbl, celltbl, {'cells'});
            compcelltbl = compcelltbl.removeInd({'cells'});

            compcartcelltbl.cartloccells = (1 : compCartG.getNumberOfCells())';
            compcartcelltbl.cartcells = compCartG.mappings.cellmap;
            compcartcelltbl = IndexArray(compcartcelltbl);

            map = TensorMap();
            map.fromTbl = compcelltbl;
            map.toTbl = compcartcelltbl;
            map.mergefds = {'cartcells'};

            ind = map.getDispatchInd();

            compCartG.mappings.ind = ind;

        end


        function inputparams = setupElectrolyte(gen, inputparams, params)

            tagdict = gen.tagdict;
            tag = gen.tag;

            inds = [tagdict('PositiveCoating'); tagdict('Separator'); tagdict('NegativeCoating')];
            cellind = ismember(tag, inds);
            params.cellind = find(cellind);

            inputparams = setupElectrolyte@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupSeparator(gen, inputparams, params)

            tagdict = gen.tagdict;
            tag = gen.tag;

            inds = tagdict('Separator');
            cellind = ismember(tag, inds);

            params.cellind = find(cellind);

            inputparams = setupSeparator@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupElectrodes(gen, inputparams, params)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            co  = 'Coating';

            tagdict = gen.tagdict;
            tag = gen.tag;

            cellind = ismember(tag, tagdict('NegativeCoating'));
            params.(ne).(co).cellind = find(cellind);
            cellind = ismember(tag, tagdict('NegativeCurrentCollector'));
            params.(ne).(cc).cellind = find(cellind);
            params.(ne).(cc).extfaces = gen.negativeExtCurrentFaces;
            params.(ne).cellind = [params.(ne).(co).cellind; params.(ne).(cc).cellind];

            cellind = ismember(tag, tagdict('PositiveCoating'));
            params.(pe).(co).cellind = find(cellind);
            cellind = ismember(tag, tagdict('PositiveCurrentCollector'));
            params.(pe).(cc).cellind = find(cellind);
            params.(pe).(cc).extfaces = gen.positiveExtCurrentFaces;
            params.(pe).cellind = [params.(pe).(co).cellind; params.(pe).(cc).cellind];

            inputparams = setupElectrodes@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupCurrentCollectorBcCoupTerm(gen, inputparams, params)

            % grid of the current collector
            G = inputparams.G;
            extfaces = params.extfaces;

            facemap = G.mappings.facemap;
            invfacemap = zeros(G.parentGrid.topology.faces.num, 1);
            invfacemap(facemap) = (1 : G.topology.faces.num)';

            clear params
            params.bcfaces = invfacemap(extfaces);
            params.bccells = sum(G.topology.faces.neighbors(params.bcfaces, :), 2);

            inputparams = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupThermalModel(gen, inputparams, ~)
        % inputparams is instance of BatteryInputParams
        %
        % We recover the external coupling terms for the current collectors

            pG = gen.parentGrid;

            couplingfaces = gen.thermalExchangeFaces;
            couplingtags  = gen.thermalExchangeFacesTag;
            couplingcells = sum(pG.topology.faces.neighbors(couplingfaces, :), 2);

            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);

            thermal = 'ThermalModel';

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

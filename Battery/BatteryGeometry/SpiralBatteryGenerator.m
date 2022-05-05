classdef SpiralBatteryGenerator < BatteryGenerator
% Setup 3D grid with tab
    
    properties
            
        nwindings % number of windings in the spiral
        rInner    % Inner Radius (for the empty space in the middle)
        widthDict % dictionary of widths for each component. The required key names for the dictionary are
                  %                 - 'Separator'
                  %                 - 'NegativeActiveMaterial'
                  %                 - 'NegativeCurrentCollector'
                  %                 - 'PositiveActiveMaterial'
                  %                 - 'PositiveCurrentCollector'
        nrDict    % dicionary with number of cell in radial direction for each component (same keys as in widthDict).

        L            % length of the battery
        nas          % number of cells in the angular direction
        nL           % number of discretization cells in the longitudonal
        angleuniform % 

        tag     % cell-valued vector giving component number (indexing is given by tagdict)
        tagdict % dictionary giving the component number
    
        tabparams % parameters for the tab on the positive current collector
                  % if no tab, set as empty.
        
        positiveExtCurrentFaces
        negativeExtCurrentFaces
        thermalExchangeFaces
        thermalExchangeFacesTag

        celltbl
        widthLayer
        nWidthLayer
        heightLayer
        nHeightLayer

        tabwidths      % computed tab width (due to discretization, we cannot enforce the tab widths)
        windingnumbers % for the tabs (implmented only for aligned tabs now)
        
        use_thermal
    end
    
    methods
        
        function gen = SpiralBatteryGenerator()
            gen = gen@BatteryGenerator();  
        end
        
        function [paramobj, gen] = updateBatteryInputParams(gen, paramobj, params)
            
            gen.nwindings = params.nwindings;
            gen.rInner    = params.rInner;
            gen.widthDict = params.widthDict;
            gen.nrDict    = params.nrDict;
            gen.nas       = params.nas;
            gen.L         = params.L;
            gen.nL        = params.nL;
            gen.tabparams = params.tabparams;
            gen.angleuniform = params.angleuniform;
            
            gen.use_thermal = paramobj.use_thermal;
            
            [paramobj, gen] = gen.setupBatteryInputParams(paramobj, []);
            
        end
        
        function [paramobj, gen] = setupGrid(gen, paramobj, params)
    
            gen = spiralGrid(gen);
            paramobj.G = gen.G;
            
        end

        function UGrids = setupUnRolledGrids(gen, paramobj)
            
            G            = gen.G;
            widthLayer   = gen.widthLayer;
            nWidthLayer  = gen.nWidthLayer;
            heightLayer  = gen.heightLayer;
            nHeightLayer = gen.nHeightLayer;
            nas          = gen.nas;
            nL           = gen.nL;
            nwindings    = gen.nwindings;
            celltbl      = gen.celltbl;
            
            cartG = cartGrid([nas*nwindings, sum(nWidthLayer), nL]);

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

            vol = G.cells.volumes;

            celltbl = celltbl.removeInd({'cells', 'indi', 'indj'});

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

            cartG.nodes.coords = [xnode, ynode, znode];
            cartG = computeGeometry(cartG);
            
            %% we setup the mappings between the different grids
            
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am   = 'ActiveMaterial';
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
            
            compG = paramobj.(pe).(cc).G;
            comptag = tagdict('PositiveCurrentCollector');
            compCartG = genSubGrid(cartG, find(tag(cartInd) == comptag));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            Gs.CurrentCollector = compCartG;
            
            compG = paramobj.(pe).(am).G;
            comptag = tagdict('PositiveActiveMaterial');
            compCartG = genSubGrid(cartG, find(tag(cartInd) == comptag));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            Gs.ActiveMaterial = compCartG;
            
            UGrids.PositiveElectrode = Gs;
            clear Gs;
            
            compG = paramobj.(ne).(cc).G;
            comptag = tagdict('NegativeCurrentCollector');
            compCartG = genSubGrid(cartG, find(tag(cartInd) == comptag));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            Gs.CurrentCollector = compCartG;
            
            compG = paramobj.(ne).(am).G;
            comptag = tagdict('NegativeActiveMaterial');
            compCartG = genSubGrid(cartG, find(tag(cartInd) == comptag));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            Gs.ActiveMaterial = compCartG;
            
            UGrids.NegativeElectrode = Gs;
            clear Gs;
            
            compG = paramobj.(elyte).G;
            comptags = [tagdict('ElectrolyteSeparator'); tagdict('NegativeActiveMaterial'); tagdict('PositiveActiveMaterial')];
            compCartG = genSubGrid(cartG, find(ismember(tag(cartInd), comptags)));
            compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl);
            UGrids.Electrolyte = compCartG;
            
            UGrids.ThermalModel = UGrids.G;

        end
        
        function compCartG = setupCompCartGrid(gen, compG, compCartG, celltbl)

            compcelltbl.loccells = (1 : compG.cells.num)';
            compcelltbl.cells = compG.mappings.cellmap;
            compcelltbl = IndexArray(compcelltbl);
            compcelltbl = crossIndexArray(compcelltbl, celltbl, {'cells'});
            compcelltbl = compcelltbl.removeInd({'cells'});

            compcartcelltbl.cartloccells = (1 : compCartG.cells.num)';
            compcartcelltbl.cartcells = compCartG.mappings.cellmap;
            compcartcelltbl = IndexArray(compcartcelltbl);

            map = TensorMap();
            map.fromTbl = compcelltbl;
            map.toTbl = compcartcelltbl;
            map.mergefds = {'cartcells'};

            ind = map.getDispatchInd();
            
            compCartG.mappings.ind = ind;
            
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
        % 
        % We recover the external coupling terms for the current collectors

            G = gen.G;
            
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

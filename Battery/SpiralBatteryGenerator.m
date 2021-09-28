classdef SpiralBatteryGenerator < BatteryGenerator
% Setup 3D grid with tab
    
    properties
            
        
        nwindings % number of windings in the spiral
        r0        % "radius" at the middle
        widths    % vector of widths for each component indexed in following order
                  %                 - positive current collector
                  %                 - positive electrode
                  %                 - electrolyte separator 
                  %                 - negative electrode
                  %                 - negative current collector
                  %                 - separator
        L         % length of the battery
        nrs       % number of cell in radial direction for each component (same ordering as above).
        nas       % number of cells in the angular direction
        nL        % number of discretization cells in the longitudonal

        tag       % cell-valued vector giving component number (indexing is given by tagdict)
        tagdict   % dictionary giving the component number
        
        positiveExtCurrentFaces
        negativeExtCurrentFaces
        
        % Utility variables computed once and then shared by methods (should not be set)
        allparams;
        invcellmap;
        
        % Heat parameters
        externalHeatTransferCoefficientTab = 1e3;
        externalHeatTransferCoefficient = 1e3;
        
    end
    
    methods
        
        function gen = BatteryGenerator()
            gen = gen@BatteryGenerator();  
        end
        
        function paramobj = updateBatteryInputParams(gen, paramobj)
            paramobj = gen.setupBatteryInputParams(paramobj, []);
        end
        
        function [paramobj, gen] = setupGrid(gen, paramobj, params)
            
            gen.nwindings = params.nwindings;
            gen.r0        = params.r0     ;
            gen.widths    = params.widths ;
            gen.nrs       = params.nrs    ;
            gen.nas       = params.nas    ;
            gen.L         = params.L      ;
            gen.nL        = params.nL;
    
            gen = spiralGrid(gen);
            
            paramobj.G = gen.G;
            
        end

        function paramobj = setupElectrolyte(gen, paramobj, params)
            
            tagdict = gen.tagdict;
            tag = gen.tag;

            inds = [tagdict('PositiveActiveMaterial'); tagdict('ElectrolyteSeparator'); tagdict('NegativeActiveMaterial')];
            cellind = ismember(tag, ind);
            params.cellind = find(cellind);
            inds = tagdict('ElectrolyteSeparator'); 
            cellind = ismember(tag, ind);            
            parames.Separator.cellind = find(cellind);
            
            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);
            
        end
        
        function paramobj = setupElectrodes(gen, paramobj, params)

            
            % shortcuts 
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            eac = 'ElectrodeActiveComponent';

            tagdict = gen.tagdict;
            tag = gen.tag;
            
            cellind = ismember(tag, tagdict('NegativeActiveMaterial'));
            params.(ne).(eac).cellind = find(cellind);
            cellind = ismember(tag, tagdict('NegativeCurrentCollector'));
            params.(ne).(cc).cellind = find(cellind);
            params.(ne).(cc).extfaces = gen.negativeExtCurrentFaces;            
            params.(ne).cellind = [params.(ne).(eac).cellind; params.(ne).(cc).cellind];
            
            cellind = ismember(tag, tagdict('PositiveActiveMaterial'));
            params.(pe).(eac).cellind = find(cellind);
            cellind = ismember(tag, tagdict('PositiveCurrentCollector'));
            params.(pe).(cc).extfaces = gen.positiveExtCurrentFaces;
            params.(pe).cellind = [params.(pe).(eac).cellind; params.(pe).(cc).cellind];

            paramobj = setupElectrodes@BatteryGenerator(gen, paramobj, params);            
            
        end

        function paramobj = setupCurrentCollectorBcCoupTerm(gen, paramobj, params)
            
            G = paramobj.G;
            
            params.bcfaces = find(abs(yf - myf) < eps*1000);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

            paramobj = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, paramobj, params);
        
        end

        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
        % 
        % We recover the external coupling terms for the current collectors

            % shortcuts
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
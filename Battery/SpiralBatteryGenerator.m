classdef SpiralBatteryGenerator < BatteryGenerator
% Setup 3D grid with tab
    
    properties
            
        nwindings % number of windings in the spiral
        r0        % "radius" at the middle
        widthDict % dictionary of widths for each component. The required key names for the dictionary are
                  %                 - 'Separator'
                  %                 - 'NegativeActiveMaterial'
                  %                 - 'NegativeCurrentCollector'
                  %                 - 'PositiveActiveMaterial'
                  %                 - 'PositiveCurrentCollector'
        nrDict    % dicionary with number of cell in radial direction for each component (same keys as in widthDict).

        L         % length of the battery
        nas       % number of cells in the angular direction
        nL        % number of discretization cells in the longitudonal
        angleuniform % 
        tag       % cell-valued vector giving component number (indexing is given by tagdict)
        tagdict   % dictionary giving the component number
        
        positiveExtCurrentFaces
        negativeExtCurrentFaces
        thermalExchangeFaces
    end
    
    methods
        
        function gen = SpiralBatteryGenerator()
            gen = gen@BatteryGenerator();  
        end
        
        function paramobj = updateBatteryInputParams(gen, paramobj, params)
            
            gen.nwindings = params.nwindings;
            gen.r0        = params.r0;
            gen.widthDict = params.widthDict;
            gen.nrDict    = params.nrDict;
            gen.nas       = params.nas;
            gen.L         = params.L;
            gen.nL        = params.nL;
            gen.angleuniform = params.angleuniform;
            paramobj = gen.setupBatteryInputParams(paramobj, []);
            
        end
        
        function [paramobj, gen] = setupGrid(gen, paramobj, params)
    
            gen = spiralGrid(gen);
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
            params.(pe).(cc).cellind = find(cellind);
            params.(pe).(cc).extfaces = gen.positiveExtCurrentFaces;
            params.(pe).cellind = [params.(pe).(eac).cellind; params.(pe).(cc).cellind];

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
            couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);
            
            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);
            
            paramobj.ThermalModel.externalHeatTransferCoefficient = paramobj.ThermalModel.externalHeatTransferCoefficient*ones(numel(couplingfaces), 1);
            
        end
        
    end
    
end
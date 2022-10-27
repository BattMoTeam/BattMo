classdef BatteryGenerator2D < BatteryGenerator
% Setup 2D grid
    
    properties
        
        xlength = 1e-6*[10; 100; 50; 80; 10]; % length of components in x direction - default values
                                              % x(1) : length of negative current collector
                                              % x(2) : length of negative activie material
                                              % x(3) : length of separator
                                              % x(4) : length of positive current collector
                                              % x(5) : length of positive activie material
        
        ylength = 1e-2; % length in y direction - default values

        sepnx  = 30; % discretization number for negative current collector - default value
        nenx   = 30; % discretization number for negative activie material  - default value
        penx   = 30; % discretization number for separator                  - default value
        ccnenx = 10; % discretization number for positive current collector - default value
        ccpenx = 10; % discretization number for positive activie material  - default value

        ny = 10; % discretization number in y direction - default values

        use_thermal % flag, true if grid for thermal model should be setup.

        externalHeatTransferCoefficientTab = 1e3;  % heat transfer coefficient at tab boundary - default value 
        externalHeatTransferCoefficient = 1e3;     % heat transfer coefficient at boundary - default value 
        
        
    end
    
    methods
        
        function gen = BatteryGenerator2D()
            gen = gen@BatteryGenerator();  
        end
        
        function paramobj = updateBatteryInputParams(gen, paramobj)
            assert(paramobj.include_current_collectors, 'This geometry includes current collectors and input data appears to be missing for those');
            gen.use_thermal = paramobj.use_thermal;
            paramobj = gen.setupBatteryInputParams(paramobj, []);
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;
            ccnenx = gen.ccnenx;
            ccpenx = gen.ccpenx;
            
            ny = gen.ny;
            
            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            nx = sum(nxs);
            
            xlength = gen.xlength;
            ylength = gen.ylength;

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = ylength/ny;
            y = rldecode(y, ny);
            y = [0; cumsum(y)];

            G = tensorGrid(x, y);
            G = computeGeometry(G);
            
            paramobj.G = G;
            gen.G = G;
            
        end
        
        function paramobj = setupElectrolyte(gen, paramobj, params)
            
            sepnx  = gen.sepnx; 
            nenx   = gen.nenx; 
            penx   = gen.penx; 
            ccnenx = gen.ccnenx; 
            ccpenx = gen.ccpenx;

            ny     = gen.ny;
            
            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            nx = sum(nxs);
            
            istart = ccnenx + 1;
            ni = nenx + sepnx + penx;
            params.cellind = pickTensorCells(istart, ni, nx, ny);

            istart = ccnenx + nenx + 1;
            ni = sepnx;
            params.Separator.cellind = pickTensorCells(istart, ni, nx, ny);
            
            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);
            
        end
        

        function paramobj = setupElectrodes(gen, paramobj, params)
        % setup grid and coupling term
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            am = 'ActiveMaterial';
            
            sepnx  = gen.sepnx; 
            nenx   = gen.nenx; 
            penx   = gen.penx; 
            ccnenx = gen.ccnenx; 
            ccpenx = gen.ccpenx;     

            ny     = gen.ny;

            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            nx = sum(nxs);
            
            %% Negative electrode
            istart = 1;
            ni = ccnenx + nenx;
            params.(ne).cellind = pickTensorCells(istart, ni, nx, ny);
            
            % Negative electrode - current collector
            istart = 1;
            ni = ccnenx;
            params.(ne).(cc).cellind = pickTensorCells(istart, ni, nx, ny);
            
            % Negative electrode - electode active component
            istart = ccnenx + 1;
            ni = nenx;
            params.(ne).(am).cellind = pickTensorCells(istart, ni, nx, ny);
        
            %% Positive electrode
            istart = ccnenx + nenx + sepnx + 1;
            ni = penx + ccpenx;
            params.(pe).cellind = pickTensorCells(istart, ni, nx, ny);
            
            % Positive electrode - current collector
            istart = ccnenx + nenx + sepnx + penx + 1;
            ni = ccpenx;
            params.(pe).(cc).cellind = pickTensorCells(istart, ni, nx, ny);

            % Positive electrode - electode active component
            istart = ccnenx + nenx + sepnx + 1;
            ni = penx;
            params.(pe).(am).cellind = pickTensorCells(istart, ni, nx, ny);

            paramobj = setupElectrodes@BatteryGenerator(gen, paramobj, params);
            
        end

        function paramobj = setupCurrentCollectorBcCoupTerm(gen, paramobj, ~)
            
            G = paramobj.G;
            
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            params.bcfaces = find(abs(yf - myf) < eps*1000);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

            paramobj = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, paramobj, params);
        
        end

        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
        % 
        % We recover the external coupling terms for the current collectors
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            cc    = 'CurrentCollector';
            
            coupterm_necc = paramobj.(ne).(cc).externalCouplingTerm;
            facemap_necc  = paramobj.(ne).(cc).G.mappings.facemap;
            coupterm_pecc = paramobj.(pe).(cc).externalCouplingTerm;
            facemap_pecc  = paramobj.(pe).(cc).G.mappings.facemap;
            
            % the cooling is done on the external faces
            G = gen.G;
            extfaces = any(G.faces.neighbors == 0, 2);
            couplingfaces = find(extfaces);
            couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);
            
            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);
            
            % We assign the values different values to the tab
            tabtbl.faces = [facemap_necc(coupterm_necc.couplingfaces); 
                            facemap_pecc(coupterm_pecc.couplingfaces)];
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

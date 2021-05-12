classdef BatteryGenerator2D < BatteryGenerator
% Setup 2D grid
    
    properties
        
        sepnx  = 30;
        nenx   = 30;
        penx   = 30;
        ccnenx = 20;
        ccpenx = 20;

        ny = 10;
        
        J = 0.1;
        
        % externalHeatTransferCoefficientTab = 1e6;
        externalHeatTransferCoefficientTab = 1e3;
        externalHeatTransferCoefficient = 1e3;
        
    end
    
    methods
        
        function gen = BatteryGenerator2D()
            gen = gen@BatteryGenerator();  
        end
        
        function paramobj = updateBatteryInputParams(gen, paramobj)
            paramobj = gen.setupBatteryInputParams(paramobj, []);
            paramobj.J = gen.J;
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
            
            xlength = 1e-6*[10; 100; 50; 80; 10];
            ylength = 1e-2;

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
            params.sep.cellind = pickTensorCells(istart, ni, nx, ny);
            
            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);
            
        end
        

        function paramobj = setupElectrodes(gen, paramobj, params)
        % setup grid and coupling term
            
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
            params.ne.cellind = pickTensorCells(istart, ni, nx, ny);
            
            % Negative electrode - current collector
            istart = 1;
            ni = ccnenx;
            params.ne.cc.cellind = pickTensorCells(istart, ni, nx, ny);
            
            % Negative electrode - electode active component
            istart = ccnenx + 1;
            ni = nenx;
            params.ne.eac.cellind = pickTensorCells(istart, ni, nx, ny);
        
            %% Positive electrode
            istart = ccnenx + nenx + sepnx + 1;
            ni = penx + ccpenx;
            params.pe.cellind = pickTensorCells(istart, ni, nx, ny);
            
            % Positive electrode - current collector
            istart = ccnenx + nenx + sepnx + penx + 1;
            ni = ccpenx;
            params.pe.cc.cellind = pickTensorCells(istart, ni, nx, ny);

            % Positive electrode - electode active component
            istart = ccnenx + nenx + sepnx + 1;
            ni = penx;
            params.pe.eac.cellind = pickTensorCells(istart, ni, nx, ny);

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
            
            coupterm_necc = paramobj.ne.cc.couplingTerm;
            facemap_necc = paramobj.ne.cc.G.mappings.facemap;
            coupterm_pecc = paramobj.pe.cc.couplingTerm;
            facemap_pecc = paramobj.pe.cc.G.mappings.facemap;
            
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
            
            paramobj.thermal.externalHeatTransferCoefficient = coef;
            
        end
        
    end
    
    
end



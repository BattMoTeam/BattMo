classdef BatteryGenerator3D < BatteryGenerator
% Setup 3D grid with tab
    
    properties
            
        % Physical dimension
        
        xlength = 3e-2*[0.1; 1; 0.1];
        ylength = 1e-2*[0.1; 1; 0.1];
        zlength = 1e-6*[10; 100; 50; 80; 10];

        % Input current
        
        J = 0.1;

        % Shortcuts used below
        % ne    : NegativeElectrode
        % pe    : NegativeElectrode
        % eac   : ElectrodActiveComponent
        % cc    : CurrentCollector
        % elyte : Electrolyte
        
        % Discretization resolution in z-direction
        
        facz = 1;
        
        sep_nz    = 6;
        ne_eac_nz = 6;
        pe_eac_nz = 6;
        ne_cc_nz  = 4;
        pe_cc_nz  = 4;
        
        % Discretization resolution in x-direction
        
        facx = 1;
        
        int_elyte_nx  = 10; 
        ne_cc_nx = 5;
        pe_cc_nx = 5;

        % Discretization resolution in y-direction

        facy = 1;
        
        ne_cc_ny = 2;
        pe_cc_ny = 2;
        elyte_ny = 4;

        % Utility variables computed once and then shared by methods (should not be set)
        dimG
        elyte_nz;        
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
            
            
            gen = gen.applyResolutionFactors();
            
            nxs = [gen.ne_cc_nx; gen.int_elyte_nx; gen.pe_cc_nx];
            nys = [gen.ne_cc_ny; gen.elyte_ny; gen.pe_cc_ny];
            nzs = [gen.ne_cc_nz; gen.ne_eac_nz; gen.sep_nz; gen.ne_eac_nz; gen.pe_cc_nz];

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = ylength./nys;
            y = rldecode(y, nys);
            y = [0; cumsum(y)];

            z = zlength./nzs;
            z = rldecode(z, nzs);
            z = [0; cumsum(z)];

            G = tensorGrid(x, y, z);

            nx = sum(nxs);
            ny = sum(nys);
            nz = sum(nzs);

            gen.dimGrid = [nx; ny; nz];
            gen.G       = G;
            paramobj.G  = G;

        end
        
        function gen = applyResolutionFactors(gen)
            
            facz = gen.facz;
            
            gen.sep_nz    = facz*gen.sep_nz;
            gen.ne_eac_nz = facz*gen.ne_eac_n;
            gen.pe_eac_nz = facz*gen.pe_eac_n;
            gen.ne_cc_nz  = facz*gen.ne_cc_nz;
            gen.pe_cc_nz  = facz*gen.pe_cc_nz;
            
            facx = gen.facx;
            
            gen.int_elyte_nx = facx*gen.int_elyte_nx;
            gen.ne_cc_nx= facx*gen.ne_cc_n;
            gen.pe_cc_nx= facx*gen.pe_cc_n;

            facy = gen.facy;
            
            gen.ne_cc_ny= facy*gen.ne_cc_n;
            gen.pe_cc_ny= facy*gen.pe_cc_n;
            gen.elyte_ny= facy*gen.elyte_n;
            
        end

        function paramobj = setupElectrolyte(gen, paramobj, params)
            
            startSubGrid = [1; NegativeCurrentCollector_ny + 1; NegativeCurrentCollector_nz + 1];
            dimSubGrid   = [nx; elny; elnz];
            Electrolyte_Cells   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
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


    end
    
    
end



classdef BatteryGenerator3D < BatteryGenerator
% Setup 3D grid with tab
    
    properties
            
        % Physical dimension
        
        xlength = 3e-2*[0.1; 1; 0.1];
        ylength = 1e-2*[0.1; 1; 0.1];
        zlength = 1e-6*[10; 100; 50; 80; 10];

        % Input current
        
        J = 1e-4;

        % Shortcuts used below
        % ne    : NegativeElectrode
        % pe    : NegativeElectrode
        % eac   : ElectrodActiveComponent
        % cc    : CurrentCollector
        % elyte : Electrolyte
        
        % Discretization resolution in z-direction
        
        facz = 3;
        
        sep_nz    = 6;
        ne_eac_nz = 6;
        pe_eac_nz = 6;
        ne_cc_nz  = 4;
        pe_cc_nz  = 4;
        
        % Discretization resolution in x-direction
        
        facx = 3;
        
        int_elyte_nx  = 10; 
        ne_cc_nx = 5;
        pe_cc_nx = 5;

        % Discretization resolution in y-direction

        facy = 3;
        
        ne_cc_ny = 2;
        pe_cc_ny = 2;
        elyte_ny = 4;

        % Utility variables computed once and then shared by methods (should not be set)
        elyte_nz;
        allparams;
        invcellmap;
        
    end
    
    methods
        
        function gen = BatteryGenerator3D()
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
            
            gen.elyte_nz = gen.sep_nz + gen.ne_eac_nz + gen.pe_eac_nz;
            
            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.elyte_nz];
            allparams.elyte.cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);
            
            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + gen.ne_eac_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.sep_nz];
            allparams.elyte.sep.cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);
            
            %% setup gen.ne_eac

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.ne_eac_nz];
            allparams.ne.eac.cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            %% setup gen.pe_eac

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + gen.ne_eac_nz + gen.sep_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.pe_eac_nz];
            allparams.pe.eac.cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            %% setup gen.ne_cc

            startSubGrid = [1; gen.ne_cc_ny + 1; 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.ne_cc_nz];
            allparams.ne.cc.cellind1 = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            % We add the tab

            startSubGrid = [1; 1; 1];
            dimSubGrid   = [gen.ne_cc_nx; gen.ne_cc_ny; gen.ne_cc_nz];
            allparams.ne.cc.cellind2 = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            allparams.ne.cc.cellind = [allparams.ne.cc.cellind1; allparams.ne.cc.cellind2];

            %% setup gen.pe_cc

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + gen.elyte_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.pe_cc_nz];
            allparams.pe.cc.cellind1 = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            % We add the tab

            startSubGrid = [gen.ne_cc_nx + gen.int_elyte_nx + 1; gen.ne_cc_ny + gen.elyte_ny + 1; gen.ne_cc_nz + gen.elyte_nz + 1];
            dimSubGrid   = [gen.pe_cc_nx; gen.pe_cc_ny; gen.pe_cc_nz];
            allparams.pe.cc.cellind2 = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            allparams.pe.cc.cellind = [allparams.pe.cc.cellind1; allparams.pe.cc.cellind2];

            cellind = [allparams.elyte.cellind; allparams.ne.eac.cellind; allparams.pe.eac.cellind; allparams.ne.cc.cellind; allparams.pe.cc.cellind];

            rcellind = setdiff((1 : G.cells.num)', cellind);

            nGlob = G.cells.num;
            [G, cellmap, facemap, nodemap] = removeCells(G, rcellind);
            invcellmap = zeros(nGlob, 1);
            invcellmap(cellmap) = (1 : G.cells.num)';

            paramobj.G = G;
            
            gen.invcellmap = invcellmap;
            gen.allparams = allparams;
            gen.G = G;

        end
        
        function gen = applyResolutionFactors(gen)
            
            facz = gen.facz;
            
            gen.sep_nz    = facz*gen.sep_nz;
            gen.ne_eac_nz = facz*gen.ne_eac_nz;
            gen.pe_eac_nz = facz*gen.pe_eac_nz;
            gen.ne_cc_nz  = facz*gen.ne_cc_nz;
            gen.pe_cc_nz  = facz*gen.pe_cc_nz;
            
            facx = gen.facx;
            
            gen.int_elyte_nx = facx*gen.int_elyte_nx;
            gen.ne_cc_nx= facx*gen.ne_cc_nx;
            gen.pe_cc_nx= facx*gen.pe_cc_nx;

            facy = gen.facy;
            
            gen.ne_cc_ny= facy*gen.ne_cc_ny;
            gen.pe_cc_ny= facy*gen.pe_cc_ny;
            gen.elyte_ny= facy*gen.elyte_ny;
            
        end

        function paramobj = setupElectrolyte(gen, paramobj, params)
            
            params = gen.allparams.elyte;
            imap = gen.invcellmap;
            params.cellind = imap(params.cellind);
            params.sep.cellind = imap(params.sep.cellind);
            
            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);
            
        end
        
        function paramobj = setupElectrodes(gen, paramobj, params)
            
            params = gen.allparams;
            imap = gen.invcellmap;
            
            params.ne.eac.cellind = imap(params.ne.eac.cellind);
            params.ne.cc.cellind = imap(params.ne.cc.cellind);
            params.ne.cc.name = 'negative';
            params.ne.cellind = [params.ne.eac.cellind; params.ne.cc.cellind];
            
            params.pe.eac.cellind = imap(params.pe.eac.cellind);
            params.pe.cc.cellind = imap(params.pe.cc.cellind);
            params.pe.cc.name = 'positive';
            params.pe.cellind = [params.pe.eac.cellind; params.pe.cc.cellind];
            
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


    end
    
end
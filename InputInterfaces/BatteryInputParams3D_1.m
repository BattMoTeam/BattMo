classdef BatteryInputParams3D_1 < BatteryInputParams2D
    
    methods

        function params = setupSubModels(params)
            
            fac = 1/5;
            sepnx  = 30*fac;
            NegativeElectrodenx   = 30*fac;
            PositiveElectrodenx   = 30*fac;
            NegativeCurrentCollectornx = 20*fac;
            PositiveCurrentCollectornx = 20*fac;

            nElectrolyte = NegativeElectrodenx + sepnx + PositiveElectrodenx;
            
            nxs = [NegativeCurrentCollectornx; NegativeElectrodenx; sepnx; PositiveElectrodenx; PositiveCurrentCollectornx];
            ny = 5;
            nz = 1;
            
            xlength = 1e-6*[10; 100; 50; 80; 10];
            ylength = 1e-2;
            zlength = 1;
            
            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = ylength/ny;
            y = rldecode(y, ny);
            y = [0; cumsum(y)];
            
            z = zlength/nz;
            z = rldecode(z, nz);
            z = [0; cumsum(z)];

            G = tensorGrid(x, y, z);
            G = computeGeometry(G);
            
            nx = sum(nxs);

            dimGlobGrid = [nx; ny; nz];

            %% setup Electrolyte

            startSubGrid = [NegativeCurrentCollectornx + 1; 1; 1];
            dimSubGrid   = [nElectrolyte; ny; nz];
            Electrolytecells   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup NegativeElectrode

            startSubGrid = [NegativeCurrentCollectornx + 1; 1; 1];
            dimSubGrid   = [NegativeElectrodenx; ny; nz];
            NegativeElectrodecells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup PositiveElectrode

            startSubGrid = [NegativeCurrentCollectornx + NegativeElectrodenx + sepnx + 1; 1; 1];
            dimSubGrid   = [PositiveElectrodenx; ny; nz];
            PositiveElectrodecells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
            %% setup sep

            startSubGrid = [NegativeCurrentCollectornx + NegativeElectrodenx + 1; 1;  1];
            dimSubGrid   = [sepnx; ny; nz];
            sepcells     = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup NegativeCurrentCollector

            startSubGrid = [1; 1; 1];
            dimSubGrid   = [NegativeCurrentCollectornx; ny; nz];
            NegativeCurrentCollectorcells    = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup PositiveCurrentCollector

            startSubGrid = [NegativeCurrentCollectornx + nElectrolyte + 1; 1; 1];
            dimSubGrid   = [PositiveCurrentCollectornx; ny; nz];
            PositiveCurrentCollectorcells    = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% assign the parameters
            
            params.G = G;

            params.Electrolyte = orgLiPF6('Electrolyte', G, Electrolytecells);
            params.NegativeElectrode    = graphiteElectrode('NegativeElectrode', G, NegativeElectrodecells);
            params.PositiveElectrode    = nmc111Electrode('PositiveElectrode', G, PositiveElectrodecells);
            params.NegativeCurrentCollector  = currentCollector('NegativeCurrentCollector', G, NegativeCurrentCollectorcells);
            params.PositiveCurrentCollector  = currentCollector('PositiveCurrentCollector', G, PositiveCurrentCollectorcells);
            params.sep   = celgard2500('sep', G, sepcells);
            
        end
        
        function coupTerm = setupNegativeCurrentCollectorBcCoupTerm(params)

            NegativeCurrentCollector = params.NegativeCurrentCollector;
            G = NegativeCurrentCollector.G;

            % We pick up the faces at the top of NegativeCurrentCollector
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf-myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'NegativeCurrentCollector'};
            coupTerm = couplingTerm('bc-NegativeCurrentCollector', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupPositiveCurrentCollectorBcCoupTerm(params)

            PositiveCurrentCollector = params.PositiveCurrentCollector;
            G = PositiveCurrentCollector.G;

            % We pick up the faces at the bottom of PositiveCurrentCollector
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces =  find(abs(yf-myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'PositiveCurrentCollector'};
            coupTerm = couplingTerm('bc-PositiveCurrentCollector', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end


    end
    
end

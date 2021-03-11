classdef BatteryInputParams3D_1 < BatteryInputParams2D
    
    methods

        function params = setupSubModels(params)
            
            fac = 1/5;
            sepnx  = 30*fac;
            nenx   = 30*fac;
            penx   = 30*fac;
            ccnenx = 20*fac;
            ccpenx = 20*fac;

            nelyte = nenx + sepnx + penx;
            
            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
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

            %% setup elyte

            startSubGrid = [ccnenx + 1; 1; 1];
            dimSubGrid   = [nelyte; ny; nz];
            elytecells   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup ne

            startSubGrid = [ccnenx + 1; 1; 1];
            dimSubGrid   = [nenx; ny; nz];
            necells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup pe

            startSubGrid = [ccnenx + nenx + sepnx + 1; 1; 1];
            dimSubGrid   = [penx; ny; nz];
            pecells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
            %% setup sep

            startSubGrid = [ccnenx + nenx + 1; 1;  1];
            dimSubGrid   = [sepnx; ny; nz];
            sepcells     = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup ccne

            startSubGrid = [1; 1; 1];
            dimSubGrid   = [ccnenx; ny; nz];
            ccnecells    = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup ccpe

            startSubGrid = [ccnenx + nelyte + 1; 1; 1];
            dimSubGrid   = [ccpenx; ny; nz];
            ccpecells    = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% assign the parameters
            
            params.G = G;

            params.elyte = orgLiPF6('elyte', G, elytecells);
            params.ne    = graphiteElectrode('ne', G, necells);
            params.pe    = nmc111Electrode('pe', G, pecells);
            params.ccne  = currentCollector('ccne', G, ccnecells);
            params.ccpe  = currentCollector('ccpe', G, ccpecells);
            params.sep   = celgard2500('sep', G, sepcells);
            
        end
        
        function coupTerm = setupCcneBcCoupTerm(params)

            ccne = params.ccne;
            G = ccne.G;

            % We pick up the faces at the top of ccne
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf-myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccne'};
            coupTerm = couplingTerm('bc-ccne', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupCcpeBcCoupTerm(params)

            ccpe = params.ccpe;
            G = ccpe.G;

            % We pick up the faces at the bottom of ccpe
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces =  find(abs(yf-myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccpe'};
            coupTerm = couplingTerm('bc-ccpe', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end


    end
    
end

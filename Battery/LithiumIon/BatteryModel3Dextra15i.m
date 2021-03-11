classdef BatteryModel3Dextra15i < BatteryModel15i

    methods
        
        function model = BatteryModel3Dextra15i(varargin)
            model = model@BatteryModel15i(varargin{:});
        end
        function validforces = getValidDrivingForces(model)
            validforces=struct('src', [], 'stopFunction', []); 
        end
        
        function model = setupBatteryComponents(model)
            
            fac = 1/5;
            sepnx  = 30*fac;
            nenx   = 30*fac;
            penx   = 30*fac;
            ccnenx = 20*fac;
            ccpenx = 20*fac;

            nelyte = nenx + sepnx + penx;
            
            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            ny = 4;
            nz = 4;
            
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

            %% 

            model.G = G;
            
            submodels = {};
            submodels{end + 1} = orgLiPF6('elyte'       , G, elytecells);
            submodels{end + 1} = graphiteElectrode('ne' , G, necells);
            submodels{end + 1} = nmc111Electrode('pe'   , G, pecells);
            submodels{end + 1} = currentCollector('ccne', G, ccnecells);
            submodels{end + 1} = currentCollector('ccpe', G, ccpecells);
            submodels{end + 1} = celgard2500('sep'      , G, sepcells);
            model.SubModels = submodels;
            
        end
        
        function coupTerm = setupCcneBcCoupTerm(model)

            ccne = model.ccne;
            G = ccne.G;

            % We pick up the faces at the top of ccne
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(yf >= (1 - eps)*myf);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccne'};
            coupTerm = couplingTerm('bc-ccne', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupCcpeBcCoupTerm(model)

            ccpe = model.ccpe;
            G = ccpe.G;

            % We pick up the faces at the bottom of ccpe
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(yf >= (1 - eps)*myf);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccpe'};
            coupTerm = couplingTerm('bc-ccpe', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        
    end
    
end

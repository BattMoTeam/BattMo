classdef BatteryModel3D < BatteryModel

    methods
        
        function model = BatteryModel3D(varargin)
            model = model@BatteryModel(varargin{:});
        end
        
        function model = setupBatteryComponents(model)
            
            sepnz  = 10;
            nenz   = 10;
            penz   = 10;

            elnz = sepnz + nenz + penz;

            ccnenx = 2;
            ccpenx = 2;
            elnx   = 10;

            ccneny = 2;
            ccpeny = 2;
            elny   = 10;

            nxs = [ccnenx; elnx; ccpenx];
            nys = [ccneny; elny; ccpeny];
            nzs = [nenz; sepnz; penz];

            xlength = 3e-2*[0.1; 1; 0.1];
            ylength = 1e-2*[0.1; 1; 0.1];
            zlength = 1e-6*[100; 50; 80];

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

            dimGlobGrid = [nx; ny; nz];

            %% setup elyte

            startSubGrid = [ccnenx + 1; ccpeny + 1; 1];
            dimSubGrid   = [elnx; elny; elnz];
            elytecells   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup ne

            startSubGrid = [ccnenx + 1; ccpeny + 1; 1];
            dimSubGrid   = [elnx; elny; nenz];
            necells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup pe

            startSubGrid = [ccnenx + 1; ccpeny + 1; nenz + sepnz + 1];
            dimSubGrid   = [elnx; elny; penz];
            pecells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
            %% setup sep

            startSubGrid = [ccnenx + 1; ccpeny + 1; nenz +  1];
            dimSubGrid   = [elnx; elny; sepnz];
            sepcells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup ccne

            startSubGrid = [1; ccpeny + 1; 1];
            dimSubGrid   = [ccnenx; elny + ccneny; nenz];
            ccnecells    = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup ccpe

            startSubGrid = [ccnenx + elnx + 1; 1; nenz + sepnz];
            dimSubGrid   = [ccpenx; elny + ccpeny; penz];
            ccpecells    = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% 
            
            cells = [elytecells; necells; pecells; ccnecells; ccpecells];

            rcells = setdiff((1 : G.cells.num)', cells);

            nGlob = G.cells.num;

            [G, cellmap, facemap, nodemap] = removeCells(G, rcells);

            model.G = G;
            
            invcellmap = zeros(nGlob, 1);
            invcellmap(cellmap) = (1 : G.cells.num)';
            
            submodels = {};
            submodels{end + 1} = orgLiPF6('elyte'       , G, invcellmap(elytecells));
            submodels{end + 1} = graphiteElectrode('ne' , G, invcellmap(necells));
            submodels{end + 1} = nmc111Electrode('pe'   , G, invcellmap(pecells));
            submodels{end + 1} = currentCollector('ccne', G, invcellmap(ccnecells));
            submodels{end + 1} = currentCollector('ccpe', G, invcellmap(ccpecells));
            submodels{end + 1} = celgard2500('sep'      , G, invcellmap(sepcells));
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
            myf = min(yf);
            faces = find(yf <= (1 + eps)*myf);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccpe'};
            coupTerm = couplingTerm('bc-ccpe', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        
    end
    
end

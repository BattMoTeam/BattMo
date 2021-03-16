classdef BatteryInputParams3D_2 < BatteryInputParams2D
    
    methods

        function params = setupSubModels(params)
            
            % physical dimensions in the three cartesian directions 
            xlength = 3e-2*[0.1; 1; 0.1];
            ylength = 1e-2*[0.1; 1; 0.1];
            zlength = 1e-6*[10; 100; 50; 80; 10];
            
            % discretization in z direction 
            facz= 1;
            sepnz  = 6*facz;
            nenz   = 6*facz;
            penz   = 6*facz;
            ccnenz = 4*facz;
            ccpenz = 4*facz;
            
            elnz = sepnz + nenz + penz;

            % discretization in x direction 
            facx = 1;
            intelnx = 10*facx; % "interior" of electrolyte
            ccnenx  = 5*facx;
            ccpenx  = 5*facx;
            
            % discretization in y direction
            facy = 1;
            ccneny = 2*facy;
            ccpeny = 2*facy;
            elny   = 4*facy;

            nxs = [ccnenx; intelnx; ccpenx];
            nys = [ccneny; elny; ccpeny];
            nzs = [ccnenz; nenz; sepnz; penz; ccpenz];

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

            startSubGrid = [1; ccneny + 1; ccnenz + 1];
            dimSubGrid   = [nx; elny; elnz];
            elytecells   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup ne

            startSubGrid = [1; ccneny + 1; ccnenz + 1];
            dimSubGrid   = [nx; elny; nenz];
            necells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup pe

            startSubGrid = [1; ccneny + 1; ccnenz + nenz + sepnz + 1];
            dimSubGrid   = [nx; elny; penz];
            pecells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
            %% setup sep

            startSubGrid = [1; ccneny + 1; ccnenz + nenz + 1];
            dimSubGrid   = [nx; elny; sepnz];
            sepcells     = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup ccne

            startSubGrid = [1; ccneny + 1; 1];
            dimSubGrid   = [nx; elny; ccnenz];
            ccnecells1   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            % we add the tab
            startSubGrid = [1; 1; 1];
            dimSubGrid   = [ccnenx; ccneny; ccnenz];
            ccnecells2   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
            ccnecells = [ccnecells1; ccnecells2];

            %% setup ccpe

            startSubGrid = [1; ccneny + 1; ccnenz + elnz + 1];
            dimSubGrid   = [nx; elny; ccpenz];
            ccpecells1   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            % we add the tab
            startSubGrid = [ccnenx + intelnx + 1; ccneny + elny + 1; ccnenz + elnz + 1];
            dimSubGrid   = [ccpenx; ccpeny; ccpenz];
            ccpecells2   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
            ccpecells = [ccpecells1; ccpecells2];
            
            %% 
            
            cells = [elytecells; necells; pecells; ccnecells; ccpecells];

            rcells = setdiff((1 : G.cells.num)', cells);

            nGlob = G.cells.num;
            [G, cellmap, facemap, nodemap] = removeCells(G, rcells);
            invcellmap = zeros(nGlob, 1);
            invcellmap(cellmap) = (1 : G.cells.num)';
            
            params.G = G;
            
            params.elyte = orgLiPF6('elyte'       , G, invcellmap(elytecells));
            params.ne    = graphiteElectrode('ne' , G, invcellmap(necells));
            params.pe    = nmc111Electrode('pe'   , G, invcellmap(pecells));
            params.ccne  = currentCollector('ccne', G, invcellmap(ccnecells));
            params.ccpe  = currentCollector('ccpe', G, invcellmap(ccpecells));
            params.sep   = celgard2500('sep'      , G, invcellmap(sepcells));
            
        end
        
        function coupTerm = setupCcneBcCoupTerm(model)

            ccne = model.ccne;
            G = ccne.G;

            % We pick up the faces at the top of ccne
            yf = G.faces.centroids(:, 2);
            myf = min(yf);
            faces = find(abs(yf - myf) < eps*1000);
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
            faces = find(abs(yf - myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccpe'};
            coupTerm = couplingTerm('bc-ccpe', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

    end
    
end

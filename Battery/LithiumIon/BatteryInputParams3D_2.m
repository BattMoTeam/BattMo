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
            NegativeElectrode_nz   = 6*facz;
            PositiveElectrode_nz   = 6*facz;
            NegativeCurrentCollector_nz = 4*facz;
            PositiveCurrentCollector_nz = 4*facz;
            
            elnz = sepnz + NegativeElectrode_nz + PositiveElectrode_nz;

            % discretization in x direction 
            facx = 1;
            intelnx = 10*facx; % "interior" of electrolyte
            NegativeCurrentCollector_nx  = 5*facx;
            PositiveCurrentCollector_nx  = 5*facx;
            
            % discretization in y direction
            facy = 1;
            NegativeCurrentCollector_ny = 2*facy;
            PositiveCurrentCollector_ny = 2*facy;
            elny   = 4*facy;

            nxs = [NegativeCurrentCollector_nx; intelnx; PositiveCurrentCollector_nx];
            nys = [NegativeCurrentCollector_ny; elny; PositiveCurrentCollector_ny];
            nzs = [NegativeCurrentCollector_nz; NegativeElectrode_nz; sepnz; PositiveElectrode_nz; PositiveCurrentCollector_nz];

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

            %% setup Electrolyte

            startSubGrid = [1; NegativeCurrentCollector_ny + 1; NegativeCurrentCollector_nz + 1];
            dimSubGrid   = [nx; elny; elnz];
            Electrolyte_Cells   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup NegativeElectrode

            startSubGrid = [1; NegativeCurrentCollector_ny + 1; NegativeCurrentCollector_nz + 1];
            dimSubGrid   = [nx; elny; NegativeElectrode_nz];
            NegativeElectrode_Cells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup PositiveElectrode

            startSubGrid = [1; NegativeCurrentCollector_ny + 1; NegativeCurrentCollector_nz + NegativeElectrode_nz + sepnz + 1];
            dimSubGrid   = [nx; elny; PositiveElectrode_nz];
            PositiveElectrode_Cells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
            %% setup sep

            startSubGrid = [1; NegativeCurrentCollector_ny + 1; NegativeCurrentCollector_nz + NegativeElectrode_nz + 1];
            dimSubGrid   = [nx; elny; sepnz];
            sepcells     = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            %% setup NegativeCurrentCollector

            startSubGrid = [1; NegativeCurrentCollector_ny + 1; 1];
            dimSubGrid   = [nx; elny; NegativeCurrentCollector_nz];
            NegativeCurrentCollector_Cells1   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            % we add the tab
            startSubGrid = [1; 1; 1];
            dimSubGrid   = [NegativeCurrentCollector_nx; NegativeCurrentCollector_ny; NegativeCurrentCollector_nz];
            NegativeCurrentCollector_Cells2   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
            NegativeCurrentCollector_Cells = [NegativeCurrentCollector_Cells1; NegativeCurrentCollector_Cells2];

            %% setup PositiveCurrentCollector

            startSubGrid = [1; NegativeCurrentCollector_ny + 1; NegativeCurrentCollector_nz + elnz + 1];
            dimSubGrid   = [nx; elny; PositiveCurrentCollector_nz];
            PositiveCurrentCollector_Cells1   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

            % we add the tab
            startSubGrid = [NegativeCurrentCollector_nx + intelnx + 1; NegativeCurrentCollector_ny + elny + 1; NegativeCurrentCollector_nz + elnz + 1];
            dimSubGrid   = [PositiveCurrentCollector_nx; PositiveCurrentCollector_ny; PositiveCurrentCollector_nz];
            PositiveCurrentCollector_Cells2   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);
            
            PositiveCurrentCollector_Cells = [PositiveCurrentCollector_Cells1; PositiveCurrentCollector_Cells2];
            
            %% 
            
            cells = [Electrolyte_Cells; NegativeElectrode_Cells; PositiveElectrode_Cells; NegativeCurrentCollector_Cells; PositiveCurrentCollector_Cells];

            rcells = setdiff((1 : G.cells.num)', cells);

            nGlob = G.cells.num;
            [G, cellmap, facemap, nodemap] = removeCells(G, rcells);
            invcellmap = zeros(nGlob, 1);
            invcellmap(cellmap) = (1 : G.cells.num)';
            
            params.G = G;
            
            params.Electrolyte = orgLiPF6('Electrolyte'       , G, invcellmap(Electrolyte_Cells));
            params.NegativeElectrode    = GraphiteElectrode('NegativeElectrode' , G, invcellmap(NegativeElectrode_Cells));
            params.PositiveElectrode    = nmc111Electrode('PositiveElectrode'   , G, invcellmap(PositiveElectrode_Cells));
            params.NegativeCurrentCollector  = CurrentCollector('NegativeCurrentCollector', G, invcellmap(NegativeCurrentCollector_Cells));
            params.PositiveCurrentCollector  = CurrentCollector('PositiveCurrentCollector', G, invcellmap(PositiveCurrentCollector_Cells));
            params.sep   = celgard2500('sep'      , G, invcellmap(sepcells));
            
        end
        
        function coupTerm = setupNegativeCurrentCollectorBcCoupTerm(model)

            NegativeCurrentCollector = model.NegativeCurrentCollector;
            G = NegativeCurrentCollector.G;

            % We pick up the faces at the top of NegativeCurrentCollector
            yf = G.faces.centroids(:, 2);
            myf = min(yf);
            faces = find(abs(yf - myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'NegativeCurrentCollector'};
            coupTerm = couplingTerm('bc-NegativeCurrentCollector', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupPositiveCurrentCollectorBcCoupTerm(model)

            PositiveCurrentCollector = model.PositiveCurrentCollector;
            G = PositiveCurrentCollector.G;

            % We pick up the faces at the bottom of PositiveCurrentCollector
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf - myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'PositiveCurrentCollector'};
            coupTerm = couplingTerm('bc-PositiveCurrentCollector', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

    end
    
end

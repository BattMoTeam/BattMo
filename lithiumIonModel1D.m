classdef lithiumIonModel1D < lithiumIonModel

    methods
        function obj = lithiumIonModel1D(varargin)
            obj = obj@lithiumIonModel(varargin{:});
        end

        function obj = setupComponents(obj, SOC)
            
            T = obj.T;
            
            %% Define battery components
            
            sepnx  = 10;
            nenx   = 10;
            penx   = 10;
            ccnenx = 5;
            ccpenx = 5;
            
            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            ny = 1;

            xlength = 1e-6*ones(5, 1);
            ylength = 1;
            
            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];
            
            y = ylength/ny;
            y = rldecode(y, ny);
            y = [0; cumsum(y)];
            
            G = tensorGrid(x, y);
            G = computeGeometry(G);
            obj.G = G;
            
            obj.componentnames = {'elyte', 'ne', 'pe', 'ccne', 'ccpe'};
            
            %% setup elyte
            nx = sum(nxs); 

            istart = ccnenx + 1;
            ni = nenx + sepnx + penx;
            cells = pickTensorCells(istart, ni, nx, ny);
            obj.elyte = orgLiPF6('elyte', 1000, T, G, cells);
            
            %% setup ne
            istart = ccnenx + 1;
            cells = pickTensorCells(istart, nenx, nx, ny);
            obj.ne = graphiteElectrode('ne', SOC, T, G, cells);
            
            %% setup pe
            istart = ccnenx + nenx + sepnx + 1;
            cells = pickTensorCells(istart, penx, nx, ny);
            obj.pe = nmc111Electrode('pe', SOC, T, G, cells);

            %% setup ccne
            istart = 1;
            cells = pickTensorCells(istart, ccnenx, nx, ny);
            obj.ccne = currentCollector('ccne', T, G, cells);
            
            %% setup ccpe
            istart = ccnenx + nenx + sepnx + penx + 1;
            cells = pickTensorCells(istart, ccpenx, nx, ny);
            obj.ccpe = currentCollector('ccpe', T, G, cells);

            %% setup sep 
            istart = ccnenx + nenx + 1;
            cells = pickTensorCells(istart, sepnx, nx, ny);
            obj.sep = celgard2500('sep', G, cells);
        end
        
        
        function coupTerm = setupCcneBcCoupTerm(obj)
            
            G = obj.ccne.Grid;
            % We pick up the faces at the top of Cccne
            xf = G.faces.centroids(:, 1);
            mxf = min(xf);
            if mxf ~= 0
                faces = find(xf < (1 + eps)*mxf);
            else
                faces = find(xf < eps);
            end
            
            cells = sum(G.faces.neighbors(faces, :), 2);
            
            compnames = {'ccne'};
            coupTerm = couplingTerm('bc-ccne', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;
            
        end
        
        function coupTerm = setupCcpeBcCoupTerm(obj)
            
            G = obj.ccpe.Grid;
            % We pick up the faces at the top of Cccpe
            xf = G.faces.centroids(:, 1);
            mxf = max(xf);
            if mxf ~= 0
                faces = find(xf > (1 - eps)*mxf);            
            else
                faces = find(xf > - eps );
            end            
            cells = sum(G.faces.neighbors(faces, :), 2);
            
            compnames = {'ccpe'};
            coupTerm = couplingTerm('bc-ccpe', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;
            
        end
        
    end
    
end


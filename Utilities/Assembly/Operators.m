classdef Operators

    properties

        G

        hT  % Half-transmissibilities
        sgn % sign for the external faces (only)

        internalConn % flag for internal faces
        
        M % cell to cellface matrix mapping
        P % interior-cell to cellface matrix mapping

        cellfluxP % Reconstruction operator from (integrated) face fluxes to cell vector fluxes (grad)
        cellfluxS % Square root operator used to compute norm of vector flux on cells

        Grad   % Grad operator (see setupOperatorsTPFA)
        AccDiv % AccDiv operator (see setupOperatorsTPFA)
        
    end
    
    methods
        
        function op = Operators(G, varargin)
            
            opts = struct('assembleCellFluxOperator', false);
            opts = merge_options(opts, varargin{:});

            op.G = G;
            
            nc = G.cells.num;
            cst = ones(nc, 1);
            rock = struct('perm', cst, 'poro', cst);
            output = setupOperatorsTPFA(G, rock);

            op.internalConn = output.internalConn;
            op.Grad         = output.Grad;
            op.AccDiv       = output.AccDiv;
            
            op.hT = computeTrans(G, rock);

            tbls = setupSimpleTables(G);
            cellfacetbl = tbls.cellfacetbl;
            celltbl     = tbls.celltbl;
            facetbl     = tbls.facetbl;
            intfacetbl  = tbls.intfacetbl;

            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = cellfacetbl;
            map.mergefds = {'cells'};
            map = map.setup();

            M = SparseTensor();
            M = M.setFromTensorMap(map);
            
            op.M = M.getMatrix;

            map = TensorMap();
            map.fromTbl = cellfacetbl;
            map.toTbl = intfacetbl;
            map.mergefds = {'faces'};
            map = map.setup();

            P = SparseTensor();
            P = P.setFromTensorMap(map);
            
            op.P = P.getMatrix;
            
            %% setup the sign for *external* faces
            cells  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
            faces  = G.cells.faces(:, 1);
            
            extfaces = find(any(G.faces.neighbors == 0, 2));
            extcells = sum(G.faces.neighbors(extfaces, :), 2);

            extsgn = 2*(extcells == G.faces.neighbors(extfaces, 1)) - 1;

            sgn = nan(G.faces.num, 1);
            sgn(extfaces) = extsgn;

            op.sgn = sgn;
            
            %% setup cell flux reconstruction operator
            if opts.assembleCellFluxOperator
                
                output = getCellFluxOperatorsAll(G);
                op.cellfluxP = output.P;
                op.cellfluxS = output.S;
                
            end

        end

        function jsq = cellFlux(op, j)
        % Compute cell-valued square of the norm of a face-valued flux (j)
            
            P = op.cellfluxP
            S = op.cellfluxS

            j   = P*j;
            jsq = j.^2;
            jsq = S*jsq;
            
        end
        
        function harmfacevalue = harmFace(op, cellvalue)

            P  = op.P;
            hT = op.hT;
            M  = op.M;
            
            harmfacevalue = 1./(P*(1./(hT.*(M*cellvalue))));

        end

        function vols = getCellVolumes(op)

            vols = op.G.cells.volumes;
            
        end

        function areas = getFaceAreas(op)

            areas = op.G.faces.areas;
            
        end
        
        function [T, cells] = harmFaceBC(op, cvalue, faces)

            G = op.G;
            
            cells = sum(G.faces.neighbors(faces, :), 2);
            cn    = sqrt(sum((G.faces.centroids(faces, :) - G.cells.centroids(cells, :)).^2, 2));
            t = G.faces.areas(faces)./cn;
            T = t.*cvalue(cells);
            
        end

    end

end

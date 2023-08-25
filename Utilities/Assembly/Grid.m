classdef Grid
%
% Grid structure that can be used for a parent grid in a sub-grid (see SubGrid class).
%
% The geometrical properties are all stored in a immutable handle (see TwoPointFiniteVolumeGeometry class).
%
% By using handle, we can share the grid as a parent grid between different sub-grid instances.
%
% By using immutable properties, we avoid any collateral effects arising from the use of handler.
%
% The assembly is done using the function listed here that are at the end. Note that, typically, we will used the
% sub-grid versions of those as this class is mainly meant to be used as parent grid
%
% - function nc                     = getNumberOfCells(grid)          % returns number of cells
% - function vols                   = getVolumes(grid)                % returns cells volumes
% - function areas                  = getFaceAreas(grid)              % return face areas
% - function u                      = getHarmFace(grid, c)            % return harmonic average on the face from cell values
% - function [bchT, bccells, bcsgn] = getBcHarmFlux(grid, u, bcfaces) % return weighted half-transmissibilities form 
% - function jsq                    = getCellFluxNorm(grid, u)        % return norms of flux (no AD compliant version for the moment!)


    properties (SetAccess = protected)

        % Description of topology using MRST grid structures: The fieds are
        % - topology.cells.facePos
        % - topology.cells.faces
        % - topology.cells.num
        % - topology.faces.nodePos
        % - topology.faces.nodes
        % - topology.faces.num
        % - topology.faces.neighbors
        % - topology.nodes.num
        % - topology.griddim
        topology

        % Instance of TwoPointFiniteVolumeGeometry (handle with immutable properties)
        tPFVgeometry

        % Helper structures that are used in the computation of the gradient, divergence, harmonic face averages, ...
        % - helpers.diffop.grad                 (sparse matrix used in getGradient)
        % - helpers.diffop.div                  (sparse matrix used in getDiv)
        % - helpers.trans.D                     (sparse matrix used in getHarmFace method)
        % - helpers.trans.P                     (sparse matrix used in getHarmFace method)
        % - helpers.trans.S                     (sparse matrix used in getHarmFace method)
        % - helpers.extfaces.faces              (index of the external faces, sub-grid indexing)
        % - helpers.extfaces.cells              (index of the corresponding cells, sub-grid indexing)
        % - helpers.extfaces.halfTransParentInd (index of the corresponding half-transmissibility values in parent grid indexing)
        % - helpers.faceextfacemap              (mapping from face to extface, sub-grid indexing)
        helpers

        % Operators to compute norm of the flux velocity at the cell centers, see getCellFluxNorm
        % - cellFluxOperators.P
        % - cellFluxOperators.S
        cellFluxOperators

        % Helpers for AD assembly
        matrixOperators
        vectorHelpers
        
    end

    methods

        function grid = Grid(G, varargin)
        % We initialize the grid using a MRST grid structure

            opt = struct('faceArea', []);
            opt = merge_options(opt, varargin{:});
            
            topology.cells.facePos   = G.cells.facePos;
            topology.cells.faces     = G.cells.faces;
            topology.cells.num       = G.cells.num;
            topology.faces.nodePos   = G.faces.nodePos;
            topology.faces.nodes     = G.faces.nodes;
            topology.faces.num       = G.faces.num;
            topology.faces.neighbors = G.faces.neighbors;
            topology.nodes.num       = G.nodes.num;
            topology.griddim         = G.griddim;
            
            grid.topology = topology;

            grid = grid.setupHelpers();

            % Setup  tPFVgeometry
            nodecoords = reshape(G.nodes.coords', [], 1);
            tPFVgeometry = grid.computeTPFVgeometry(nodecoords, ...
                                                    'faceArea', opt.faceArea);
            grid = grid.assignTPFVgeometry(tPFVgeometry);
            
        end

        function grid = assignTPFVgeometry(grid, tPFVgeometry)
        % Assign given TwoPointFiniteVolumeGeometry object given by tPFVgeometry to grid
            
            grid.tPFVgeometry = TwoPointFiniteVolumeGeometry(tPFVgeometry);
            
        end
        
        function [tPFVgeometry, grid] = computeTPFVgeometry(grid, nodecoords, varargin)
        % Given node coordinates (and faceArea - only required for 1D problem), we return the tPFVgeometry as a
        % TwoPointFiniteVolumeGeometry instance

            opt = struct('faceArea', [], ...
                         'assemblyType', 'MRST');
            opt = merge_options(opt, varargin{:});
            
            switch opt.assemblyType

              case 'MRST'
                
                tPFVgeometry = computeTPFVgeometryMRST(grid, nodecoords, opt.faceArea);
                
              case 'AD'                

                if isempty(grid.matrixOperators)

                    grid = grid.setupADhelpers();
                    
                end
                
                tPFVgeometry = computeTPFVgeometryAD(grid, nodecoords, opt.faceArea);

              otherwise
                
                error('assembly type not recognized');
                
            end
            
        end
        
        function tPFVgeometry = computeTPFVgeometryMRST(grid, nodecoords, faceArea)
        % Compute the geometrical properties using genuine MRST functions. The input can only be vector of double (no AD support)
        % The argument faceArea is only needed for 1D case
            
            mrstG = grid.topology;
            d = mrstG.griddim;
            
            if  mrstG.griddim == 1 & ~isempty(faceArea)
                farea = faceArea;
            else
                farea = 1;
            end
            
            mrstG.nodes.coords = reshape(nodecoords, d, [])';

            mrstG.type = 'generic';

            % Compute geometrical data
            mrstG = computeGeometry(mrstG);
            
            rock.poro = ones(mrstG.cells.num, 1);
            rock.perm = ones(mrstG.cells.num, 1);

            % Compute half transmissibilities
            hT = computeTrans(mrstG, rock);

            cells.centroids = reshape(mrstG.cells.centroids', [], 1);
            cells.volumes   = farea*mrstG.cells.volumes;
            
            faces.centroids = reshape(mrstG.faces.centroids', [], 1);
            faces.normals   = farea*reshape(mrstG.faces.normals', [], 1);
            faces.areas     = farea*mrstG.faces.areas;

            nodes.coords = nodecoords;
            
            tPFVgeometry.cells = cells; 
            tPFVgeometry.faces = faces;
            tPFVgeometry.nodes = nodes;
            tPFVgeometry.hT    = farea*hT;
            
            if  mrstG.griddim == 1
                tPFVgeometry.faceArea = farea;
            end

            tPFVgeometry = TwoPointFiniteVolumeGeometry(tPFVgeometry);
            
        end        

        function tPFVgeometry = computeTPFVgeometryAD(grid, nodecoords, faceArea)
        % Compute the geometrical properties using function that can handle AD input for nodecoords and faceArea
        % The argument faceArea is only needed in 1D case.
            
            tp = grid.topology;
            
            switch tp.griddim
              case 1
                if nargin < 3 || isempty(faceArea)
                    % we use the faceArea stored in current tPFVgeometry
                    faceArea = grid.tPFVgeometry.faceArea;
                end
                tPFVgeometry = computeTPFVgeometryAD_1D(grid, nodecoords, faceArea);
              case 3
                tPFVgeometry = computeTPFVgeometryAD_3D(grid, nodecoords);
              otherwise
                error('Grid dimension not implemented yet.')
            end
            
        end
        
        function grid = setupHelpers(grid)
        % set the helpers structure.

            tp = grid.topology;
            
            tbls = setupTables(tp, 'includetbls', {'intfacetbl', 'extfacetbl'});
            intfacetbl     = tbls.intfacetbl;
            cellintfacetbl = tbls.cellintfacetbl;
            cellfacetbl    = tbls.cellfacetbl;
            celltbl        = tbls.celltbl;
            facetbl        = tbls.facetbl;
            extfacetbl     = tbls.extfacetbl;

            sgn = ones(cellintfacetbl.num, 1);
            f = cellintfacetbl.get('faces');
            c = cellintfacetbl.get('cells');
            sgn(grid.topology.faces.neighbors(f, 2) == c) = -1;

            prod = TensorProd();
            prod.tbl1 = cellintfacetbl;
            prod.tbl2 = intfacetbl;
            prod.tbl3 = celltbl;
            prod.reducefds = {'faces'};
            prod = prod.setup;

            divM = SparseTensor();
            divM = divM.setFromTensorProd(sgn, prod);
            divM = divM.getMatrix();

            diffop.div  = divM;
            diffop.grad = -divM';
            
            map = TensorMap();
            map.fromTbl  = cellfacetbl;
            map.toTbl    = cellintfacetbl;
            map.mergefds = {'cells', 'faces'};
            map = map.setup();
            
            P = map.getMatrix();

            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = cellintfacetbl;
            map.mergefds = {'cells'};
            map = map.setup();

            D = map.getMatrix();
            
            map = TensorMap();
            map.fromTbl = cellintfacetbl;
            map.toTbl = intfacetbl;
            map.mergefds = {'faces'};
            map = map.setup();

            S = map.getMatrix();

            trans = struct('P', P, ...
                           'D', D, ...
                           'S', S);
            
            cellextfacetbl = crossIndexArray(cellfacetbl, extfacetbl, {'faces'});

            extfaces.faces = cellextfacetbl.get('faces');
            extfaces.cells = cellextfacetbl.get('cells');

            map = TensorMap();
            map.fromTbl  = cellfacetbl;
            map.toTbl    = cellextfacetbl;
            map.mergefds = {'cells', 'faces'};
            
            extfaces.halfTransParentInd = map.getDispatchInd();

            sgn = ones(cellextfacetbl.num, 1);
            f = extfaces.faces;
            c = extfaces.cells;
            sgn(grid.topology.faces.neighbors(f, 2) == c) = -1;

            extfaces.sgn = sgn;

            faceextfacemap = zeros(facetbl.num, 1);
            faceextfacemap(extfacetbl.get('faces')) = (1 : extfacetbl.num)';

            grid.helpers = struct('diffop'        , diffop        , ...
                                  'trans'         , trans         , ...
                                  'extfaces'      , extfaces      , ...
                                  'faceextfacemap', faceextfacemap);
        end

        function grid = setupCellFluxOperators(grid)
        % Helper operators for CellFluxNorm. Note : those have no AD compliant version for the moment.
            G = getMRSTgrid(grid);
            grid.cellFluxOperators = getCellFluxOperatorsAll(G);
            
        end

        function grid = setupADhelpers(grid)
        % Setup the AD helper structures matrixOperators and vectorHelpers.
        % Those are used to compute the geometrical properties in computeTPFVgeometryAD.
            
            tp = grid.topology;
            
            switch tp.griddim
              case 1
                grid = grid.setupADhelpers1D();
              case 3
                grid = grid.setupADhelpers3D();               
              otherwise
                error('Grid dimension not implemented yet.')
            end
            
        end
        

        function grid = setupADhelpers1D(grid)
        % AD helpers in 1D
            
            tp = grid.topology;

            tbls = setupTables(tp);

            celltbl     = tbls.celltbl;
            facetbl     = tbls.facetbl;
            cellfacetbl = tbls.cellfacetbl;

            nc = celltbl.num;
            
            cellfacetbl2 = sortIndexArray(cellfacetbl, {'cells', 'faces'});
            faces = cellfacetbl2.get('faces');
            faces = reshape(faces, 2, [])';
            
            cellface12tbl.cells = (1 : nc)';
            cellface12tbl.faces1 = faces(:, 1);
            cellface12tbl.faces2 = faces(:, 2);
            cellface12tbl = IndexArray(cellface12tbl);
            
            % faceCentroids = nodeCoords

            prod = TensorProd();
            prod.tbl1 = cellfacetbl;
            prod.tbl2 = facetbl;
            prod.tbl3 = celltbl;
            prod.reducefds = {'faces'};
            prod = prod.setup();

            p = 0.5*ones(cellfacetbl.num, 1);

            tens = SparseTensor();
            tens = tens.setFromTensorProd(p, prod);

            matrixop.cellCentroids = tens.getMatrix();
            % To run:
            % cellCentroids = prod.eval(p, faceCentroids)

            %% Compute volumes(c, f1, f2) = abs(faceCentroids(f2) - faceCentroids(f1))*faceArea
            
            map1 = TensorMap();
            map1.fromTbl = facetbl;
            map1.toTbl = cellface12tbl;
            map1.replaceFromTblfds = {{'faces', 'faces1'}};
            map1.mergefds = {'faces1'};
            map1 = map1.setup();

            map2 = TensorMap();
            map2.fromTbl = facetbl;
            map2.toTbl = cellface12tbl;
            map2.replaceFromTblfds = {{'faces', 'faces2'}};
            map2.mergefds = {'faces2'};
            map2 = map2.setup();

            map3 = TensorMap();
            map3.fromTbl = cellface12tbl;
            map3.toTbl = celltbl;
            map3.mergefds = {'cells'};
            map3 = map3.setup();

            matrixop.vols1 = map1.getMatrix();
            matrixop.vols2 = map2.getMatrix();
            matrixop.vols3 = map3.getMatrix();
            % To run:
            % vols = map3.eval(map2.eval(faceCentroids) - map1.eval(faceCentroids))*faceArea

            %% Compute cellfacedist(c, f) = faceCentroids(f) - cellCentroids(c)
            
            map1 = TensorMap();
            map1.fromTbl = celltbl;
            map1.toTbl = cellfacetbl;
            map1.mergefds = {'cells'};
            map1 = map1.setup();
                        
            map2 = TensorMap();
            map2.fromTbl = facetbl;
            map2.toTbl = cellfacetbl;
            map2.mergefds = {'faces'};
            map2 = map2.setup();

            matrixop.cellfacedist1 = map1.getMatrix();
            matrixop.cellfacedist2 = map2.getMatrix();
            % To run:
            % cellfacedist = abs(map2.eval(faceCentroids) - map1.eval(cellCentroids))
            % hT = 1./cellfacedist*faceArea

            grid.matrixOperators = matrixop;
            
        end
        
        function grid = setupADhelpers3D(grid)
        % AD helpers in 3D
            
            tp = grid.topology;

            tbls = setupTables(tp, 'includetbls', {'vectbl', 'facenode12tbl'});

            celltbl       = tbls.celltbl;
            facetbl       = tbls.facetbl;
            vectbl        = tbls.vectbl;
            nodetbl       = tbls.nodetbl;
            cellfacetbl   = tbls.cellfacetbl;
            facenodetbl   = tbls.facenodetbl;
            facenode12tbl = tbls.facenode12tbl;
            
            % Triangles on the faces : facenode12tbl

            % sub-tetrahedra

            cellfacenode12tbl = crossIndexArray(cellfacetbl, facenode12tbl, {'faces'});

            % Variants that include spatial index

            facevectbl     = tbls.facevectbl;
            nodevectbl     = tbls.nodevectbl;
            cellvectbl     = tbls.cellvectbl;
            cellfacevectbl = tbls.cellfacevectbl;
            facenodevectbl = tbls.facenodevectbl;
            
            facenode12vectbl     = crossIndexArray(facenode12tbl    , vectbl       , {}, 'optpureproduct', true);
            cellfacenode12vectbl = crossIndexArray(cellfacenode12tbl, vectbl       , {}, 'optpureproduct', true);
            
            % Setup pseudo face centroids : pseudoFaceCentroids(f) = noodeCoords(n, f)/number(n in f)

            map = TensorMap();
            map.fromTbl = nodevectbl;
            map.toTbl = facenodevectbl;
            map.mergefds = {'nodes', 'vec'};
            map = map.setup();

            matrixop.nodecoords_facenode = map.getMatrix();
            % To run:
            % nodecoords_facenode = map.eval(nodeCoords);  

            map = TensorMap();
            map.fromTbl = facenodevectbl;
            map.toTbl = facevectbl;
            map.mergefds = {'faces', 'vec'};
            map = map.setup();

            matrixop.facevec = map.getMatrix();
            % To run:
            % facevec = map.eval(nodecoords_facenode); % to be modified later
            
            map = TensorMap();
            map.fromTbl = facenodetbl;
            map.toTbl = facetbl;
            map.mergefds = {'faces'};
            map = map.setup();

            nnode_per_face = map.eval(ones(facenodetbl.num, 1));

            map = TensorMap();
            map.fromTbl = facetbl;
            map.toTbl = facevectbl;
            map.mergefds = {'faces'};
            map = map.setup();

            vechelps.nnode_per_face = map.eval(nnode_per_face);

            % To run:
            % pseudoFaceCentroids = facevec./nnode_per_face;

            %% Setup nfVector = nodeCoords - faceCentroids,  in triangle, belongs to nodefacevec

            map1          = TensorMap();
            map1.fromTbl  = facevectbl;
            map1.toTbl    = facenodevectbl;
            map1.mergefds = {'faces', 'vec'};
            map1          = map1.setup();

            map2          = TensorMap();
            map2.fromTbl  = nodevectbl;
            map2.toTbl    = facenodevectbl;
            map2.mergefds = {'nodes', 'vec'};
            map2          = map2.setup();

            matrixop.nfVector_1 = map1.getMatrix();
            matrixop.nfVector_2 = map2.getMatrix();
            % To Run:
            % nfVector = map2.eval(nodeCoords) - map1.eval(pseudoFaceCentroids); % belongs to facenodevectbl

            %% Setup normals, centers and areas of the triangles

            %% Setup the cross-product operator

            val = [2, 3, 1,  1
                   3, 2, 1, -1
                   1, 3, 2, -1
                   3, 1, 2,  1
                   1, 2, 3,  1
                   2, 1, 3, -1
                  ];

            permvec123tbl.vec1 = val(:, 1);
            permvec123tbl.vec2 = val(:, 2);
            permvec123tbl.vec3 = val(:, 3);
            permvec123tbl = IndexArray(permvec123tbl);

            sigma = val(:, 4); % belongs to permvec123tbl

            %% We start with the normals of the triangle : triNormals(n1, n2, f, i) = 1/2*(sum(sigma(i, j, k)*nfVector(n1, f, i)*nfVector(n2, f, j))

            map = TensorMap();
            map.fromTbl = facenodevectbl;
            map.toTbl = facenode12vectbl;
            map.replaceFromTblfds = {{'nodes', 'nodes1'}};
            map.mergefds = {'faces', 'nodes1', 'vec'};
            map = map.setup();

            matrixop.nfVector1 = map.getMatrix();
            % To run:
            % nfVector1 = map.eval(nfVector);
            
            map = TensorMap();
            map.fromTbl = facenodevectbl;
            map.toTbl = facenode12vectbl;
            map.replaceFromTblfds = {{'nodes', 'nodes2'}};
            map.mergefds = {'faces', 'nodes2', 'vec'};
            map = map.setup();

            matrixop.nfVector2 = map.getMatrix();
            % To run:
            % nfVector2 = map.eval(nfVector);

            prod = TensorProd();
            prod.tbl1 = permvec123tbl;
            prod.tbl2 = facenode12vectbl;
            prod.replacefds1 = {{'vec1', 'vec'}};
            prod.reducefds = {'vec'};
            prod = prod.setup();

            facenode12vec13tbl = prod.tbl3;

            M = SparseTensor();
            matrixop.triNormals1 = M.setFromTensorProd(sigma, prod);
            matrixop.triNormals1 = matrixop.triNormals1.getMatrix();
            % To Run:
            % triNormals = prod.eval(sigma, nfVector1); % temporary value

            prod = TensorProd();
            prod.tbl1 = facenode12vec13tbl;
            prod.tbl2 = facenode12vectbl;
            prod.tbl3 = facenode12vectbl;
            prod.replacefds1 = {{'vec3', 'vec'}};
            prod.replacefds2 = {{'vec', 'vec2'}};
            prod.reducefds = {'vec2'};
            prod.mergefds = {'faces', 'nodes1', 'nodes2'};
            prod = prod.setup();

            matrixop.triNormals = prod.getMatrices();
            % To run:
            % triNormals = 0.5*prod.eval(triNormals, nfVector2);

            %% We compute now the areas of the triangle : triAreas(n1, n2, f) = ( sum ( triNormals(n1, n2, f, i)*triNormals(n1, n2, f, i) ) )^(1/2)

            prod = TensorProd();
            prod.tbl1 = facenode12vectbl;
            prod.tbl2 = facenode12vectbl;
            prod.tbl3 = facenode12tbl;
            prod.reducefds = {'vec'};
            prod.mergefds = {'faces', 'nodes1', 'nodes2'};
            prod = prod.setup();

            matrixop.triAreas = prod.getMatrices();
            % To run:
            % triAreas = prod.eval(triNormals, triNormals);
            % triAreas = triAreas.^0.5;

            %% We compute the centroid of the triangle

            map = TensorMap();
            map.fromTbl = nodevectbl;
            map.toTbl = facenode12vectbl;
            map.replaceFromTblfds = {{'nodes', 'nodes1'}};
            map.mergefds = {'nodes1', 'vec'};
            map = map.setup();

            matrixop.n1Coords = map.getMatrix();
            % To run:
            % n1Coords = map.eval(nodeCoords);
            
            map = TensorMap();
            map.fromTbl = nodevectbl;
            map.toTbl = facenode12vectbl;
            map.replaceFromTblfds = {{'nodes', 'nodes2'}};
            map.mergefds = {'nodes2', 'vec'};
            map = map.setup();

            matrixop.n2Coords = map.getMatrix();
            % To run:
            % n2Coords = map.eval(nodeCoords);
            
            map = TensorMap();
            map.fromTbl = facevectbl;
            map.toTbl = facenode12vectbl;
            map.mergefds = {'faces', 'vec'};
            map = map.setup();

            matrixop.pseudoFaceCentroids = map.getMatrix();
            % To run:
            % pseudoFaceCentroids = map.eval(pseudoFaceCentroids);
            % and
            % triCentroids = (n1Coords + n2Coords + pseudoFaceCentroids)/3;

            %% We compute faceNormals(f, i) = sum( triNormals(n1, n2, f, i))

            map = TensorMap();
            map.fromTbl = facenode12vectbl;
            map.toTbl = facevectbl;
            map.mergefds = {'faces', 'vec'};
            map = map.setup();

            matrixop.faceNormals = map.getMatrix();
            % To run:
            % faceNormals = map.eval(triNormals);
            
            %% We compute faceAreas(f) = sum( triAreas(n1, n2, f))

            map = TensorMap();
            map.fromTbl = facenode12tbl;
            map.toTbl = facetbl;
            map.mergefds = {'faces'};
            map = map.setup();

            matrixop.faceAreas = map.getMatrix();
            % To run:
            % faceAreas = map.eval(triAreas);


            %% We compute the face centroids faceCentroids(f, i) = triCentroids(n1, n2, f, i)*triAreas(n1, n2, f)/faceAreas(f)

            prod = TensorProd();
            prod.tbl1 = facenode12vectbl;
            prod.tbl2 = facenode12tbl;
            prod.tbl3 = facevectbl;
            prod.mergefds = {'faces'};
            prod.reducefds = {'nodes1', 'nodes2'};
            prod = prod.setup();

            map = TensorMap();
            map.fromTbl = facetbl;
            map.toTbl = facevectbl;
            map.mergefds = {'faces'};
            map = map.setup();

            matrixop.faceCentroids1 = prod.getMatrices();
            matrixop.faceCentroids2 = map.getMatrix();
            % To run:
            % faceCentroids = prod.eval(triCentroids, triAreas)./map.eval(faceAreas);
            
            %% We compute pseudoCellCentroids(c, i) = sum(faceCentroids(f, i))/number_faces_per_cell

            map = TensorMap();
            map.fromTbl = facevectbl;
            map.toTbl = cellfacevectbl; 
            map.mergefds = {'faces', 'vec'};
            map = map.setup();

            matrixop.pseudoCellCentroids1 = map.getMatrix();
            % To run:
            % pseudoCellCentroids = map.eval(faceCentroids); % temporary value
            
            map = TensorMap();
            map.fromTbl = cellfacevectbl; 
            map.toTbl = cellvectbl;
            map.mergefds = {'cells', 'vec'};
            map = map.setup();

            matrixop.pseudoCellCentroids = map.getMatrix();
            % To run:
            % pseudoCellCentroids = map.eval(pseudoCellCentroids); % temporary value
            
            number_faces_per_cell = ones(cellfacevectbl.num, 1); % temporary values
            number_faces_per_cell = map.eval(number_faces_per_cell);

            vechelps.number_faces_per_cell = number_faces_per_cell;
            % To run:
            % pseudoCellCentroids = pseudoCellCentroids./number_faces_per_cell;

            %% We compute the scalar product : scalprod(f, c) = sum((faceCentroids(f, i) - pseudoCellCentroids(c, i))*faceNormals(f, i))

            map1 = TensorMap();
            map1.fromTbl = facevectbl;
            map1.toTbl = cellfacevectbl;
            map1.mergefds = {'faces', 'vec'};
            map1 = map1.setup();

            map2 = TensorMap();
            map2.fromTbl = cellvectbl;
            map2.toTbl = cellfacevectbl;
            map2.mergefds = {'cells', 'vec'};
            map2 = map2.setup();

            matrixop.scalprod1 = map1.getMatrix();
            matrixop.scalprod2 = map2.getMatrix();
            % To run
            % scalprod = (map1.eval(faceCentroids) - map2.eval(pseudoCellCentroids)).*map1.eval(faceNormals); % temporary value

            map = TensorMap();
            map.fromTbl = cellfacevectbl;
            map.toTbl = cellfacetbl;
            map.mergefds = {'cells', 'faces'};
            map = map.setup();

            matrixop.scalprod = map.getMatrix();
            % To run:
            % scalprod = map.eval(scalprod);
            % and
            % cellfacesign = (scalprod > 0) - (scalprod < 0);

            %% We compute the volumes of the tetrahedra :
            %% tetraVolumes(c, n1, n2, f)  = 1/3*sum((triCentroids(n1, n2, f, i) - pseudoCellCentroids(c, i))*cellfacesign(c, f)*triNormals(n1, n2, f, i))

            map1 = TensorMap();
            map1.fromTbl = facenode12vectbl;
            map1.toTbl = cellfacenode12vectbl;
            map1.mergefds = {'faces', 'nodes1', 'nodes2', 'vec'};
            map1 = map1.setup();

            map2 = TensorMap();
            map2.fromTbl = cellvectbl;
            map2.toTbl = cellfacenode12vectbl;
            map2.mergefds = {'cells', 'vec'};
            map2 = map2.setup();

            map3 = TensorMap();
            map3.fromTbl = cellfacetbl;
            map3.toTbl = cellfacenode12vectbl;
            map3.mergefds = {'cells', 'faces'};
            map3 = map3.setup();

            matrixop.tetraVolumes1 = map1.getMatrix();
            matrixop.tetraVolumes2 = map2.getMatrix();
            matrixop.tetraVolumes3 = map3.getMatrix();
            % To run:
            % tetraVolumes = 1/3*(map1.eval(triCentroids) - map2.eval(pseudoCellCentroids)).*map1.eval(triNormals).*map3.eval(cellfacesign);

            map = TensorMap();
            map.fromTbl = cellfacenode12vectbl;
            map.toTbl = cellfacenode12tbl;
            map.mergefds = {'cells', 'faces', 'nodes1', 'nodes2'};
            map = map.setup();

            matrixop.tetraVolumes = map.getMatrix();
            % To run:
            % tetraVolumes = map.eval(tetraVolumes);
            
            % We compute cell volumes cellVolumes(c) = sum( tetraVolumes(c, f, n1, n2))

            map = TensorMap();
            map.fromTbl = cellfacenode12tbl;
            map.toTbl = celltbl;
            map.mergefds = {'cells'};
            map = map.setup();

            matrixop.cellVolumes = map.getMatrix();
            % To run:
            % cellVolumes = map.eval(tetraVolumes);
            
            %% We adjust the centroids with relativeTetraCentroids(c, f, n1, n2, i) = 3/4*sum( triCentroids(f, n1, n2, i) - pseudoCellCentroids(c, i) )

            map1 = TensorMap();
            map1.fromTbl = facenode12vectbl;
            map1.toTbl = cellfacenode12vectbl;
            map1.mergefds = {'faces', 'nodes1', 'nodes2', 'vec'};
            map1 = map1.setup();

            map2 = TensorMap();
            map2.fromTbl = cellvectbl;
            map2.toTbl = cellfacenode12vectbl;
            map2.mergefds = {'cells', 'vec'};
            map2 = map2.setup();

            matrixop.relativeTetraCentroids1 = map1.getMatrix();
            matrixop.relativeTetraCentroids2 = map2.getMatrix();
            % To run:
            % relativeTetraCentroids = 3/4*(map1.eval(triCentroids) - map2.eval(pseudoCellCentroids));
            
            %% We define cellCentroids(c, i) = pseudoCellCentroids(c, i) + sum(relativeTetraCentroids(c, f, n1, n2, i)*tetraVolumes(c, f, n1, n2))/V(c)

            map1 = TensorMap();
            map1.fromTbl = cellfacenode12tbl;
            map1.toTbl = cellfacenode12vectbl;
            map1.mergefds = {'cells', 'faces', 'nodes1', 'nodes2'};
            map1 = map1.setup();

            map2 = TensorMap();
            map2.fromTbl = cellfacenode12vectbl;
            map2.toTbl = cellvectbl;
            map2.mergefds = {'cells', 'vec'};
            map2 = map2.setup();

            map3 = TensorMap();
            map3.fromTbl = celltbl;
            map3.toTbl = cellvectbl;
            map3.mergefds = {'cells'};
            map3 = map3.setup();

            matrixop.cellCentroids1 = map1.getMatrix();
            matrixop.cellCentroids2 = map2.getMatrix();
            matrixop.cellCentroids3 = map3.getMatrix();
            % To run:
            % cellCentroids = pseudoCellCentroids + map2.eval(relativeTetraCentroids.*map1.eval(tetraVolumes))./map3.eval(cellVolumes);

            %% Compute the half-transmissibility

            %% We compute d(c,f,i) = faceCentroids(f,i) - cellCentroids(c,i)

            map1 = TensorMap();
            map1.fromTbl = facevectbl;
            map1.toTbl = cellfacevectbl;
            map1.mergefds = {'faces', 'vec'};
            map1 = map1.setup();

            map2 = TensorMap();
            map2.fromTbl = cellvectbl;
            map2.toTbl = cellfacevectbl;
            map2.mergefds = {'cells', 'vec'};
            map2 = map2.setup();

            matrixop.d1 = map1.getMatrix();
            matrixop.d2 = map2.getMatrix();
            % To run:
            % d = map1.eval(faceCentroids) - map2.eval(cellCentroids)


            %% We compute the signed normals n(c, f, i) = cellfacesign(c, f)*faceNormals(f, i)
            
            prod = TensorProd();
            prod.tbl1 = cellfacetbl;
            prod.tbl2 = facevectbl;
            prod.tbl3 = cellfacevectbl;
            prod.mergefds = {'faces'};
            prod = prod.setup();

            matrixop.cellfaceNormals = prod.getMatrices();
            % To run:
            % cellfaceNormals = prod.eval(cellfacesign, faceNormals);
            
            %% We compute:
            %% dscaln(c, f) = d(c, f, i)*cellfaceNormals(c, f, i),
            %% dsq(c, f)    = d(c, f, i)*d(c, f, i),
            %% hT(c, f)     = dKn(c, f)/dsq(c, f)
            
            prod = TensorProd();
            prod.tbl1 = cellfacevectbl;
            prod.tbl2 = cellfacevectbl;
            prod.tbl3 = cellfacetbl;
            prod.mergefds = {'cells', 'faces'};
            prod.reducefds = {'vec'};
            prod = prod.setup();

            matrixop.dscaln = prod.getMatrices();
            matrixop.dsq    = matrixop.dscaln; % for simplicity
                                               % To run:
                                               % dscaln = prod.eval(d, cellfaceNormals);
                                               % dsq    = prod.eval(d, d);
                                               % hT     = dscaln./dsq;

            grid.matrixOperators = matrixop;
            grid.vectorHelpers   = vechelps;
            
        end

        
        function tPFVgeometry = computeTPFVgeometryAD_1D(grid, nodecoords, faceArea)
        % version of computeTPFVgeometryAD in 1D

            if isempty(faceArea)
                faceArea = 1;
            end
            
            m  = grid.matrixOperators;
            tp = grid.topology;

            nf = tp.faces.num;

            faceCentroids = nodecoords;

            cellCentroids = m.cellCentroids*faceCentroids;
            
            cellVolumes = m.vols3*(abs(m.vols2*faceCentroids - m.vols1*faceCentroids))*faceArea;
            
            tPFVgeometry.cells.volumes   = cellVolumes;
            tPFVgeometry.cells.centroids = cellCentroids;
            tPFVgeometry.faces.areas     = faceArea*ones(nf, 1);
            tPFVgeometry.faces.normals   = faceArea*ones(nf, 1);
            tPFVgeometry.faces.centroids = faceCentroids;
            tPFVgeometry.nodes.coords    = nodecoords;

            cellfacedist = abs(m.cellfacedist2*faceCentroids - m.cellfacedist1*cellCentroids);

            tPFVgeometry.hT = 1./cellfacedist*faceArea;

            tPFVgeometry.faceArea = faceArea;
            
            tPFVgeometry = TwoPointFiniteVolumeGeometry(tPFVgeometry);

        end

        function tPFVgeometry = computeTPFVgeometryAD_3D(grid, nodecoords)
        % version of computeTPFVgeometryAD in 1D
            
            tp = grid.topology;
            m  = grid.matrixOperators;
            v  = grid.vectorHelpers;

            n = tp.griddim;
            
            nodecoords_facenode = m.nodecoords_facenode*nodecoords;
            facevec             = m.facevec*nodecoords_facenode;
            pseudoFaceCentroids = facevec./v.nnode_per_face;
            
            nfVector  = m.nfVector_2*nodecoords - m.nfVector_1*pseudoFaceCentroids; 

            nfVector1 = m.nfVector1*nfVector;
            nfVector2 = m.nfVector2*nfVector;
            
            triNormals = m.triNormals1*nfVector1;
            triNormals = 0.5*grid.applyProduct(m.triNormals, triNormals, nfVector2);

            triAreas = grid.applyProduct(m.triAreas, triNormals, triNormals);
            triAreas = triAreas.^0.5;
            
            n1Coords       = m.n1Coords*nodecoords;
            n2Coords       = m.n2Coords*nodecoords;
            pseudoFaceCentroids = m.pseudoFaceCentroids*pseudoFaceCentroids;

            triCentroids = (n1Coords + n2Coords + pseudoFaceCentroids)/3;

            faceNormals = m.faceNormals*triNormals;
            faceAreas   = m.faceAreas*triAreas;

            D1 = m.faceCentroids1.D1;
            D2 = m.faceCentroids1.D2;
            S  = m.faceCentroids1.S;            

            faceCentroids = (S*((D1*triCentroids).*(D2*triAreas)))./(m.faceCentroids2*faceAreas);
            
            pseudoCellCentroids = m.pseudoCellCentroids1*faceCentroids;
            pseudoCellCentroids = m.pseudoCellCentroids*pseudoCellCentroids;
            pseudoCellCentroids = pseudoCellCentroids./v.number_faces_per_cell;

            scalprod = (m.scalprod1*faceCentroids - m.scalprod2*pseudoCellCentroids).*(m.scalprod1*faceNormals);
            scalprod = m.scalprod*scalprod;
            scalprod = value(scalprod);
            cellfacesign = (scalprod > 0) - (scalprod < 0);
            
            tetraVolumes           = 1/3*(m.tetraVolumes1*triCentroids - m.tetraVolumes2*pseudoCellCentroids).*(m.tetraVolumes1*triNormals).*(m.tetraVolumes3*cellfacesign);
            tetraVolumes           = m.tetraVolumes*tetraVolumes;
            cellVolumes            = m.cellVolumes*tetraVolumes;
            relativeTetraCentroids = 3/4*(m.relativeTetraCentroids1*triCentroids - m.relativeTetraCentroids2*pseudoCellCentroids);
            
            cellCentroids = pseudoCellCentroids + (m.cellCentroids2*(relativeTetraCentroids.*(m.cellCentroids1*tetraVolumes)))./(m.cellCentroids3*cellVolumes);

            tPFVgeometry.cells.volumes   = cellVolumes;
            tPFVgeometry.cells.centroids = cellCentroids;
            tPFVgeometry.faces.areas     = faceAreas;
            tPFVgeometry.faces.normals   = faceNormals;
            tPFVgeometry.faces.centroids = faceCentroids;
            tPFVgeometry.nodes.coords    = nodecoords;

            % Computation of transmissibilities
            
            d = m.d1*faceCentroids - m.d2*cellCentroids;

            cellfaceNormals = grid.applyProduct(m.cellfaceNormals, cellfacesign, faceNormals);
            
            dscaln = grid.applyProduct(m.dscaln, d, cellfaceNormals);
            dsq    = grid.applyProduct(m.dsq, d, d);

            tPFVgeometry.hT = dscaln./dsq;
            
            tPFVgeometry = TwoPointFiniteVolumeGeometry(tPFVgeometry);
            
        end

        function G = getMRSTgrid(grid)
        % Returns MRST grid that can be used for plotting
            
            G  = grid.topology;
            tg = grid.tPFVgeometry;

            d = G.griddim;

            G.type = 'generic';

            G.cells.centroids = reshape(tg.cells.centroids, d, [])';
            G.cells.volumes   = tg.cells.volumes;
            G.faces.centroids = reshape(tg.faces.centroids, d, [])';
            G.faces.normals   = reshape(tg.faces.normals, d, [])';
            G.faces.areas     = tg.faces.areas;
            G.nodes.coords    = reshape(tg.nodes.coords, d, [])';
            
        end


        function vols = getVolumes(grid)
            
            vols = grid.tPFVgeometry.cells.volumes;
            
        end

        function areas = getFaceAreas(grid)
            
            areas = grid.tPFVgeometry.faces.areas;
            
        end
        
        function u = getHarmFace(grid, c)
        % Returns fluxes for each internal faces for the cell-valued vector c

            op = grid.helpers.trans;
            hT = grid.tPFVgeometry.hT;
            
            u = 1 ./ (op.S * ( 1 ./ (op.D*c .* op.P*hT)));
            
        end
        
        function [bchT, bccells, bcsgn] = getBcHarmFlux(grid, u, bcfaces)
        % Returns half transmissibilities and cell indexing for the given boundary faces

            hT   = grid.tPFVgeometry.hT;
            exf  = grid.helpers.extfaces;

            extfaceind = grid.helpers.faceextfacemap(bcfaces);
            
            bccells = exf.cells(extfaceind);
            bcsgn   = exf.sgn(extfaceind);
            bchT    = hT(exf.halfTransParentInd(extfaceind));
            
        end

        function jsq = getCellFluxNorm(grid, u)
    
            P = grid.cellFluxOperators.P;
            S = grid.cellFluxOperators.S;
    
            j = P*u;
            jsq = j.^2;
            jsq = S*jsq;
            
        end

        function nc = getNumberOfCells(grid)

            nc = grid.topology.cells.num;
            
        end        
    end

    methods(Static)

        function value = applyProduct(matrixOperator, arg1, arg2)
        % Function that assembly product structure (two dispatches followed by one reduction)
            
            D1 = matrixOperator.D1;
            D2 = matrixOperator.D2;
            S  = matrixOperator.S;

            value = S*((D1*arg1).*(D2*arg2));
            
        end
    end    
    
end

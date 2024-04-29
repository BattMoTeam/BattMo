classdef StokesSolver

% Setup stokes solver    
    
%% Reference
% @inbook{da_Veiga_2014,
%         title={Foundations of mimetic finite difference method},
%         ISBN={9783319026633},
%         url={http://dx.doi.org/10.1007/978-3-319-02663-3_2},
%         DOI={10.1007/978-3-319-02663-3_2},
%         booktitle={The Mimetic Finite Difference Method for Elliptic Problems},
%         publisher={Springer International Publishing},
%         author={da Veiga,
%         Lourenço Beirão and Lipnikov,
%         Konstantin and Manzini,
%         Gianmarco},
%         year={2014},
%         pages={41–65} }

    properties

        G
        
        viscosity % (scalar value)
        
        dim = 2 % dimension of the grid (only 2d support for now)
        
        dirichletBc % Dirichlet boundary specification struct with fields
                    %   - nodes
                    %   - faces
        
        neumannFaces   % Index of the Dirichlet boundary faces
        
        operators % struct with fields
                  % - A    : System matrix (see definition below)
                  % - M    : Stiffness matrix
                  % - Div  : Divergence operator
                  % - Flux : Flux operator
                  % - Diri : Dirichlet operator matrix (used in system matrix and to setup dirichletValues)
                  % - Neum : Neuman operator to setup external viscous forces term.
                  % 
                  % The Stokes system takes the form
                  %
                  %      A*x = b
                  %
                  % with
                  %
                  %      [u 
                  %  x =  p 
                  %       l]
                  %
                  %      [ M     Div' Diri'      
                  %  A=    Div   0    0          
                  %        Diri  0    0    ]     
                  %
                  %      [viscforces      
                  %  b =  source          
                  %       dirichletValues]
                  %
                  % where
                  %        u : velocities degrees of freedom (dim per node and 1 per face)
                  %        p : pressure degrees of freedom (1 per cell)
                  %        l : pseudo-lagrangian multiplier to include Dirichlet conditions in a generic way (at cost to extending the syste)
                  %
                  %        M    : Stiffness matrich
                  %        Div  : Divergence operator
                  %        Diri : Diri matrix
                  %        
                  %        viscforces      : vector of same size of velocity u using the same degree of freedom.
                  %                          It corresponds by duality to the viscous forces occuring at that the Neuman
                  %                          boundary faces. It is setup using the Neum operator
                  %        source          : source term with one value per cell
                  %        dirichletValues : They are given through the degrees of freedom (same as velocity u) for the
                  %                          Dirichlet boundary faces. We use Diri to setup those

        
        tbls % Index array tables generated in setup process
             %
             % - cellNodeTbl
             % - vecTbl
             % - nodeTbl
             % - nodeVecTbl
             % - cellTbl
             % - faceTbl
             % - cellVecTbl
             % - faceVecTbl
             % - cellFaceTbl
             % - cellFaceVecTbl
             % - faceNodeVecTbl
             % - faceNodeTbl
             % - extFaceTbl
             % - cellGtypeGindTbl
             % - nodeVecGtypeGindTbl
             % - faceGtypeGindTbl

        
        helpers % struct with fields
                % - R : Matrix as in eq 8.35 in Reference
                % - N : Matrix as in eq 8.33 in Reference
                % - dirOp : operators for the dirichlet boundary conditions with fields
                %      - gPn : global full-dof to nodal-vector-dof at Dirichlet
                %      - gPf : global full-dof to face-dof at Dirichlet
                %      - dPn : dirichlet from full-dof to nodal-vector-dof 
                %      - dPf : dirichlet from full-dof to face-dof
    end

    methods

        function solver = StokesSolver(G, dirichletBc, neumannFaces, viscosity)

            solver.G            = G;
            solver.dirichletBc  = dirichletBc;
            solver.neumannFaces = neumannFaces;
            solver.viscosity    = viscosity;
            
            dim = solver.dim;
            
            tbls = setupTables(G, 'includetbls', {'cellnodetbl', 'vectbl', 'extfacetbl'});

            cellNodeTbl    = tbls.cellnodetbl;
            vecTbl         = tbls.vectbl;
            nodeTbl        = tbls.nodetbl;
            nodeVecTbl     = tbls.nodevectbl;
            cellTbl        = tbls.celltbl;
            faceTbl        = tbls.facetbl;
            cellVecTbl     = tbls.cellvectbl;
            faceVecTbl     = tbls.facevectbl;
            cellFaceTbl    = tbls.cellfacetbl;
            cellFaceVecTbl = tbls.cellfacevectbl;
            faceNodeVecTbl = tbls.facenodevectbl;
            faceNodeTbl    = tbls.facenodetbl;
            extFaceTbl     = tbls.extfacetbl;

            cellNodeVecTbl = crossIndexArray(cellNodeTbl, vecTbl, {});
            cellNodeVecTbl = sortIndexArray(cellNodeVecTbl, {'cells', 'vec', 'nodes'});

            nodeVecGindTbl          = nodeVecTbl.addInd('gind', (1 : nodeVecTbl.num)');
            nodeVecGtypeGindTbl     = nodeVecGindTbl.addInd('gtype', ones(nodeVecGindTbl.num, 1));
            cellNodeVecGtypeGindTbl = crossIndexArray(cellNodeTbl, nodeVecGtypeGindTbl, {'nodes'});

            faceGindTbl          = faceTbl.addInd('gind', (1 : faceTbl.num)');
            faceGtypeGindTbl     = faceGindTbl.addInd('gtype', 2*ones(faceGindTbl.num, 1));
            cellFaceGtypeGindTbl = crossIndexArray(cellFaceTbl, faceGtypeGindTbl, {'faces'});

            % table for our degrees of freedom
            gtypeGindTbl1 = projIndexArray(nodeVecGtypeGindTbl, {'gtype', 'gind'});
            gtypeGindTbl2 = projIndexArray(faceGtypeGindTbl, {'gtype', 'gind'});
            gtypeGindTbl = concatIndexArray(gtypeGindTbl1, gtypeGindTbl2, {});

            % cell decompositions
            cellGtypeGindTbl1 = projIndexArray(cellNodeVecGtypeGindTbl, {'cells', 'gtype', 'gind'});
            cellGtypeGindTbl2 = projIndexArray(cellFaceGtypeGindTbl, {'cells', 'gtype', 'gind'});
            cellGtypeGindTbl = concatIndexArray(cellGtypeGindTbl1, cellGtypeGindTbl2, {});
            cellGtypeGindTbl = sortIndexArray(cellGtypeGindTbl, {'cells', 'gtype', 'gind'});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %% Assembly of the stiffness matrix M %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Assembly of matrix N (see after eq 8.33 in ref)

            clear polTbl
            polTbl.pol = (1 : dim + 1)';
            polTbl = IndexArray(polTbl);

            cellNodePolTbl = crossIndexArray(cellNodeTbl, polTbl, {});

            nodeVec = reshape(G.nodes.coords', [], 1);
            cellVec = reshape(G.cells.centroids', [], 1);

            map1 = TensorMap;
            map1.fromTbl  = nodeVecTbl;
            map1.toTbl    = cellNodeVecTbl;
            map1.mergefds = {'nodes', 'vec'};
            map1 = map1.setup();

            map2 = TensorMap;
            map2.fromTbl  = cellVecTbl;
            map2.toTbl    = cellNodeVecTbl;
            map2.mergefds = {'cells', 'vec'};
            map2 = map2.setup();

            cellNodeVec = map1.eval(nodeVec) - map2.eval(cellVec);

            clear linTbl
            linTbl.lin = (1 : dim*(dim + 1))';
            linTbl = IndexArray(linTbl);

            clear nk_linTbl
            nk_linTbl.lin = (dim*(dim + 1)/2 + 1 : linTbl.num)';
            nk_linTbl = IndexArray(nk_linTbl);

            cellLinTbl     = crossIndexArray(cellTbl, linTbl, {});
            nk_cellLinTbl = crossIndexArray(cellTbl, nk_linTbl, {});

            A = [[1 1 1  1; ...
                  2 2 1  1; ...
                  3 1 3  1; ...
                  3 2 2 -1; ...
                  4 1 2  1; ...
                  5 2 3  1; ...
                  6 1 3  1; ...
                  6 2 2  1]];

            clear linVecPolTbl
            linVecPolTbl.lin =  A(:, 1);
            linVecPolTbl.vec =  A(:, 2);
            linVecPolTbl.pol =  A(:, 3);
            linVecPolTbl = IndexArray(linVecPolTbl);

            linVecPol = A(:, 4);

            A = [[1 2;
                  2 3]];

            clear vecPolTbl
            vecPolTbl.vec = A(:, 1);
            vecPolTbl.pol = A(:, 2);
            vecPolTbl = IndexArray(vecPolTbl);

            prod = TensorProd();
            prod.tbl1 = vecPolTbl;
            prod.tbl2 = linVecPolTbl;
            prod.replacefds1 = {{'vec', 'vec2'}};
            prod.replacefds2 = {{'vec', 'vec1'}};
            prod.reducefds = {'pol'};
            prod = prod.setup();

            delta = ones(vecPolTbl.num, 1);
            firstLinVec12    = prod.eval(delta, linVecPol);
            firstLinVec12Tbl = prod.tbl3;

            cellNodeVecLinTbl = crossIndexArray(cellNodeVecTbl, linTbl, {});

            prod = TensorProd();
            prod.tbl1 = firstLinVec12Tbl;
            prod.tbl2 = cellNodeVecTbl;
            prod.tbl3 = cellNodeVecLinTbl;
            prod.replacefds1 = {{'vec1', 'vec'}};
            prod.replacefds2 = {{'vec', 'vec2'}};
            prod.reducefds = {'vec2'};
            prod = prod.setup();

            cellNodeVecLin1 = prod.eval(firstLinVec12, cellNodeVec);

            clear constPolTbl
            constPolTbl.pol = 1;
            constPolTbl = IndexArray(constPolTbl);

            constLinVecTbl = crossIndexArray(constPolTbl, linVecPolTbl, {'pol'});
            constLinVecTbl = projIndexArray(constLinVecTbl, {'lin', 'vec'});

            prod = TensorProd();
            prod.tbl1 = linVecPolTbl;
            prod.tbl2 = constPolTbl;
            prod.tbl3 = constLinVecTbl;
            prod.reducefds = {'pol'};
            prod = prod.setup();

            delta = ones(constPolTbl.num, 1);

            constLinVec = prod.eval(linVecPol, delta);

            map = TensorMap();
            map.fromTbl = constLinVecTbl;
            map.toTbl = cellNodeVecLinTbl;
            map.mergefds = {'lin', 'vec'};
            map = map.setup();

            cellNodeVecLin2 = map.eval(constLinVec);

            cellNodeVecLin = cellNodeVecLin1 + cellNodeVecLin2;

            [cellNodeVecLinGtypeGindTbl, indstruct] = crossIndexArray(cellNodeVecLinTbl, nodeVecGtypeGindTbl, {'nodes', 'vec'});

            direct = true;
            if direct
                cellNodeVecLinGtypeGind = cellNodeVecLin(indstruct{1}.inds);
            else
                map = TensorMap();
                map.fromTbl = cellNodeVecLinTbl;
                map.toTbl = cellNodeVecLinGtypeGindTbl;
                map.mergefds = {'cells', 'nodes', 'vec', 'lin'};
                map = map.setup();

                cellNodeVecLinGtypeGind = map.eval(cellNodeVecLin);
            end

            prod = TensorProd();
            prod.tbl1 = cellNodeVecLinGtypeGindTbl;
            prod.tbl2 = cellLinTbl;
            prod.tbl3 = cellGtypeGindTbl;
            prod.reducefds = {'lin'};
            prod.mergefds  = {'cells'};
            prod = prod.setup();

            tens = SparseTensor();
            tens = tens.setFromTensorProd(cellNodeVecLinGtypeGind, prod);

            N = tens.getMatrix();

            prod = TensorProd();
            prod.tbl1 = cellNodeVecLinGtypeGindTbl;
            prod.tbl2 = nk_cellLinTbl;
            prod.tbl3 = cellGtypeGindTbl;
            prod.reducefds = {'lin'};
            prod.mergefds  = {'cells'};
            prod = prod.setup();

            tens = SparseTensor();
            tens = tens.setFromTensorProd(cellNodeVecLinGtypeGind, prod);

            nk_N = tens.getMatrix();

            %% Setup of matrix R (see eq 8.35 in reference)

            faceVec = reshape(G.faces.normals', [], 1);

            cells = cellFaceTbl.get('cells');
            faces = cellFaceTbl.get('faces');
            sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1; % sgn is in cellfacetbl

            prod = TensorProd();
            prod.tbl1 = cellFaceTbl;
            prod.tbl2 = faceVecTbl;
            prod.tbl3 = cellFaceVecTbl;
            prod.mergefds = {'faces'};
            prod = prod.setup();

            cellFaceVec = prod.eval(sgn, faceVec);

            prod = TensorProd();
            prod.tbl1 = faceNodeTbl;
            prod.tbl2 = cellFaceVecTbl;
            prod.reducefds = {'faces'};
            prod = prod.setup();

            delta = ones(faceNodeTbl.num, 1);

            cellNodeVecTbl1 = prod.tbl3;
            cellNodeVec1    = prod.eval(delta, cellFaceVec);

            A = [[4 1 1 1; ...
                  5 2 2 1; ...
                  6 1 2 1; ...
                  6 2 1 1]];

            clear linVec12Tbl
            linVec12Tbl.lin  = A(:, 1);
            linVec12Tbl.vec1 = A(:, 2);
            linVec12Tbl.vec2 = A(:, 3);
            linVec12Tbl = IndexArray(linVec12Tbl);

            linVec12 = A(:, 4);

            prod = TensorProd();
            prod.tbl1 = linVec12Tbl;
            prod.tbl2 = cellNodeVecTbl1;
            prod.replacefds1 = {{'vec1', 'vec'}};
            prod.replacefds2 = {{'vec', 'vec2'}};
            prod.reducefds = {'vec2'};
            prod = prod.setup();

            cellNodeVecLinTbl1 = prod.tbl3;
            cellNodeVecLin1    = prod.eval(0.5*linVec12, cellNodeVec1);

            [cellNodeVecLinGtypeGindTbl1, indstruct] = crossIndexArray(cellNodeVecLinTbl1, nodeVecGtypeGindTbl, {'nodes', 'vec'});
            cellNodeVecLinGtypeGind1 = cellNodeVecLin1(indstruct{1}.inds);

            prod = TensorProd();
            prod.tbl1 = cellNodeVecLinGtypeGindTbl1;
            prod.tbl2 = cellLinTbl;
            prod.tbl3 = cellGtypeGindTbl;
            prod.reducefds = {'lin'};
            prod.mergefds = {'cells'};
            prod = prod.setup();

            tens = SparseTensor();
            tens = tens.setFromTensorProd(cellNodeVecLinGtypeGind1, prod);

            R1 = tens.getMatrix();

            prod = TensorProd();
            prod.tbl1 = cellNodeVecLinGtypeGindTbl1;
            prod.tbl2 = nk_cellLinTbl;
            prod.tbl3 = cellGtypeGindTbl;
            prod.reducefds = {'lin'};
            prod.mergefds = {'cells'};
            prod = prod.setup();

            tens = SparseTensor();
            tens = tens.setFromTensorProd(cellNodeVecLinGtypeGind1, prod);

            nk_R1 = tens.getMatrix();

            prod = TensorProd();
            prod.tbl1 = linVec12Tbl;
            prod.tbl2 = cellFaceVecTbl;
            prod.replacefds1 = {{'vec1', 'vec'}};
            prod.replacefds2 = {{'vec', 'vec2'}};
            prod.reducefds = {'vec2'};
            prod = prod.setup();

            cellFaceVecLinTbl1 = prod.tbl3;
            cellFaceVecLin1    = prod.eval(linVec12, cellFaceVec);

            prod = TensorProd();
            prod.tbl1 = faceVecTbl;
            prod.tbl2 = cellFaceVecLinTbl1;
            prod.reducefds = {'vec'};
            prod.mergefds = {'faces'};
            prod = prod.setup();

            cellFaceLinTbl1 = prod.tbl3;
            cellFaceLin1    = prod.eval(faceVec, cellFaceVecLin1);

            prod = TensorProd();
            prod.tbl1 = faceTbl;
            prod.tbl2 = cellFaceLinTbl1;
            prod.tbl3 = cellFaceLinTbl1;
            prod.mergefds = {'faces'};
            prod = prod.setup();

            cellFaceLin1 = prod.eval(1./G.faces.areas, cellFaceLin1);

            [cellFaceLinGtypeGindTbl1, indstruct] = crossIndexArray(cellFaceLinTbl1, faceGtypeGindTbl, {'faces'});
            cellFaceLinGtypeGind1 = cellFaceLin1(indstruct{1}.inds);

            prod = TensorProd();
            prod.tbl1 = cellFaceLinGtypeGindTbl1;
            prod.tbl2 = cellLinTbl;
            prod.tbl3 = cellGtypeGindTbl;
            prod.reducefds = {'lin'};
            prod.mergefds = {'cells'};
            prod = prod.setup();

            tens = SparseTensor();
            tens = tens.setFromTensorProd(cellFaceLinGtypeGind1, prod);

            R2 = tens.getMatrix();

            prod = TensorProd();
            prod.tbl1 = cellFaceLinGtypeGindTbl1;
            prod.tbl2 = nk_cellLinTbl;
            prod.tbl3 = cellGtypeGindTbl;
            prod.reducefds = {'lin'};
            prod.mergefds = {'cells'};
            prod = prod.setup();

            tens = SparseTensor();
            tens = tens.setFromTensorProd(cellFaceLinGtypeGind1, prod);

            nk_R2 = tens.getMatrix();

            R = R1 + R2;

            nk_R = nk_R1 + nk_R2;

            %% finishes assembly of M (see eq. 8.36 and 8.37 in reference)

            % NOTE : We should use block inversion

            invNtR = inv(nk_N'*nk_R);
            invNNt = inv(N'*N);

            M1 = (nk_R)*invNtR*(nk_R');

            % add regularisation
            gamma = 1; % NOTE : this should be replaced by proper expression (see eq. 8.37 in reference)
            M2    = gamma*(diag(ones(size(N, 1), 1)) - N*invNNt*(N'));

            M = M1 + M2;

            map = TensorMap();
            map.fromTbl  = gtypeGindTbl;
            map.toTbl    = cellGtypeGindTbl;
            map.mergefds = {'gtype', 'gind'};
            map = map.setup();

            P = SparseTensor();
            P = P.setFromTensorMap(map);
            P = P.getMatrix();

            M = P'*M*P;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Setup of flux operator %
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            prod = TensorProd();
            prod.tbl1 = faceNodeTbl;
            prod.tbl2 = faceVecTbl;
            prod.tbl3 = faceNodeVecTbl;
            prod.mergefds = {'faces'};
            prod = prod.setup();

            delta = ones(faceNodeTbl.num, 1);
            fluxFaceNodeVec = 0.5*prod.eval(delta, faceVec); % faceVec corresponds to the normals

            [faceNodeVecGindGtypeTbl, indstruct] = crossIndexArray(faceNodeVecTbl, nodeVecGtypeGindTbl, {'nodes', 'vec'});

            fluxFaceNodeVecGindGtype = fluxFaceNodeVec(indstruct{1}.inds);

            prod = TensorProd();
            prod.tbl1 = faceNodeVecGindGtypeTbl;
            prod.tbl2 = gtypeGindTbl;
            prod.tbl3 = faceTbl;
            prod.reducefds = {'gind', 'gtype'};
            prod = prod.setup();

            F1 = SparseTensor();
            F1 = F1.setFromTensorProd(fluxFaceNodeVecGindGtype, prod);
            F1 = F1.getMatrix();

            % note that faceGtypeGindTbl by construction has the same indexing as faceTbl
            prod = TensorProd();
            prod.tbl1 = faceGtypeGindTbl;
            prod.tbl2 = gtypeGindTbl;
            prod.tbl3 = faceTbl;
            prod.reducefds = {'gind', 'gtype'};
            prod = prod.setup();

            F2 = SparseTensor();
            F2 = F2.setFromTensorProd(G.faces.areas, prod);
            F2 = F2.getMatrix();

            Flux = F1 + F2;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Setup Divergence operator %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            prod = TensorProd();
            prod.tbl1 = cellFaceTbl;
            prod.tbl2 = faceTbl;
            prod.tbl3 = cellTbl;
            prod.reducefds = {'faces'};
            prod = prod.setup();

            C = SparseTensor();
            C = C.setFromTensorProd(sgn, prod);
            C = C.getMatrix();

            Div = C*Flux;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Setup boundary operators %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Boundary operators for Dirichlet

            dirichletNodeTbl.nodes       = solver.dirichletBc.nodes;
            dirichletNodeTbl             = IndexArray(dirichletNodeTbl);
            dirichletNodeVecTbl          = crossIndexArray(dirichletNodeTbl, vecTbl, {});
            dirichletNodeVecGindGtypeTbl = crossIndexArray(dirichletNodeVecTbl, nodeVecGtypeGindTbl, {'nodes', 'vec'});

            dirichletFaceTbl.faces    = solver.dirichletBc.faces;
            dirichletFaceTbl          = IndexArray(dirichletFaceTbl);
            dirichletFaceGindGtypeTbl = crossIndexArray(dirichletFaceTbl, faceGtypeGindTbl, {'faces'});

            dirichletGindGypeTbl1 = projIndexArray(dirichletNodeVecGindGtypeTbl, {'gtype', 'gind'});
            dirichletGindGypeTbl2 = projIndexArray(dirichletFaceGindGtypeTbl, {'gtype', 'gind'});

            % index for dirichlet degrees of freedom
            dirichletGindGypeTbl = concatIndexArray(dirichletGindGypeTbl1, dirichletGindGypeTbl2, {});

            map = TensorMap();
            map.fromTbl  = gtypeGindTbl;
            map.toTbl    = dirichletNodeVecGindGtypeTbl;
            map.mergefds = {'gind', 'gtype'};
            map = map.setup();

            gPn = SparseTensor();
            gPn = gPn.setFromTensorMap(map);
            gPn = gPn.getMatrix();

            map = TensorMap();
            map.fromTbl  = faceTbl;
            map.toTbl    = dirichletFaceTbl;
            map.mergefds = {'faces'};
            map = map.setup();

            gPf = SparseTensor();
            gPf = gPf.setFromTensorMap(map);
            gPf = gPf.getMatrix();

            gPf = gPf*Flux;

            map = TensorMap();
            map.fromTbl  = dirichletGindGypeTbl;
            map.toTbl    = dirichletNodeVecGindGtypeTbl;
            map.mergefds = {'gind', 'gtype'};
            map = map.setup();

            dPn = SparseTensor();
            dPn = dPn.setFromTensorMap(map);
            dPn = dPn.getMatrix();

            map = TensorMap();
            map.fromTbl  = dirichletGindGypeTbl;
            map.toTbl    = dirichletFaceGindGtypeTbl;
            map.mergefds = {'gind', 'gtype'};
            map = map.setup();

            dPf = SparseTensor();
            dPf = dPf.setFromTensorMap(map);
            dPf = dPf.getMatrix();

            Diri = dPn'*gPn + dPf'*gPf;

            dirOp = struct('gPn',gPn, ...
                           'gPf',gPf, ...
                           'dPn',dPn, ...
                           'dPf',dPf);

            helpers.dirOp = dirOp;

            solver.helpers = helpers;
            
            %% Boundary operators for Neumann

            neumannFaceTbl.faces = solver.neumannFaces;
            neumannFaceTbl = IndexArray(neumannFaceTbl);

            % dof space for neumann bc
            neumannFaceVecTbl = crossIndexArray(neumannFaceTbl, vecTbl, {});

            neumannNodeFaceVecTbl = crossIndexArray(neumannFaceVecTbl, faceNodeTbl, {'faces'});
            neumannNodeFaceVecGtypeGindTbl = crossIndexArray(neumannNodeFaceVecTbl, nodeVecGtypeGindTbl, {'nodes', 'vec'});

            prod = TensorProd();
            prod.tbl1 = neumannNodeFaceVecGtypeGindTbl;
            prod.tbl2 = neumannFaceVecTbl;
            prod.tbl3 = gtypeGindTbl;
            prod.reducefds = {'faces', 'vec'};
            prod = prod.setup();

            Neum = SparseTensor();
            Neum = Neum.setFromTensorProd(0.5*ones(neumannNodeFaceVecGtypeGindTbl.num, 1), prod);
            Neum = Neum.getMatrix();

            operators = struct('M'   , M   , ...
                               'Div' , Div , ...
                               'Flux', Flux, ...
                               'Diri', Diri, ...
                               'Neum', Neum);

            solver.operators = operators;

            n2 = size(Div, 1);
            n3 = size(Diri, 1);
            Z22 = sparse(n2, n2);
            Z23 = sparse(n2, n3);
            Z32 = sparse(n3, n2);
            Z33 = sparse(n3, n3);

            A = [[M   , Div', Diri']; ...
                 [Div , Z22 , Z23]  ; ...
                 [Diri, Z32 , Z33]];

            solver.operators.A = A;

            tbls = struct('cellNodeTbl'        , cellNodeTbl        , ...
                          'vecTbl'             , vecTbl             , ...
                          'nodeTbl'            , nodeTbl            , ...
                          'nodeVecTbl'         , nodeVecTbl         , ...
                          'cellTbl'            , cellTbl            , ...
                          'faceTbl'            , faceTbl            , ...
                          'cellVecTbl'         , cellVecTbl         , ...
                          'faceVecTbl'         , faceVecTbl         , ...
                          'cellFaceTbl'        , cellFaceTbl        , ...
                          'cellFaceVecTbl'     , cellFaceVecTbl     , ...
                          'faceNodeVecTbl'     , faceNodeVecTbl     , ...
                          'faceNodeTbl'        , faceNodeTbl        , ...
                          'extFaceTbl'         , extFaceTbl         , ...
                          'gtypeGindTbl'       , gtypeGindTbl       , ...
                          'cellGtypeGindTbl'   , cellGtypeGindTbl   , ...
                          'nodeVecGtypeGindTbl', nodeVecGtypeGindTbl, ...
                          'faceGtypeGindTbl'   , faceGtypeGindTbl);

            solver.tbls = tbls;
            
        end

        function u = getNodalVelocity(solver, x)

            dim = solver.dim;

            tbls                = solver.tbls;
            gtypeGindTbl        = tbls.gtypeGindTbl;
            nodeVecGtypeGindTbl = tbls.nodeVecGtypeGindTbl;

            x = x(1 : gtypeGindTbl.num);
            
            % note that nodeVecGtypeGindTbl has same indexing as nodeVecTbl;
            map = TensorMap();
            map.fromTbl = gtypeGindTbl;
            map.toTbl   = nodeVecGtypeGindTbl;
            map.mergefds = {'gtype', 'gind'};
            map = map.setup();

            u = map.eval(x);

            u = reshape(u, dim, [])';
            
        end

        function p = getPressure(solver, x)

            
            tbls         = solver.tbls;
            gtypeGindTbl = tbls.gtypeGindTbl;
            cellTbl      = tbls.cellTbl;
            
            nu = gtypeGindTbl.num;
            np = cellTbl.num;

            p = x(nu + (1 : np));
            
        end

        
    end

    
end


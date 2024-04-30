dim = 2;

close all
clear all

gridtype = 'radial';

viscosity = 1;

switch gridtype

  case 'debugging'

    lx = 3;
    ly = 1;
    G = cartGrid([5, 5], [lx, ly]);
    G = computeGeometry(G);

    extfaces = find(any(G.faces.neighbors == 0, 2));

    % Setup some Dirichlet boundary faces
    
    dirichletNodes1 = find(G.nodes.coords(:, 1) == 0);
    dirichletNodes2 = find(G.nodes.coords(:, 1) == lx);

    dirichletBc.nodes = vertcat(dirichletNodes1, dirichletNodes2);

    dirichletFaces1 = find(G.faces.centroids(:, 1) == 0);
    dirichletFaces2 = find(G.faces.centroids(:, 1) == lx);

    dirichletBc.faces = vertcat(dirichletFaces1, dirichletFaces2);

    % Setup some Neumann boundary faces
    ind = ismember(extfaces, dirichletBc.faces);
    neumannFaces = extfaces(~ind);

    
    op = solver.helpers.dirOp;

    v = reshape(repmat([1, 1], numel(solver.dirichletBc.nodes), 1)', [], 1); % in nodevec format
    v = op.dPn'*v;
    
  case 'pebi'

    mrstModule add upr
    
    l = {[0.1,0.42; 0.4,.55; 0.7,0.65], ...
         [0.8,0.13; 0.6,0.4; 0.55,0.6],...
         [0.42,1.08; 0.45,0.9; 0.5,0.8; 0.58,0.6]};
    
    gS = [1/30,1/30];
    
    %% Create grid
    % We use the routine compositePebiGrid2D to create the grid. We use the
    % preset options. 
    G = compositePebiGrid2D(gS, [1,1.15], 'faceConstraints',l, 'useMrstPebi', true);
    G = computeGeometry(G);
    
    extfaces = find(any(G.faces.neighbors == 0, 2));

    % Setup some Dirichlet boundary faces
    
    dirichletNodes1 = find(G.nodes.coords(:, 1) == 0);
    dirichletNodes2 = find(G.nodes.coords(:, 1) == 1);

    dirichletBc.nodes = vertcat(dirichletNodes1, dirichletNodes2);

    dirichletFaces1 = find(G.faces.centroids(:, 1) == 0);
    dirichletFaces2 = find(G.faces.centroids(:, 1) == 1);

    dirichletBc.faces = vertcat(dirichletFaces1, dirichletFaces2);

    % Setup some Neumann boundary faces
    ind = ismember(extfaces, dirichletBc.faces);
    neumannFaces = extfaces(~ind);

    op = solver.helpers.dirOp;

    v = reshape(repmat([1, 1], numel(solver.dirichletBc.nodes), 1)', [], 1); % in nodevec format
    v = op.dPn'*v;

  case 'radial'

    mrstModule add nwm
    
    %% Make a layered radial-Cartesian hybrid grid with inclination
    % Build the Cartesian grid
    G = cartGrid([30, 30], [4, 4]);
    G = computeGeometry(G);

    tol = 1e-6;
    % Define the well region by logical indices
    cI = find(G.cells.centroids(:, 1)   - 1 >= - tol   ...
              & G.cells.centroids(:, 1) - 3 <= tol ...
              & G.cells.centroids(:, 2) - 1 >= - tol   ...
              & G.cells.centroids(:, 2) - 3 <= tol ...
             );


    % Place the well at the region center
    pW  = [2, 2];

    % Define radial parameters
    [nR, rW, rM] = deal(4, 0.5, 1);

    % Get the hybrid grid
    G = radCartHybridGrid(G, cI, rW, rM, nR, pW);

    tbls = setupTables(G, 'includetbls', {'extfacetbl', 'vectbl'});
    extfaceTbl  = tbls.extfacetbl;
    faceNodeTbl = tbls.facenodetbl;
    vecTbl      = tbls.vectbl;
    
    extfaces = extfaceTbl.get('faces');
    c = G.faces.centroids(extfaces, :);

    cr = bsxfun(@plus, c, -pW);

    clear radfaceTbl
    radfaceTbl.faces = extfaces(sum(cr.^2, 2) <= 1.5);
    radfaceTbl = IndexArray(radfaceTbl);
    radfaceIndTbl = radfaceTbl.addInd('ind', ones(radfaceTbl.num, 1));
    
    clear hdirfaceTbl
    hdirfaceTbl.faces = extfaces( c(:, 1) <= 0 + tol |   c(:, 1) >= 4 - tol  );
    hdirfaceTbl = IndexArray(hdirfaceTbl);
    hdirfaceIndTbl = hdirfaceTbl.addInd('ind', 2*ones(hdirfaceTbl.num, 1));

    dirfaceIndTbl = concatIndexArray(radfaceIndTbl, hdirfaceIndTbl);
    
    dirfaceNodeIndTbl = crossIndexArray(dirfaceIndTbl, faceNodeTbl, {'faces'});
    dirnodeIndTbl     = projIndexArray(dirfaceNodeIndTbl, {'nodes', 'ind'});
    
    clear dirichletBc
    dirichletBc.faces = dirfaceIndTbl.get('faces');
    dirichletBc.nodes = dirnodeIndTbl.get('nodes');
    
    doplot = true;
    if doplot
        plotGrid(G);
        plotFaces(G, radfaceTbl.get('faces'), 'edgecolor', 'red')
        plotFaces(G, hdirfaceTbl.get('faces'), 'edgecolor', 'green')
    end
    
    % Setup some Neumann boundary faces
    ind = ismember(extfaces, dirichletBc.faces);
    neumannFaces = extfaces(~ind);

    solver = StokesSolver(G           , ...
                          dirichletBc , ...
                          neumannFaces, ...
                          viscosity);
    
    dirichletNodeVecGindGtypeTbl = solver.tbls.dirichletNodeVecGindGtypeTbl;
    
    indVecTbl1 = vecTbl;
    indVecTbl1 = indVecTbl1.addInd('ind', 2);
    
    prod = TensorProd();
    prod.tbl1 = indVecTbl1;
    prod.tbl2 = dirnodeIndTbl;
    prod.reducefds = {'ind'};
    prod = prod.setup();

    dirnodeVecTbl1 = prod.tbl3;
    v = prod.eval([1; 0], ones(dirnodeIndTbl.num, 1));

    map = TensorMap();
    map.fromTbl = dirnodeVecTbl1;
    map.toTbl = dirichletNodeVecGindGtypeTbl;
    map.mergefds = {'nodes', 'vec'};
    map = map.setup();

    v = map.eval(v);

    op = solver.helpers.dirOp;
    v = op.dPn'*v;
    
    A    = solver.operators.A;
    Diri = solver.operators.Diri;

    b = zeros(size(A, 1), 1);
    b(end - size(Diri, 1) + 1 : end) = v;

    x = A\b;

  otherwise
    
    error('case not recognized');
    
end



u = solver.getNodalVelocity(x);
p = solver.getPressure(x);


%% Example for the stokes solver.
% It requires the following extra modules to setup the irregular grid
% nwm : https://bitbucket.org/LinZhao9/nwm.git
% upr : https://github.com/rbe051/UPR.git
% Download those and put them in your path using : addpath(genpath('.'))
%
% The option debugging works without those (not that we get a constant pressure in this example)

dim = 2;

close all
clear all

gridtype = 'debugging';

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

    solver = StokesSolver(G           , ...
                          dirichletBc , ...
                          neumannFaces, ...
                          viscosity);

    op = solver.helpers.dirOp;

    v = reshape(repmat([1, 1], numel(solver.dirichletBc.nodes), 1)', [], 1); % in nodevec format
    v = op.dPn'*v;
    
  case 'pebi'

    mrstModule add upr
    
    l = {[0.1, 0.42; 0.4,.55; 0.7, 0.65]             , ...
         [0.8, 0.13; 0.6, 0.4; 0.55, 0.6]            , ...
         [0.42, 0.91; 0.45, 0.9; 0.5, 0.8; 0.58, 0.6], ...
         [0.1, 0.1; 0.8, 0.5; 0.9, 0.9]
        };
    
    gS = [1/30,1/30];
    
    %% Create grid
    % We use the routine compositePebiGrid2D to create the grid. We use the
    % preset options. 
    G = compositePebiGrid2D(gS, [1, 1], 'faceConstraints',l, 'useMrstPebi', true);
    G = computeGeometry(G);

    plotGrid(G);
    
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

    solver = StokesSolver(G           , ...
                          dirichletBc , ...
                          neumannFaces, ...
                          viscosity);

    op = solver.helpers.dirOp;

    v = reshape(repmat([1, 1], numel(solver.dirichletBc.nodes), 1)', [], 1); % in nodevec format
    v = op.dPn'*v;

  case 'radial'

    mrstModule add nwm

    background_grid = 'pebi';
    
    lx = 1;
    ly = 1;

    switch background_grid

      case 'cart'

        nx = 30;
        ny = 30;

        G = cartGrid([nx, nx], [lx, ly]);
        G = computeGeometry(G);

        xldir = 0;
        
      case 'pebi'

        mrstModule add upr

        nx = 30;
        ny = 30;

        G1 = cartGrid([nx, nx], [lx, ly]);
        G1 = computeGeometry(G1);
        
        wl = {[0.2,0.8;0.8,0.2]};
        G2  = pebiGrid2D(1/10,[1,1],'cellConstraints',wl);
        G2.nodes.coords = bsxfun(@plus, G2.nodes.coords, -[1, 0]);

        G = glue2DGrid(G1, G2);
        G.cells = rmfield(G.cells, 'indexMap');

        G = computeGeometry(G);
        
        xldir = -1;
        
      otherwise
        
        error('background_grid not recognized');
        
    end

    nR = 10;
    rW = 0.2;
    rM = 0.3;
    rE = 0.31;
    
    pW = [lx/2, ly/2];
    
    tol = 1e-6;
    % Define the well region by logical indices
    cI = find(G.cells.centroids(:, 1)   - (pW(1) - rE) >= - tol   ...
              & G.cells.centroids(:, 1) - (pW(1) + rE) <= tol ...
              & G.cells.centroids(:, 2) - (pW(2) - rE) >= - tol   ...
              & G.cells.centroids(:, 2) - (pW(2) + rE) <= tol ...
             );
    % plotGrid(G)
    % plotGrid(G, cI, 'facecolor', 'red')

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
    radfaceTbl.faces = extfaces(sum(cr.^2, 2) <= rW + tol);
    radfaceTbl = IndexArray(radfaceTbl);
    radfaceIndTbl = radfaceTbl.addInd('ind', ones(radfaceTbl.num, 1));
    
    clear hdirfaceTbl
    hdirfaceTbl.faces = extfaces( c(:, 1) <= xldir + tol |   c(:, 1) >= lx - tol  );
    % hdirfaceTbl.faces = extfaces( c(:, 1) <= xldir + tol);
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
        plotFaces(G, radfaceTbl.get('faces'), 'edgecolor', 'red', 'linewidth', 3)
        plotFaces(G, hdirfaceTbl.get('faces'), 'edgecolor', 'green', 'linewidth', 3)
        title('boundary conditions')
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
    

  otherwise
    
    error('case not recognized');
    
end

%%

A    = solver.operators.A;
Diri = solver.operators.Diri;

b = zeros(size(A, 1), 1);
b(end - size(Diri, 1) + 1 : end) = v;

x = A\b;

u  = solver.getNodalVelocity(x);
p  = solver.getPressure(x);

xn = G.nodes.coords;

figure
hold on
plotGrid(G, 'facecolor', 'none');
quiver(xn(:, 1), xn(:, 2), u(:, 1), u(:, 2), 'linewidth', 2);
title('velocity field')

figure
hold on
plotCellData(G, p);
title('pressure')
colorbar

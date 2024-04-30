dim = 2;
close all
clear all
gridtype = 'pebi';

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

    viscosity = 1;
    
  case 'complex'

    load('G.mat');

  case 'pebi'

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

    viscosity = 1;

  case 'radial'

    
    %% Make a layered radial-Cartesian hybrid grid with inclination
    % Build the Cartesian grid
    GC = cartGrid([15, 15], [300, 300]);
    GC = computeGeometry(GC);

    % Define the well region by logical indices
    ij = gridLogicalIndices(GC);
    idxI = ij{1} >= 6 & ij{1} <= 10 & ij{2} >= 6 & ij{2} <= 10;
    cI = find(idxI); 

    % Place the well at the region center
    pW  = [150, 150];

    % Define radial parameters
    [nR, rW, rM] = deal(10, 0.2, 35);

    % Get the hybrid grid
    G = radCartHybridGrid(GC, cI, rW, rM, nR, pW);

  otherwise
    error('case not recognized');
end

return
solver = StokesSolver(G           , ...
                      dirichletBc , ...
                      neumannFaces, ...
                      viscosity);

op = solver.helpers.dirOp;

v = reshape(repmat([1, 1], numel(solver.dirichletBc.nodes), 1)', [], 1); % in nodevec format
v = op.dPn'*v;

A    = solver.operators.A;
Diri = solver.operators.Diri;

b = zeros(size(A, 1), 1);
b(end - size(Diri, 1) + 1 : end) = v;

x = A\b;

u = solver.getNodalVelocity(x);
p = solver.getPressure(x);


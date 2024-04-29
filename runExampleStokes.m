dim = 2;
close all
gridtype = 'debugging';

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
    ind = ismember(extfaces, dirichletFaces);
    neumannFaces = extfaces(~ind);

    viscosity = 1;
    
  case 'complex'
    load('G.mat');
  otherwise
    error('case not recognized');
end


solver = StokesSolver(G           , ...
                      dirichletBc , ...
                      neumannFaces, ...
                      viscosity);

op = solver.helpers.dirOp;

v = reshape(repmat([1, 0], numel(solver.dirichletBc.nodes), 1)', [], 1); % in nodevec format
vn = op.dPn'*v;

v = ones(numel(solver.dirichletBc.faces), 1);
vf = op.dPf'*v;

v = vn + vf;

A    = solver.operators.A;
Diri = solver.operators.Diri;

b = zeros(size(A, 1), 1);
b(1 : size(Diri', 1)) = Diri'*v;

x = A\b;

u = solver.getNodalVelocity(x);



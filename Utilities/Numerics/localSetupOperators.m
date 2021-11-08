function operators = localSetupOperators(G, varargin)
    opts = struct('assembleCellFluxOperator', false);
    opts = merge_options(opts, varargin{:});
    
    nc = G.cells.num;
    cst = ones(nc, 1);
    rock = struct('perm', cst, 'poro', cst);
    operators = setupOperatorsTPFA(G, rock);
    operators.allDiv = getAllDiv(G);
    %operators.harmFace = getFaceHarmMean(G);
    operators.harmFace =@(cellvalue) getTrans(G).*(1./operators.faceAvg(1./cellvalue));
    operators.harmFaceBC = @(cvalue, faces) getFaceHarmBC(G, cvalue, faces);
    
    if opts.assembleCellFluxOperator
        operators.cellFluxOp = getCellFluxOperator(G);
    end
end

function [T, cells] = getFaceHarmBC(G, cvalue, faces)
    cells = sum(G.faces.neighbors(faces, :), 2);
    cn = sqrt(sum((G.faces.centroids(faces, :) - G.cells.centroids(cells, :)).^2, 2));
    t = G.faces.areas(faces)./cn;
    T = t.*cvalue(cells);
end

function hm = getFaceHarmMean(G)
    internal = all(G.faces.neighbors>0, 2);
    N = G.faces.neighbors(internal, :);
    ni = sum(internal);
    cd = sqrt(sum((G.cells.centroids(N(:, 1), :) - G.cells.centroids(N(:, 2), :)).^2, 2)); % NB
    t = G.faces.areas(internal)./cd;
    A = sparse([[1:ni]'; [1:ni]'], N, 1, ni, G.cells.num);
    hm = @(cellvalue) 2.*t./(A*(1./cellvalue));
end

function t = getTrans(G)
    internal = all(G.faces.neighbors>0, 2);
    N = G.faces.neighbors(internal, :);
    cd = sqrt(sum((G.cells.centroids(N(:, 1), :) - G.cells.centroids(N(:, 2), :)).^2, 2)); % NB
    t = G.faces.areas(internal)./cd;
end


function tp = getTwoPointOperator(G)
% Mappings from cells to its faces
    cells = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    faces = G.cells.faces(:, 1);
    % Vector from cell to face centroid
    C = G.faces.centroids(faces, :) - G.cells.centroids(cells, :);
    % Oriented normals
    sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1;
    N = bsxfun(@times, sgn, G.faces.normals(faces, :));
    % Make function
    cn = sum(C.*N, 2)./sum(C.*C, 2);
    tp = @(lambda) cn.*lambda(cells);
end

function ha = getHarmonicAvgOperator(G)
% Harmonic averaging operator
    faces = G.cells.faces(:, 1);
    M = sparse(faces, 1:numel(faces), 1, G.faces.num, numel(faces));
    ha = @(T) 1./(M*(1./T));
end

function allDiv = getAllDiv(G)
    nc = G.cells.num;
    nf = G.faces.num;
    Nall = G.faces.neighbors;
    internal = all(Nall>0, 2);
    ifn = find(internal);
    efn = find(~internal);
    inf = numel(ifn);
    N = Nall(internal, :);
    Nb = Nall(~internal, :);
    signb = 2*(Nb(:, 1)>0) - 1;
    Nb = sum(Nb, 2);
    C = sparse([ifn; ifn], N, ones(inf, 1) * [1,- 1], nf, nc);
    C = C + sparse(efn, Nb, signb, nf, nc);
    allDiv = @(x) (C'*x)./G.cells.volumes;
end

dim = 2;

G = cartGrid([1, 1], [1, 3]);
G = computeGeometry(G);

tbls = setupTables(G, 'includetbls', {'cellnodetbl', 'vectbl'});

cellNodeTbl = tbls.cellnodetbl;
vecTbl      = tbls.vectbl;
nodeTbl     = tbls.nodetbl;
nodeVecTbl  = tbls.nodevectbl;
cellTbl     = tbls.celltbl;
cellVecTbl  = tbls.cellvectbl;

cellNodeVecTbl = crossIndexArray(cellNodeTbl, vecTbl, {});
cellNodeVecTbl = sortIndexArray(cellNodeVecTbl, {'cells', 'vec', 'nodes'});

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
firstLinVec12 = prod.eval(delta, linVecPol);

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

prod = TensorProd();
prod.tbl1 = cellNodeVecLinTbl;
prod.tbl2 = linTbl;
prod.tbl3 = cellNodeVecTbl;
prod.reducefds = {'lin'};
prod = prod.setup();

tens = SparseTensor();
tens = tens.setFromTensorProd(cellNodeVecLin, prod);

M = tens.getMatrix();

% full(M)



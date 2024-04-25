dim = 2;

G = cartGrid([1, 1], [2, 3]);
G = computeGeometry(G);

tbls = setupTables(G, 'includetbls', {'cellnodetbl', 'vectbl'});

cellNodeTbl = tbls.cellnodetbl;
vecTbl      = tbls.vectbl;
nodeTbl     = tbls.nodetbl;
nodeVecTbl  = tbls.nodevectbl;
cellTbl     = tbls.celltbl;
cellVecTbl  = tbls.cellvectbl;

cellNodeVecTbl = crossIndexArray(cellNodeTbl, vecTbl, {});

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
prod.tbl1 = linVecPoTbl;
prod.tbl2 = vecPolTbl;
prod.tbl3 = linVecTbl;
prod.reducefds = {'pol'};
prod = prod.setup();

delta = ones(vecPolTbl.num, 1);

cellNodePol = prod.eval(linVecPol, delta);

clear constPolTbl
constPolTbl.pol = 1;
constPolTbl = IndexArray(constPolTbl);

map = TensorMap;
map.fromTbl = constPolTbl;
map.toTbl = cellNodePolTbl;
map.mergefds =  {'pol'};
map = map.setup;

cellNodePol = map.eval(1) + cellNodePol;

cellNodeLinVecPolTbl = crossIndexArray(cellNodePolTbl, linVecPolTbl, {'pol'});

map1 = TensorMap();
map1.fromTbl  = cellNodePolTbl;
map1.toTbl    = cellNodeLinVecPolTbl;
map1.mergefds = {'cells', 'nodes', 'pol'};
map1 = map1.setup();

map2 = TensorMap();
map2.fromTbl  = linVecPolTbl;
map2.toTbl    = cellNodeLinVecPolTbl;
map2.mergefds = {'lin', 'vec', 'pol'};
map2 = map2.setup();

cellNodeLinVecPol = map1.eval(cellNodePol).*map2.eval(linVecPol);


gen = CrossIndexArrayGenerator();
gen.tbl1 = cellNodeVecTbl;
gen.tbl2 = linVecPolTbl;
gen.replacefds1 = {{'vec', 'vec1'}};
gen.replacefds2 = {{'vec', 'vec2'}};
gen.mergefds = {};

cellNodeVec1LinVec2PolTbl = gen.eval();

prod = TensorProd();
prod.tbl1 = vecPolTbl;
prod.tbl2 = cellNodeLinVecPolTbl;
prod.tbl3 = cellNodeVec1LinVec2PolTbl;
prod.replacefds1 = {{'vec', 'vec1'}};
prod.replacefds2 = {{'vec', 'vec2'}};
prod.mergefds = {'pol'};
prod = prod.setup();

u = ones(vecPolTbl.num, 1);
N = prod.eval(u, cellNodeLinVecPol);


prod = TensorProd();
prod.tbl1 = cellNodeVec1LinVec2PolTbl;
prod.tbl2 = cellNodeVecTbl;
prod.tbl3 = cellLinVecPolTbl;
prod.replacefds2 = {{'vec', 'vec1'}};
prod.replacefds3 = {{'vec', 'vec2'}};
prod.mer

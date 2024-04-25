dim = 2;

G = cartGrid([2, 1], [1, 3]);
G = computeGeometry(G);

tbls = setupTables(G, 'includetbls', {'cellnodetbl', 'vectbl'});

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
faceNodeTbl    = tbls.facenodetbl;

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

clear inv_linTbl
inv_linTbl.lin = (dim*(dim + 1)/2 + 1 : linTbl.num)';
inv_linTbl = IndexArray(inv_linTbl);

cellLinTbl     = crossIndexArray(cellTbl, linTbl, {});
inv_cellLinTbl = crossIndexArray(cellTbl, inv_linTbl, {});

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
prod.tbl2 = inv_cellLinTbl;
prod.tbl3 = cellGtypeGindTbl;
prod.reducefds = {'lin'};
prod.mergefds  = {'cells'};
prod = prod.setup();

tens = SparseTensor();
tens = tens.setFromTensorProd(cellNodeVecLinGtypeGind, prod);

inv_N = tens.getMatrix();

%% Setup of R

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
prod.tbl2 = inv_cellLinTbl;
prod.tbl3 = cellGtypeGindTbl;
prod.reducefds = {'lin'};
prod.mergefds = {'cells'};
prod = prod.setup();

tens = SparseTensor();
tens = tens.setFromTensorProd(cellNodeVecLinGtypeGind1, prod);

inv_R1 = tens.getMatrix();

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
prod.tbl2 = inv_cellLinTbl;
prod.tbl3 = cellGtypeGindTbl;
prod.reducefds = {'lin'};
prod.mergefds = {'cells'};
prod = prod.setup();

tens = SparseTensor();
tens = tens.setFromTensorProd(cellFaceLinGtypeGind1, prod);

inv_R2 = tens.getMatrix();

R = R1 + R2;

inv_R = inv_R1 + inv_R2;

full(inv_N'*inv_R)








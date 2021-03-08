close all
clear all

sepnz  = 10;
nenz   = 10;
penz   = 10;

elnz = sepnz + nenz + penz;

ccnenx = 2;
ccpenx = 2;
elnx   = 10;

ccneny = 2;
ccpeny = 2;
elny   = 10;

nxs = [ccnenx; elnx; ccpenx];
nys = [ccneny; elny; ccpeny];
nzs = [nenz; sepnz; penz];

xlength = 3e-2*[0.1; 1; 0.1];
ylength = 1e-2*[0.1; 1; 0.1];
zlength = 1e-6*[100; 50; 80];

x = xlength./nxs;
x = rldecode(x, nxs);
x = [0; cumsum(x)];

y = ylength./nys;
y = rldecode(y, nys);
y = [0; cumsum(y)];

z = zlength./nzs;
z = rldecode(z, nzs);
z = [0; cumsum(z)];

G = tensorGrid(x, y, z);

nx = sum(nxs);
ny = sum(nys);
nz = sum(nzs);

dimGlobGrid = [nx; ny; nz];

%% setup elyte
startSubGrid = [ccnenx + 1; ccpeny + 1; 1];
dimSubGrid   = [elnx; elny; elnz];
elytecells   = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

%% setup ne
startSubGrid = [ccnenx + 1; ccpeny + 1; 1];
dimSubGrid   = [elnx; elny; nenz];
necells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

%% setup pe

startSubGrid = [ccnenx + 1; ccpeny + 1; nenz + sepnz + 1];
dimSubGrid   = [elnx; elny; penz];
pecells      = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

%% setup ccne

startSubGrid = [1; ccpeny + 1; 1];
dimSubGrid   = [ccnenx; elny + ccneny; nenz];
ccnecells    = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

%% setup ccpe

startSubGrid = [ccnenx + elnx + 1; 1; nenz + sepnz];
dimSubGrid   = [ccpenx; elny + ccpeny; penz];
ccpecells    = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobGrid);

%% 

cellid = zeros(G.cells.num);
cellid(elytecells) = 1;
cellid(necells)    = 2;
cellid(pecells)    = 3;
cellid(ccnecells)  = 4;
cellid(ccpecells)  = 5;

cells = [elytecells; necells; pecells; ccnecells; ccpecells];

rcells = setdiff((1 : G.cells.num)', cells);

nGlob = G.cells.num;

[G, cellmap, facemap, nodemap] = removeCells(G, rcells);

globtbl.gind = (1 : nGlob)';
globtbl = IndexArray(globtbl);

globloctbl.gind = cellmap;
globloctbl.lind = (1 : G.cells.num)';
globloctbl = IndexArray(globloctbl);

map = TensorMap();
map.fromTbl = globtbl;
map.toTbl = globloctbl;
map.mergefds = {'gind'};
map = map.setup();

cellid = map.eval(cellid);

plotCellData(G, cellid);


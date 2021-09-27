clear all
close all

mrstModule add ad-core multimodel mrst-gui battery mpfa

%% Parameters 

% number of layers
nlayers = 4;

% "radius" at the middle
r0 = 0.5;
% widths of each component ordered as
%  - negative current collector
%  - negative electrode
%  - electrolyte separator 
%  - positive electrode
%  - positive current collector
widths = [1; 2; 3; 2; 1];

% number of cell in radial direction for each component (same ordering as above).
nrs = [3; 3; 3; 3; 3]; 
% number of cells in the angular direction
nas = 20; 


params = struct('nlayers', nlayers, ...
                'r0'     , r0     , ...
                'widths' , widths , ...
                'nrs'    , nrs    , ...
                'nas'    , nas);

%% Grid setup



layerwidth = sum(widths);

w = widths./nrs;
w = rldecode(w, nrs);

w = repmat(w, [nlayers, 1]);
w = [0; cumsum(w)];

h = linspace(0, 2*pi*r0, nas*nlayers + 1);

nperlayer = sum(nrs);

G = tensorGrid(h, w);

m = numel(w) - 1;
n = numel(h) - 1;

% plotGrid(G)

% We roll the domain into a spirale
dotransform = true;
if dotransform
    x = G.nodes.coords(:, 1);
    y = G.nodes.coords(:, 2);

    theta = x./r0;

    G.nodes.coords(:, 1) = (r0 + y + (theta/(2*pi))*layerwidth).*cos(theta);
    G.nodes.coords(:, 2) = (r0 + y + (theta/(2*pi))*layerwidth).*sin(theta);
end

tbls = setupSimpleTables(G);

% We add cartesian indexing for the nodes
nodetbl.nodes = (1 : G.nodes.num)';
nodetbl.indi = repmat((1 : (n + 1))', m + 1, 1);
nodetbl.indj = rldecode((1 : (m + 1))', (n + 1)*ones(m + 1, 1));
nodetbl = IndexArray(nodetbl);

% We add cartesian indexing for the vertical faces (in original cartesian block)
vertfacetbl.faces = (1 : (n + 1)*m)';
vertfacetbl.indi = repmat((1 : (n + 1))', m, 1);
vertfacetbl.indj = rldecode((1 : m)', (n + 1)*ones(m, 1));
vertfacetbl = IndexArray(vertfacetbl);

% Add structure to merge the nodes
node2tbl.indi1 = ones(m - nperlayer + 1, 1);
node2tbl.indj1 = ((nperlayer + 1) : (m + 1))';
node2tbl.indi2 = (n + 1)*ones(m - nperlayer + 1, 1);
node2tbl.indj2 = (1 : (m - nperlayer + 1))';
node2tbl = IndexArray(node2tbl);

gen = CrossIndexArrayGenerator();
gen.tbl1 = nodetbl;
gen.tbl2 = node2tbl;
gen.replacefds1 = {{'indi', 'indi1'}, {'indj', 'indj1'}, {'nodes', 'nodes1'}};
gen.mergefds = {'indi1', 'indj1'};

node2tbl = gen.eval();

gen = CrossIndexArrayGenerator();
gen.tbl1 = nodetbl;
gen.tbl2 = node2tbl;
gen.replacefds1 = {{'indi', 'indi2'}, {'indj', 'indj2'}, {'nodes', 'nodes2'}};
gen.mergefds = {'indi2', 'indj2'};

node2tbl = gen.eval();

node2tbl = sortIndexArray(node2tbl, {'nodes1', 'nodes2'});

% Add structure to merge the faces
face2tbl.indi1 = ones(m - nperlayer, 1);
face2tbl.indj1 = nperlayer + (1 : (m - nperlayer))';
face2tbl.indi2 = (n + 1)*ones(m - nperlayer, 1);
face2tbl.indj2 = (1 : (m - nperlayer))';
face2tbl = IndexArray(face2tbl);

gen = CrossIndexArrayGenerator();
gen.tbl1 = vertfacetbl;
gen.tbl2 = face2tbl;
gen.replacefds1 = {{'indi', 'indi1'}, {'indj', 'indj1'}, {'faces', 'faces1'}};
gen.mergefds = {'indi1', 'indj1'};

face2tbl = gen.eval();

gen = CrossIndexArrayGenerator();
gen.tbl1 = vertfacetbl;
gen.tbl2 = face2tbl;
gen.replacefds1 = {{'indi', 'indi2'}, {'indj', 'indj2'}, {'faces', 'faces2'}};
gen.mergefds = {'indi2', 'indj2'};

face2tbl = gen.eval();


%% We setup the new indexing for the nodes

nodetoremove = node2tbl.get('nodes2');
newnodes = (1 : G.nodes.num)';
newnodes(nodetoremove) = [];

newnodetbl.newnodes = (1 : numel(newnodes))';
newnodetbl.nodes = newnodes;
newnodetbl = IndexArray(newnodetbl);

gen = CrossIndexArrayGenerator();
gen.tbl1 = node2tbl;
gen.tbl2 = newnodetbl;
gen.replacefds2 = {{'nodes', 'nodes1'}};
gen.mergefds = {'nodes1'};

node2tbl = gen.eval();

newnodes = [newnodetbl.get('newnodes'); node2tbl.get('newnodes')];
nodes = [newnodetbl.get('nodes'); node2tbl.get('nodes2')];

clear newnodetbl;
newnodetbl.newnodes = newnodes;
newnodetbl.nodes = nodes;
newnodetbl = IndexArray(newnodetbl);

%% We setup the new indexing for the faces

facetoremove = face2tbl.get('faces2');
newfaces = (1 : G.faces.num)';
newfaces(facetoremove) = [];

clear facetbl
newfacetbl.newfaces = (1 : numel(newfaces))';
newfacetbl.faces = newfaces;
newfacetbl = IndexArray(newfacetbl);

gen = CrossIndexArrayGenerator();
gen.tbl1 = face2tbl;
gen.tbl2 = newfacetbl;
gen.replacefds2 = {{'faces', 'faces1'}};
gen.mergefds = {'faces1'};

face2tbl = gen.eval();

newfaces = [newfacetbl.get('newfaces'); face2tbl.get('newfaces')];
faces = [newfacetbl.get('faces'); face2tbl.get('faces2')];

allnewfacetbl.newfaces = newfaces;
allnewfacetbl.faces = faces;
allnewfacetbl = IndexArray(allnewfacetbl);

%% we maps from old to new

cellfacetbl = tbls.cellfacetbl;
% we store the order previous to mapping. Here we just assumed that the grid is cartesian for simplicity
cellfacetbl = cellfacetbl.addInd('order', repmat((1 : 4)', G.cells.num, 1));

cellfacetbl = crossIndexArray(cellfacetbl, allnewfacetbl, {'faces'});
cellfacetbl = sortIndexArray(cellfacetbl, {'cells', 'order' 'newfaces'});
cellfacetbl = replacefield(cellfacetbl, {{'newfaces', 'faces'}});

facenodetbl = tbls.facenodetbl;
% facenodetbl = facenodetbl.addInd('order', repmat((1 : 2)', G.faces.num, 1));
facenodetbl = crossIndexArray(facenodetbl, newfacetbl, {'faces'});
facenodetbl = crossIndexArray(facenodetbl, newnodetbl, {'nodes'});

facenodetbl = sortIndexArray(facenodetbl, {'newfaces',  'newnodes'});
facenodetbl = replacefield(facenodetbl, {{'newfaces', 'faces'}, {'newnodes', 'nodes'}});

clear nodes
nodes.coords = G.nodes.coords(newnodetbl.get('nodes'), :);
nodes.num = size(nodes.coords, 1);

clear faces
[~, ind] = rlencode(facenodetbl.get('faces'));
faces.nodePos = [1; 1 + cumsum(ind)];
faces.nodes = facenodetbl.get('nodes');
faces.num = newfacetbl.num;

clear cells
[~, ind] = rlencode(cellfacetbl.get('cells'));
cells.facePos = [1; 1 + cumsum(ind)];
cells.faces = cellfacetbl.get('faces');
cells.num = G.cells.num;


newG.cells = cells;
newG.faces = faces;
newG.nodes = nodes;
newG.griddim = 2;

newG = computeGeometry(newG);

% plotGrid(newG)

comptag = rldecode((1 : 5)', nrs);
comptag = repmat(comptag, [nlayers, 1]);

comptagtbl.tag = comptag;
comptagtbl.indj = (1 : (sum(nrs)*nlayers))';
comptagtbl = IndexArray(comptagtbl);

celltbl.cells = (1 : G.cells.num)';
celltbl.indi = repmat((1 : nas*nlayers)', [sum(nrs)*nlayers, 1]);
celltbl.indj = rldecode((1 : sum(nrs)*nlayers)', nas*nlayers*ones(sum(nrs)*nlayers, 1));
celltbl = IndexArray(celltbl);

celltagtbl = crossIndexArray(celltbl, comptagtbl, {'indj'});
celltagtbl = sortIndexArray(celltagtbl, {'cells', 'tag'});

tag = celltagtbl.get('tag');

plotGrid(newG);
plotCellData(newG, tag);


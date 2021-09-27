clear all
close all

mrstModule add ad-core multimodel mrst-gui battery mpfa

nelyte = 3;
npe    = 3;
nne    = 3;
nccpe  = 3;
nccne  = 3;

nlayers = 2;
r0 = 0.5;
n = 20; 

ns = [nccne; nne; nelyte; nccpe; npe];
compinds = rldecode((1 : 5)', ns);
widths = [1; 2; 3; 2; 1];

layerwidth = sum(widths);

w = widths./ns;
w = rldecode(w, ns);

w = repmat(w, [nlayers, 1]);
w = [0; cumsum(w)];

h = linspace(0, 2*pi*r0, n);

nperlayer = sum(ns);

G = tensorGrid(h, w);

n = numel(h);
m = numel(w);

plotGrid(G)

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

plotGrid(newG)



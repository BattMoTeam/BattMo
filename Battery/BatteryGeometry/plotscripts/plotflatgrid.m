close all
clear all

% load MRST modules
try
    require ad-core multimodel mrst-gui battery mpfa
catch
    mrstModule add ad-core multimodel mrst-gui battery mpfa
end

%% Parameters 

% number of layers in the spiral
nwindings = 3;

% "radius" at the middle
depth = 2*milli*meter;

% widths of each component ordered as
%  - positive current collector
%  - positive electrode
%  - electrolyte separator 
%  - negative electrode
%  - negative current collector
widthDict = containers.Map( ...
    {'ElectrolyteSeparator', ... 
     'NegativeActiveMaterial', ...
     'NegativeCurrentCollector', ...
     'PositiveActiveMaterial', ...
     'PositiveCurrentCollector'}, ...
    1e-6*[25, 47.8, 15, 42.6, 15]);

% length of the battery
L = 65*milli*meter;

% number of cell in radial direction for each component (same ordering as above).

nrDict = containers.Map( ...
    {'ElectrolyteSeparator', ... 
     'NegativeActiveMaterial', ...
     'NegativeCurrentCollector', ...
     'PositiveActiveMaterial', ...
     'PositiveCurrentCollector'}, ...
    [3, 3, 3, 3, 3]);

% number of cells in the angular direction
nas = 1; 
% number of discretization cells in the longitudonal
nL = 1;

params = struct('nwindings', nwindings, ...
                'depth'    , depth    , ...
                'widthDict', widthDict, ...
                'nrDict'   , nrDict   , ...
                'nas'      , nas      , ...
                'L'        , L        , ...
                'nL'       , nL);


output = flatGrid(params);

G                       = output.G;
tag                     = output.tag;
tagdict                 = output.tagdict;
positiveExtCurrentFaces = output.positiveExtCurrentFaces;
negativeExtCurrentFaces = output.negativeExtCurrentFaces;
thermalExchangeFaces    = output.thermalExchangeFaces;

k = tagdict.keys();
v = tagdict.values();

figure
c = hot(5);
cdict = containers.Map( ...
    {'NegativeCurrentCollector', ...
     'NegativeActiveMaterial', ...
     'ElectrolyteSeparator', ...
     'PositiveActiveMaterial', ...
     'PositiveCurrentCollector'}, ...
    {c(1, :), c(2, :), c(3, :), c(4, :), c(5, :)});


leg = {};
for ind = 1 : tagdict.length
    name = k{ind};
    cdict(name);
    plotGrid(G, v{ind} == tag, 'facecolor', cdict(name));
    leg{end + 1} = name;
end

legend(leg);

for ind = 1 : tagdict.length
    name = k{ind};
    cdict(name);
    figure
    plotGrid(G, v{ind} == tag, 'facecolor', cdict(name));
    title(name);
end


%%

figure
plotGrid(G, 'facecolor', 'none');
plotFaces(G, thermalExchangeFaces);
title('thermal faces')

%% 

figure
plotGrid(G, tag == tagdict('NegativeCurrentCollector'), 'facecolor', 'none', 'edgealpha', 0.1);
plotFaces(G, negativeExtCurrentFaces);
title('NegativeCurrentCollector connection');

figure
plotGrid(G, tag == tagdict('PositiveCurrentCollector'), 'facecolor', 'none', 'edgealpha', 0.1);
plotFaces(G, positiveExtCurrentFaces);
title('PositiveCurrentCollector connection');
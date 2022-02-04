try
    require ad-core multimodel mrst-gui battery mpfa
catch
    mrstModule add ad-core multimodel mrst-gui battery mpfa
end

mrstDebug(20);

%% Indexing of the different cases

%% Tab cases
tabparamscases = {'no tab', '1 tab', '3 tabs'};

%% Cooling cases
uniform_coolingcases = {true, false};

%% battery cases
batterycases = {'18650', '4680'};


%% Choice of case index

batterycase         = batterycases{2};
tabparamscase       = tabparamscases{2};
uniform_coolingcase = uniform_coolingcases{1};

%% Start setup

r0 = 2*milli*meter; 

% widths of each component ordered as
% - positive current collector
% - positive electrode
% - electrolyte separator 
% - negative electrode
% - negative current collector
widths = 1e-6*[25, 64, 15, 57, 15]; 

widthDict = containers.Map( ...
    {'ElectrolyteSeparator',... 
     'NegativeActiveMaterial',...
     'NegativeCurrentCollector',...
     'PositiveActiveMaterial',...
     'PositiveCurrentCollector'},...
    widths); 

nwidths = [widthDict('PositiveActiveMaterial');...
           widthDict('PositiveCurrentCollector');...
           widthDict('PositiveActiveMaterial');...
           widthDict('ElectrolyteSeparator');...
           widthDict('NegativeActiveMaterial');...
           widthDict('NegativeCurrentCollector');...
           widthDict('NegativeActiveMaterial');...
           widthDict('ElectrolyteSeparator')]; 

dr = sum(nwidths);

switch batterycase
  case '18650'
    rOuter    = 18*milli*meter/2;
  case '4680'
    rOuter    = 46*milli*meter/2;
end

dR        = rOuter - r0; 
nwindings = ceil(dR/dr);

%% 
% length of the battery
switch batterycase
  case '18650'
    L = 65*milli*meter; 
  case '4680'
    L = 80*milli*meter; 
end

% number of cell in radial direction for each component (same ordering as above).

nrDict = containers.Map( ...
    {'ElectrolyteSeparator',... 
     'NegativeActiveMaterial',...
     'NegativeCurrentCollector',...
     'PositiveActiveMaterial',...
     'PositiveCurrentCollector'},...
    [3, 3, 3, 3, 3]); 

% number of cells in the angular direction
nas = 20; 
 
% number of discretization cells in the longitudonal
nL = 2; 

tabparams.tabcase = tabparamscase;
switch tabparamscase
  case 'no tab'
  case '1 tab'
    tabparams.width = 3*milli*meter;
  case '3 tabs'
    tabparams.widths = [2*pi*r0/5; 3*milli*meter; 3*milli*meter];
  otherwise
    error('tabparamscase not recognized');
end

params = struct('nwindings'   , nwindings, ...
                'r0'          , r0       , ...
                'widthDict'   , widthDict, ...
                'nrDict'      , nrDict   , ...
                'nas'         , nas      , ...
                'L'           , L        , ...
                'nL'          , nL       , ...
                'tabparams'   , tabparams, ...
                'angleuniform', false); 

output = spiralGrid(params);

G = output.G;
positiveExtCurrentFaces = output.positiveExtCurrentFaces;

%%

close all
plotGrid(G);
plotFaces(G, positiveExtCurrentFaces, 'edgecolor', 'red', 'linewidth', 3);
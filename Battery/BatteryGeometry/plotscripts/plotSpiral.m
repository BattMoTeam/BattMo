try
    require ad-core mrst-gui battery mpfa
catch
    mrstModule add ad-core mrst-gui mpfa
end

mrstDebug(20);

%% Indexing of the different cases

%% Tab cases
tabparamscases = {'no tab', 'aligned tabs', '1 tab', '3 tabs'};

%% Cooling cases
coolingparamscases = {};
coolingparams = struct('top', 10, ...
                       'side', 10);
coolingparamscases{end + 1} = coolingparams;

%% battery cases
batterycases = {'4680', '18650'};

%% Choice of case index
batterycase       = batterycases{1};
tabparamscase     = tabparamscases{2};
coolingparamscase = coolingparamscases{1};

%% Start setup

rInner = 2*milli*meter; 

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
    rOuter = 18*milli*meter/2;
    L      = 65*milli*meter; 
  case '4680'
    rOuter = 46*milli*meter/2;
    L      = 80*milli*meter; 
end

dR        = rOuter - rInner; 
nwindings = ceil(dR/dr);

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
  case 'aligned tabs'
    tabparams.width = 3*milli*meter;
    tabparams.fractions = linspace(0.01, 0.9, 6);
  case '1 tab'
    tabparams.width = 3*milli*meter;
  case '3 tabs'
    tabparams.widths = [2*pi*rInner/5; 3*milli*meter; 3*milli*meter];
  otherwise
    error('tabparamscase not recognized');
end

params = struct('nwindings'   , nwindings, ...
                'rInner'          , rInner       , ...
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
view([-2, -90]);

output.tabwidths


%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}

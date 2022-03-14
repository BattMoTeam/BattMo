try
    require ad-core mrst-gui battery mpfa
catch
    mrstModule add ad-core mrst-gui mpfa
end

clear all
dataSintefDirectory =  '/home/xavier/Sintef-Computer/Matlab/mrst-bitbucket/mrst-core/output/JellyRole';
dataDirectory =  '/home/xavier/Matlab/mrst-bitbucket/mrst-core/output/JellyRole';

tabcases = {true, false};

tabcase = true;

if tabcase
    % tab case
    dataFolder = '8e81eec5fcd8a999d7559f47d1ef17fa';
else
    % no tab case
    dataFolder = 'ce5e6c16a88ba329c04fa8a81c3fe77f';
end

uselocal = true;
if ~uselocal
    dataDirectory = dataSintefDirectory;
end

%%

filename = fullfile(dataDirectory, dataFolder, 'simulationInput.mat');
load(filename);

paramobj = BatteryInputParams(jsonstruct);
gen = SpiralBatteryGenerator(); 

params = options.params;
[paramobj, gen] = gen.updateBatteryInputParams(paramobj, params); 

model = Battery(paramobj, ...
                'use_thermal'        , options.use_thermal, ...
                'use_solid_diffusion', options.use_solid_diffusion);

states = ResultHandler();
states.dataDirectory = dataDirectory;
states.dataFolder = dataFolder;
states = states(:);

%%

UGrids = gen.setupUnRolledGrids(paramobj);

%% 

for ind = 1 : numel(states)
    states{ind} = addProperties(model, states{ind});
end

%%

cartInd = UGrids.NegativeElectrode.ActiveMaterial.mappings.ind;
cartG =  UGrids.NegativeElectrode.ActiveMaterial;

clear newstates
for ind = 1 : numel(states)
    clear state
    state.phi = states{ind}.NegativeElectrode.ActiveMaterial.phi(cartInd);
    state.c = states{ind}.NegativeElectrode.ActiveMaterial.c(cartInd);
    newstates{ind} = state;
end

close all
plotToolbar(cartG, newstates);





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

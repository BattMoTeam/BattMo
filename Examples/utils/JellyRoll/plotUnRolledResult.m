try
    require ad-core multimodel mrst-gui battery mpfa
catch
    mrstModule add ad-core multimodel mrst-gui mpfa
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

cartInd = UGrids.NegativeElectrode.ElectrodeActiveComponent.mappings.ind;
cartG =  UGrids.NegativeElectrode.ElectrodeActiveComponent;

clear newstates
for ind = 1 : numel(states)
    clear state
    state.phi = states{ind}.NegativeElectrode.ElectrodeActiveComponent.phi(cartInd);
    state.c = states{ind}.NegativeElectrode.ElectrodeActiveComponent.c(cartInd);
    newstates{ind} = state;
end

close all
plotToolbar(cartG, newstates);



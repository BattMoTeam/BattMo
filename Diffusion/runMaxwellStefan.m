filename = fullfile(battmoDir(), 'Diffusion', 'maxwellstefan.json');
jsonstruct = parseBattmoJson(filename);
inputparams = MaxwellStefanDiffusionInputParams(jsonstruct);

model = MaxwellStefanGasDiffusion(inputparams);

cgp = model.cgp;
cgt = model.cgt;

cgt.printRootVariables;
cgt.printTailVariables;

return

%% plotting

close all

figure
cgp.plot;

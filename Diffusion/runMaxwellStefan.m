filename = fullfile(battmoDir(), 'Diffusion', 'gasdiffusioncell.json');
jsonstruct = parseBattmoJson(filename);
inputparams = GasDiffusionCellInputParams(jsonstruct);

model = GasDiffusionCell(inputparams);

cgt = model.cgt;
cgp = model.cgp;

cgt.printRootVariables;
cgt.printTailVariables;

return

%% plotting

close all

figure
cgp.plot;

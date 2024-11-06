filename = fullfile(battmoDir(), 'Diffusion', 'Parameters', 'gasdiffusioncell.json');
jsonstruct = parseBattmoJson(filename);
inputparams = GasDiffusionCellInputParams(jsonstruct);

gen = GasDiffusionCellGridGenerator1D(1, 10);
inputparams = gen.updateGasDiffusionCellInputParams(inputparams);

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

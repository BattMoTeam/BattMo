filename = 'CO2membrane/co2membrane.json';
jsonstruct = parseBattmoJson(filename);

inputparams = CO2membraneInputParams(jsonstruct);

gen = CO2membraneGridGenerator();

inputparams = gen.updateInputParams(inputparams);

model = CO2membrane(inputparams);

return

cgp = model.cgp
cgt = model.cgt

close all
figure
cgp.plot;

cgt.printRootVariables;
cgt.printTailVariables;

mrstModule add ad-core mpfa


jsonstring = fileread('/home/xavier/Matlab/Projects/battmo/Electrolyser/Parameters/alkalineElectrolyser.json');
jsonstruct = jsondecode(jsonstring);
paramobj = ElectrolyserInputParams(jsonstruct);


xs = ones(5, 1);
nxs = 5 + (1 : 5)';

gen = ElectrolyserGridGenerator1D(xs, nxs);

paramobj = gen.updateElectrolyserInputParams(paramobj);

inm = 'IonomerMembrane';
her = 'HydrogenEvolutionElectrode';
oer = 'OxygenEvolutionElectrode';
ptl = 'PorousTransportLayer';
exl = 'ExchangeLayer';
ctl = 'CatalystLayer';

model = Electrolyser(paramobj);

doplotgraph = true;
if doplotgraph
    cgt = ComputationalGraphTool(model);
    g = cgt.getComputationalGraph();
    close all
    plot(g);
end


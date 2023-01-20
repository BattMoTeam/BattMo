mrstModule add ad-core mpfa

jsonstring = fileread('/home/xavier/Matlab/Projects/battmo/Electrolyser/Parameters/alkalineElectrolyser.json');
jsonstruct = jsondecode(jsonstring);
paramobj = ElectrolyserInputParams(jsonstruct);

jsonstring = fileread('/home/xavier/Matlab/Projects/battmo/Electrolyser/Parameters/electrolysergeometry1d.json');
jsonstruct = jsondecode(jsonstring);

paramobj = setupElectrolyserGridFromJson(paramobj, jsonstruct);

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


mrstModule add ad-core mpfa matlab_bgl

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

doplotgraph = false;
if doplotgraph
    cgt = ComputationalGraphTool(model);
    g = cgt.getComputationalGraph();
    close all
    plot(g);
    return
end

model = model.validateModel();
[model, initstate] = model.setupBcAndInitialState();

cgt = model.computationalGraph;

total = 1*minute;
n  = 100;
dt = total/n;

controlI = -1e-2; % if negative, H2  is produced
tup = total/2; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time) rampupControl(time, tup, controlI);
control = struct('src', srcfunc);


step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
schedule = struct('control', control, 'step', step);

nls = NonLinearSolver();
nls.verbose = true;
model.verbose = true;

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls); 


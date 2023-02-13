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
cgt = model.computationalGraph;

[model, initstate] = model.setupBcAndInitialState();

% total = 10*hour;
total = 5*minute;
n  = 100;
dt = total/n;

controlI = -30000; % if negative, O2 and H2  are produced

tup = total; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time) rampupControl(time, tup, controlI, 'rampupcase', 'linear');
control = struct('src', srcfunc);

step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
schedule = struct('control', control, 'step', step);

nls = NonLinearSolver();
nls.verbose = false;
nls.errorOnFailure = false;

model.verbose = false;

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls, 'OutputMiniSteps', true);


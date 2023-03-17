mrstModule add ad-core mpfa matlab_bgl

mrstDebug(0);

jsonfilename = '/home/xavier/Matlab/Projects/battmo/Electrolyser/Parameters/alkalineElectrolyser.json';
jsonstruct_base = parseBattmoJson(jsonfilename);

jsonfilename = '/home/xavier/Matlab/Projects/battmo/Electrolyser/Parameters/dissolution.json';
jsonstruct_diss = parseBattmoJson(jsonfilename);

jsonstruct = mergeJsonStructs({jsonstruct_base, ...
                               jsonstruct_diss});
                              
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
total = 36000;
n  = 100;
dt = total/n;
dts = rampupTimesteps(total, dt, 5);

controlI = -30000; % if negative, O2 and H2  are produced

tup = total; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time) rampupControl(time, tup, controlI, 'rampupcase', 'linear');
control = struct('src', srcfunc);

% dts = dts(1 : 45);
step = struct('val', dts, 'control', ones(numel(dts), 1));
schedule = struct('control', control, 'step', step);

nls = NonLinearSolver();
nls.verbose = false;
nls.errorOnFailure = false;

model.verbose = false;

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls, 'OutputMiniSteps', true);

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);
for istate = 1 : numel(states)
    states{istate} = model.addVariables(states{istate});
end
time = cellfun(@(state) state.time, states);
E = cellfun(@(state) state.(oer).(ptl).E, states);
I = cellfun(@(state) state.(oer).(ctl).I, states);

close all

figure
plot(time/hour, E)
xlabel('time [hour]');
ylabel('voltage');

figure
plot(time/hour, -I/(1/(centi*meter)^2));
xlabel('time [hour]');
ylabel('Current [A/cm^2]');

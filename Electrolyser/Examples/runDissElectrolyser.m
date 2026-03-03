mrstModule add ad-core mpfa

mrstDebug(20);

jsonfilename = fullfile('Electrolyser','Parameters','alkalineElectrolyser.json');
jsonstruct_base = parseBattmoJson(jsonfilename);

jsonfilename = fullfile('Electrolyser','Parameters','dissolution.json');
jsonstruct_diss = parseBattmoJson(jsonfilename);

jsonstruct = mergeJsonStructs({jsonstruct_base, ...
                               jsonstruct_diss});

inputparams = ElectrolyserInputParams(jsonstruct);

jsonstring = fileread(fullfile('Electrolyser','Parameters','electrolysergeometry1d.json'));
jsonstruct = battMojsondecode(jsonstring);

inputparams = setupElectrolyserGridFromJson(inputparams, jsonstruct);

inm = 'IonomerMembrane';
her = 'HydrogenEvolutionElectrode';
oer = 'OxygenEvolutionElectrode';
ptl = 'PorousTransportLayer';
exr = 'ExchangeReaction';
ctl = 'CatalystLayer';
dm  = 'DissolutionModel';

model = Electrolyser(inputparams);

doplotgraph = false;
if doplotgraph
    cg = ComputationalGraph(model);
    g = cg.getComputationalGraph();
    close all
    plot(g);
    return
end

model = model.validateModel();

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

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls, 'OutputMiniSteps', true);

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



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}

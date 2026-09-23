el      = 'ElectronicModel';
thermal = 'ThermalModel';
ctrl    = 'Control';

filename = fullfile(battmoDir(), 'Examples', 'Experimental', 'ConductorBlock', 'conductorBlock.json');
jsonstruct = parseBattmoJson(filename);

jsonstruct.(ctrl).Imax = 0.8;

inputparams = ConductorBlockInputParams(jsonstruct);

gen = ConductorBlockGridGenerator();

gen.xlength = 1.1;
gen.ylength = 1.3;
gen.zlength = 1.5;

gen.nx = 5;
gen.ny = 4;
gen.nz = 3;

inputparams = gen.updateGridInputParams(inputparams, jsonstruct);

model = ConductorBlock(inputparams);

initstate = model.setupInitialState();

totalTime = 1*minute;
N = 10;

step = struct('val'    , totalTime/N*ones(N, 1), ...
              'control', ones(N, 1));

control = struct('src', []);

schedule = struct('step'   , step, ...
                  'control', control);

%% Running the simulation

model.verbose = true;
[~, states, report] = simulateScheduleAD(initstate, model, schedule);

%% plotting

close all

state = states{end};
G = model.(el).grid;

figure
plotToolbar(model.(el).grid, state.(el).phi);

figure
plotToolbar(model.(thermal).grid, state.(thermal).T);

state = model.addVariables(state);
figure
plotToolbar(model.(thermal).grid, state.(thermal).jHeatSource);

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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

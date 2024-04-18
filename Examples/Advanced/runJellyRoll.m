if mrstPlatform('octave')
    error('This demo cannot be run from Octave since it does not support the linear solvers');
end

clear
close all

% mrstDebug(20);

% Setup mrst modules
mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

%% We setup the geometrical parameters for a 4680 battery.
%% Those will be gathered in structure spiralparams (see below) and used by SpiralBatteryGenerator to generate the spiral layered geometry of the jelly roll

% Inner radius of the jelly roll
rInner = 2*milli*meter;

% widths of each component ordered as
% - positive current collector
% - positive electrode
% - electrolyte separator
% - negative electrode
% - negative current collector

widths = [25, 64, 15, 57, 15]*micro*meter;

widthDict = containers.Map(...
    {'Separator'               , ...
     'NegativeCoating'         , ...
     'NegativeCurrentCollector', ...
     'PositiveCoating'         , ...
     'PositiveCurrentCollector'}, ...
    widths);

nwidths = [widthDict('PositiveCoating')          ;...
           widthDict('PositiveCurrentCollector') ;...
           widthDict('PositiveCoating')          ;...
           widthDict('Separator')                ;...
           widthDict('NegativeCoating')          ;...
           widthDict('NegativeCurrentCollector') ;...
           widthDict('NegativeCoating')          ;...
           widthDict('Separator')];

dr = sum(nwidths);

% Outer radius of the jelly roll
rOuter = 22.36*milli*meter;
% Height of the jelly roll
L = 80*milli*meter;

dR = rOuter - rInner;
% Computed number of windings
nwindings = ceil(dR/dr);

% number of discretization cells in radial direction for each component.
nrDict = containers.Map({'Separator'               , ...
                         'NegativeCoating'         , ...
                         'NegativeCurrentCollector', ...
                         'PositiveCoating'         , ...
                         'PositiveCurrentCollector'}, ...
                        [3, 3, 3, 3, 3]);

% Number of discretization cells in the angular direction
nas = 10;

% Number of discretization cells in the longitudonal
nL = 3;

% structure that describes the tab setups (see SpiralBatteryGenerator)
tabparams.tabcase   = 'aligned tabs';
tabparams.width     = 3*milli*meter;
tabparams.fractions = linspace(0.01, 0.9, 6);

testing = true;
if testing
    fprintf('We setup a smaller case for quicker testing\n');
    rOuter = rInner + 1*milli*meter;
dR = rOuter - rInner;
% Computed number of windings
nwindings = ceil(dR/dr);
nL = 2;
end

spiralparams = struct('nwindings'   , nwindings, ...
                      'rInner'      , rInner   , ...
                      'widthDict'   , widthDict, ...
                      'nrDict'      , nrDict   , ...
                      'nas'         , nas      , ...
                      'L'           , L        , ...
                      'nL'          , nL       , ...
                      'tabparams'   , tabparams, ...
                      'angleuniform', true);

% The input material parameters given in json format are used to populate the inputparams object.
jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
co      = 'Coating';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

jsonstruct.include_current_collectors = true;

jsonstruct.use_thermal = true;

diffusionModelType = 'full';

jsonstruct.(pe).(co).(am).diffusionModelType = diffusionModelType;
jsonstruct.(ne).(co).(am).diffusionModelType = diffusionModelType;

simcase = 'CCDischarge';

jsonstruct.(ctrl).controlPolicy = simcase;

inputparams = BatteryInputParams(jsonstruct);

% inputparams.(ne).(am).(sd).N = 5;
% inputparams.(pe).(am).(sd).N = 5;

DRate = 0.1;
inputparams.(ctrl).lowerCutoffVoltage = 3;
inputparams.(ctrl).DRate              = DRate;
inputparams.(ctrl).rampupTime         = 0.1/DRate;

% th = 'ThermalModel';
% inputparams.(th).externalHeatTransferCoefficientSideFaces = 100*watt/meter^2;
% inputparams.(th).externalHeatTransferCoefficientTopFaces = 10*watt/meter^2;

gen = SpiralBatteryGenerator();

inputparams = gen.updateBatteryInputParams(inputparams, spiralparams);

model = Battery(inputparams);

%% Setup schedule

DRate = model.(ctrl).DRate;
fac   = 2;
total = 1.1*hour/DRate;
n     = 10;
dt0   = total*1e-6;
times = getTimeSteps(dt0, n, total, fac);

% times = times(times < 4200);

%% We compute the cell capacity, which used to compute schedule from DRate

step = struct('val', diff(times), 'control', ones(numel(times) - 1, 1));

control = model.(ctrl).setupScheduleControl();

schedule  = struct('control', control, 'step', step);

switch simcase

  case 'CCDischarge'

    %% We setup the initial state
    initstate = model.setupInitialState();
    initElytec = 1*mol/litre;
    initstate.(elyte).c = initElytec*ones(model.(elyte).G.getNumberOfCells(), 1);

  case 'CCCharge'

    model.SOC = 0.01;
    initstate = model.setupInitialState();

  otherwise

    error('simcase not recognized')

end

%% Setup the properties of the nonlinear solver
nls = NonLinearSolver();

clear setup

casenumber = 2;
beVerbose = false;

switch casenumber

  case 1

    jsonstruct_solver = parseBattmoJson('Utilities/Linearsolvers/JsonDataFiles/default_direct_linear_solver.json');

  case 2

    switch model.(ne).(co).(am).diffusionModelType

      case  'full'
        
        jsonstruct_solver = parseBattmoJson('Utilities/Linearsolvers/JsonDataFiles/default_linear_solver_setup.json');
        
        if ~beVerbose
            
            jsonstruct_solver.NonLinearSolver.verbose = false;
            jsonstruct_solver.NonLinearSolver.LinearSolver.linearSolverSetup.verbose = 0;
            prcs = jsonstruct_solver.NonLinearSolver.LinearSolver.linearSolverSetup.preconditioners
            for iprc = 1 : numel(prcs)
                prcs(iprc).solver.verbose = 0;
                prcs(iprc).solver.solver.verbose = false;
            end
            jsonstruct_solver.NonLinearSolver.LinearSolver.linearSolverSetup.preconditioners = prcs;

        end
            
      case 'simple'
        
        % jsonstruct_solver = parseBattmoJson('Utilities/Linearsolvers/JsonDataFiles/default_linear_solver_setup_simple_diffusion.json');
        jsonstruct_solver = parseBattmoJson('Utilities/Linearsolvers/JsonDataFiles/temp.json');
        
      otherwise

        error('diffusionModelType not recognized');
        
    end
    
  otherwise

    error('case number not recognized');

end


[model, nls] = setupNonLinearSolverFromJson(model, jsonstruct_solver);

model.verbose = true;

dopacked = true;

if dopacked
    
    % Run simulation
    dataFolder = 'BattMo';
    problem = packSimulationProblem(initstate, model, schedule, dataFolder, ...
                                    'Name'           , 'jellyroll', ...
                                    'NonLinearSolver', nls);
    problem.SimulatorSetup.OutputMinisteps = true;

    clearSimulation = true;
    if clearSimulation
        %% clear previously computed simulation
        clearPackedSimulatorOutput(problem, 'prompt', false);
    end
    simulatePackedProblem(problem);
    [globvars, states, reports] = getPackedSimulatorOutput(problem);

else

    fn = afterStepConvergencePlots(nls);

    [~, states, reports] = simulateScheduleAD(initstate, model, schedule, ...
                                              'NonLinearSolver', nls, ...
                                              'afterStepFn'    , fn);

    reports = reports.ControlstepReports;

end

%% Process output and recover the output voltage and current from the output states.

E = cellfun(@(x) x.Control.E, states);
I = cellfun(@(x) x.Control.I, states);
time = cellfun(@(x) x.time, states);

figure
plot(time, E, 'linewidth', 3);
set(gca, 'fontsize', 18);
title('Cell Voltage / V')
xlabel('time')


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

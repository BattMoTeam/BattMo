function  output = GenerateModelJson(jsonstruct, varargin)
    
    mrstModule add ad-core mrst-gui mpfa
    
    % We define some shorthand names for simplicity.
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    am      = 'ActiveMaterial';
    itf     = 'Interface';
    sd      = 'SolidDiffusion';
    ctrl    = 'Control';
    cc      = 'CurrentCollector';

    paramobj = BatteryInputParams(jsonstruct);

    paramobj = paramobj.validateInputParams();

    [paramobj, gridGenerator] = setupBatteryGridFromJson(paramobj, jsonstruct);

    %% model parameter required for initialization if initializationSetup = "given SOC";
    % The initial state of the model is setup using the model.setupInitialState() method.

    if isfield(jsonstruct, 'initializationSetup')
        initializationSetup = jsonstruct.initializationSetup;
    else
        % default is given SOC
        initializationSetup = "given SOC";
    end

    switch initializationSetup
      case "given SOC"
        paramobj.initT = jsonstruct.initT;
        paramobj.SOC = jsonstruct.SOC;
      case "given input"
        eval(jsonstruct.loadStateCmd);
      otherwise
        error('initializationSetup not recognized');
    end
    
    %%  Initialize the battery model. 
    % The battery model is initialized by sending paramobj to the Battery class
    % constructor. see :class:`Battery <Battery.Battery>`.

    model = Battery(paramobj);
    model.AutoDiffBackend= AutoDiffBackend();
    
    %% Compute the nominal cell capacity and choose a C-Rate
    % The nominal capacity of the cell is calculated from the active materials.
    % This value is then combined with the user-defined C-Rate to set the cell
    % operational current. 

    CRate = model.Control.CRate;

    %% Setup the time step schedule

    total = jsonstruct.TimeStepping.totalTime;
    n = jsonstruct.TimeStepping.N;

    dt = total/n;
    dts = rampupTimesteps(total, dt, 5);
    
    step = struct('val', dts, 'control', ones(numel(dts), 1));

    % we setup the control by assigning a source and stop function.
    % control = struct('CCCV', true); 
    %  !!! Change this to an entry in the JSON with better variable names !!!

    switch model.Control.controlPolicy
      case 'IEswitch'
        tup = 0.1; % rampup value for the current function, see rampupSwitchControl
        srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                    model.Control.Imax, ...
                                                    model.Control.lowerCutoffVoltage);
        % we setup the control by assigning a source and stop function.
        control = struct('src', srcfunc, 'IEswitch', true);
      case 'CCCV'
        control = struct('CCCV', true);
      case 'powerControl'
        control = struct('powerControl', true);
      case 'CC'
        tup = dt;
        srcfunc = @(time) model.(ctrl).rampupControl(time, tup);
        switch model.(ctrl).initialControl
          case 'discharging'
            stopFunc = @(model, state, state_prev) (state.(ctrl).E < model.(ctrl).lowerCutoffVoltage);
          case 'charging'
            stopFunc = @(model, state, state_prev) (state.(ctrl).E > model.(ctrl).upperCutoffVoltage);
          otherwise
            error('initial control not recognized');
        end
        control = struct('CC', true);
        control.src = srcfunc;
        control.stopFunction = stopFunc;
      otherwise
        error('control policy not recognized');
    end
    switch initializationSetup
      case "given SOC"
        initstate = model.setupInitialState();
      case "given input"
        % allready handled
      otherwise
        error('initializationSetup not recognized');
    end
    % This control is used to set up the schedule
    schedule = struct('control', control, 'step', step); 
    output = struct('model', model, ...
                    'schedule',schedule,...
                    'initstate', initstate);
end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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



% 
% 
%     switch initializationSetup
%       case "given SOC"
%         initstate = model.setupInitialState();
%       case "given input"
%         % allready handled
%       otherwise
%         error('initializationSetup not recognized');
%     end
%     
%     
%     %% Setup the properties of the nonlinear solver 
%     nls = NonLinearSolver();
% 
%     % load from json file or use default value. 
%     if ~isfield(jsonstruct, 'NonLinearSolver')
%         
%         % setup default values
%         jsonstruct.NonLinearSolver.maxIterations = 10;
%         jsonstruct.NonLinearSolver.nonlinearTolerance = 1e-5;
%         jsonstruct.NonLinearSolver.verbose = false;
%         
%         linearSolverSetup.library = "matlab";
%         linearSolverSetup.method = "direct";
% 
%     else
% 
%         linearSolverSetup = jsonstruct.NonLinearSolver.LinearSolver.linearSolverSetup;
%         
%     end
%     
%     nls.maxIterations = jsonstruct.NonLinearSolver.maxIterations;
%     if isfield(jsonstruct.NonLinearSolver, 'verbose')
%         nls.verbose = jsonstruct.NonLinearSolver.verbose;
%     else
%         nls.verbose = false;
%     end
%     
%     % Change default behavior of nonlinear solver, in case of error
%     nls.errorOnFailure = false;
%     nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
% 
%     nls.LinearSolver = BatteryLinearSolver('linearSolverSetup', linearSolverSetup);
%     
%     % Change default tolerance for nonlinear solver
%     % For the moment, this is hard-coded here
%     if ~isempty(model.Control.Imax)
%         model.nonlinearTolerance = 1e-3*model.Control.Imax;
%     else
%         model.nonlinearTolerance = jsonstruct.NonLinearSolver.nonlinearTolerance;
%     end
%     
%     % Set verbosity
%     model.verbose = jsonstruct.NonLinearSolver.verbose;
% 
%     if isfield(linearSolverSetup, 'reduction') && linearSolverSetup.reduction.doReduction
%         model = model.setupSelectedModel('reduction', linearSolverSetup.reduction);
%     end
%     
%         
%     %% Run the simulation
%     if isfield(jsonstruct, 'Output') && isfield(jsonstruct.Output, 'saveOutput') && jsonstruct.Output.saveOutput
%         saveOptions = jsonstruct.Output.saveOptions;
%         
%         outputDirectory = saveOptions.outputDirectory;
%         name            = saveOptions.name;
%         clearSimulation = saveOptions.clearSimulation;
%         
%         problem = packSimulationProblem(initstate, model, schedule, [], ...
%                                         'Directory'      , outputDirectory, ...
%                                         'Name'           , name      , ...
%                                         'NonLinearSolver', nls);
%         problem.SimulatorSetup.OutputMinisteps = true; 
% 
%         if clearSimulation
%             %% clear previously computed simulation
%             clearPackedSimulatorOutput(problem, 'prompt', false);
%         end
%         simulatePackedProblem(problem);
%         [globvars, states, reports] = getPackedSimulatorOutput(problem);
%         
%     else
%         [~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);
%     end
% 
%     %% Process output and recover the output voltage and current from the output states.
%     ind = cellfun(@(x) not(isempty(x)), states); 
%     states = states(ind);
%     
%     E    = cellfun(@(state) state.Control.E, states); 
%     I    = cellfun(@(state) state.Control.I, states);
%     time = cellfun(@(state) state.time, states); 
% 
%     output = struct('model'    , model    , ...
%                     'schedule' , schedule , ... 
%                     'initstate', initstate, ...
%                     'time'     , time     , ... % Unit : s
%                     'E'        , E        , ... % Unit : V
%                     'I'        , I); ... % Unit : A
% 
%     output.states = states;
%     
%     if isfield(jsonstruct, 'Output') ...
%         && isfield(jsonstruct.Output, 'variables') ...
%         && any(ismember({'energy', 'energyDensity', 'specificEnergy'}, jsonstruct.Output.variables))
%         
%         % TODO : We could fine tuned output (at the moment we output all the possible extrav variables)
%         mass = computeCellMass(model);
%         vol = sum(model.G.cells.volumes);
%         
%         [E, I, energyDensity, specificEnergy, energy] = computeEnergyDensity(E, I, time, vol, mass);
% 
%         output.E              = E;
%         output.I              = I;
%         output.energy         = energy;          % Unit : J
%         output.energyDensity  = energyDensity;   % Unit : J/L
%         output.specificEnergy = specificEnergy;  % Unit : J/kg
%         
%     end

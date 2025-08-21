function  output = runBatteryJson(jsonstruct, varargin)

    opt = struct('runSimulation'       , true , ...
                 'includeGridGenerator', false, ...
                 'initstate'           , []   , ...
                 'validateJson'        , false, ...
                 'verbose'             , true);
    opt = merge_options(opt, varargin{:});

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

    %% Validate
    if opt.validateJson
        validateJsonStruct(jsonstruct);
    end

    %% model parameter required for initialization if initializationSetup = "given SOC";
    % The initial state of the model is setup using the model.setupInitialState() method.

    if ~isempty(opt.initstate)
        jsonstruct = setJsonStructField(jsonstruct, {'initializationSetup'}, 'given matlab object', 'handleMisMatch', 'warn');
    else
        jsonstruct = setDefaultJsonStructField(jsonstruct, {'initializationSetup'}, 'given SOC');
    end

    %%  Initialize the battery model.
    % The battery model is setup using :battmo:`setupModelFromJson`

    [model, inputparams, jsonstruct, gridGenerator] = setupModelFromJson(jsonstruct);

    initializationSetup = getJsonStructField(jsonstruct, {'initializationSetup'});
    
    switch initializationSetup
      case "given SOC"
        % nothing to do
      case "given input"
        eval(jsonstruct.loadStateCmd);
      case "given matlab object"
        initstate = opt.initstate;
      otherwise
        error('initializationSetup not recognized');
    end

    %% Setup the time step schedule
    %

    if isAssigned(jsonstruct, 'TimeStepping')
        timeSteppingParams = jsonstruct.TimeStepping;
    else
        timeSteppingParams = [];
    end

    step    = model.(ctrl).setupScheduleStep(timeSteppingParams);
    control = model.(ctrl).setupScheduleControl();

    % This control is used to set up the schedule
    schedule = struct('control', control, 'step', step);

    switch initializationSetup
      case "given SOC"
        if isfield(jsonstruct, 'Electrolyte') && isfield(jsonstruct.Electrolyte, 'initialConcentration')
            jsonstructInit.Electrolyte.initialConcentration = jsonstruct.Electrolyte.initialConcentration;
        else
            jsonstructInit = [];
        end
        initstate = model.setupInitialState(jsonstructInit);
      case {"given input", 'given matlab object'}
        % allready handled
      otherwise
        error('initializationSetup not recognized');
    end

    [model, nls, jsonstruct] = setupNonLinearSolverFromJson(model, jsonstruct);

    simsetupInput = struct('model'          , model    , ...
                           'initstate'      , initstate, ...
                           'schedule'       , schedule , ...
                           'NonLinearSolver', nls);

    simsetup = SimulationSetup(simsetupInput);
    
    model.verbose = opt.verbose;

    %% Run the simulation
    %

    if opt.runSimulation

        if isfield(jsonstruct, 'Output') && isfield(jsonstruct.Output, 'saveOutput') && jsonstruct.Output.saveOutput
            saveOptions = jsonstruct.Output.saveOptions;

            outputDirectory = saveOptions.outputDirectory;
            name            = saveOptions.name;
            clearSimulation = saveOptions.clearSimulation;

            problem = packSimulationProblem(initstate, model, schedule, [], ...
                                            'Directory'      , outputDirectory, ...
                                            'Name'           , name      , ...
                                            'NonLinearSolver', nls);
            problem.SimulatorSetup.OutputMinisteps = simsetup.OutputMinisteps;

            if clearSimulation
                %% clear previously computed simulation
                clearPackedSimulatorOutput(problem, 'prompt', false);
            end
            simulatePackedProblem(problem);
            [globvars, states, reports] = getPackedSimulatorOutput(problem);

        else
            [states, globvars, reports] = simsetup.run();
        end
    else
        output = struct('model'      , model      , ...
                        'inputparams', inputparams, ...
                        'simsetup'   , simsetup   , ...
                        'schedule'   , schedule);
        if opt.includeGridGenerator
            output.gridGenerator = gridGenerator;
        end

        return;
    end

    %% Process output and recover the output voltage and current from the output states.
    ind = cellfun(@(x) not(isempty(x)), states);
    states   = states(ind);
    globvars = globvars(ind);

    E    = cellfun(@(state) state.Control.E, states);
    I    = cellfun(@(state) state.Control.I, states);
    time = cellfun(@(state) state.time, states);

    output = struct('model'      , model      , ...
                    'inputparams', inputparams, ...
                    'simsetup'   , simsetup   , ...
                    'jsonstruct' , jsonstruct , ...
                    'time'       , time       , ... % Unit : s
                    'E'          , E          , ... % Unit : V
                    'I'          , I); ... % Unit : A

    output.globvars = globvars;
    output.states   = states;

    if isAssigned(jsonstruct, {'Output', 'variables'}) ...
        && any(ismember({'energy', 'energyDensity', 'specificEnergy'}, jsonstruct.Output.variables))

        if ismember('specificEnergy', jsonstruct.Output.variables)
            mass = computeCellMass(model);
        else
            mass = [];
        end
        vol = sum(model.G.getVolumes());

        [energyDensity, specificEnergy, energy] = computeEnergyDensity(E, I, time, vol, mass);

        output.energy         = energy;          % Unit : J
        output.energyDensity  = energyDensity;   % Unit : J/L
        output.specificEnergy = specificEnergy;  % Unit : J/kg

    end

    if opt.includeGridGenerator
        output.gridGenerator = gridGenerator;
    end

end



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

function  output = runBatteryJson(jsonstruct, varargin)
    
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

    switch jsonstruct.Geometry.case

      case '1D'

        gen = BatteryGenerator1D();
        
        xlength = gen.xlength;
        xlength(2) = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        xlength(3) = jsonstruct.Electrolyte.Separator.thickness;
        xlength(4) = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
        if paramobj.NegativeElectrode.include_current_collectors
            xlength(1) = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        end
        if paramobj.PositiveElectrode.include_current_collectors
            xlength(5) = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
        end
        if isfield(jsonstruct.Geometry, 'faceArea')
            gen.faceArea = jsonstruct.Geometry.faceArea;
        end

        gen.sepnx  = jsonstruct.NegativeElectrode.ActiveMaterial.N;
        gen.nenx   = jsonstruct.Electrolyte.Separator.N;
        gen.penx   = jsonstruct.PositiveElectrode.ActiveMaterial.N;
        
        % Now, we update the paramobj with the properties of the mesh. 
        paramobj = gen.updateBatteryInputParams(paramobj);

      case {'2D-demo', '3d-demo'}

        error('not yet implemented');

      case 'jellyRoll'

        % Prepare input for SpiralBatteryGenerator.updateBatteryInputParams using json input
        
        widthDict = containers.Map();
        widthDict('ElectrolyteSeparator')     = jsonstruct.Electrolyte.Separator.thickness;
        widthDict('NegativeActiveMaterial')   = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        widthDict('NegativeCurrentCollector') = jsonstruct.NegativeElectrode.CurrentCollector.thickness;
        widthDict('PositiveActiveMaterial')   = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
        widthDict('PositiveCurrentCollector') = jsonstruct.PositiveElectrode.CurrentCollector.thickness;
        
        nwidths = [widthDict('PositiveActiveMaterial');...
                   widthDict('PositiveCurrentCollector');...
                   widthDict('PositiveActiveMaterial');...
                   widthDict('ElectrolyteSeparator');...
                   widthDict('NegativeActiveMaterial');...
                   widthDict('NegativeCurrentCollector');...
                   widthDict('NegativeActiveMaterial');...
                   widthDict('ElectrolyteSeparator')]; 
        dr = sum(nwidths);
        
        rOuter = jsonstruct.Geometry.rOuter;
        rInner = jsonstruct.Geometry.rInner;
        L      = jsonstruct.Geometry.L;
        nL     = jsonstruct.Geometry.nL;
        nas    = jsonstruct.Geometry.nas;

        nrDict = containers.Map();
        nrDict('ElectrolyteSeparator')     = jsonstruct.Electrolyte.Separator.N;
        nrDict('NegativeActiveMaterial')   = jsonstruct.NegativeElectrode.ActiveMaterial.N;
        nrDict('NegativeCurrentCollector') = jsonstruct.NegativeElectrode.CurrentCollector.N;
        nrDict('PositiveActiveMaterial')   = jsonstruct.PositiveElectrode.ActiveMaterial.N;
        nrDict('PositiveCurrentCollector') = jsonstruct.PositiveElectrode.CurrentCollector.N;

        % compute numbers of winding (this is input for spiralGrid) from outer and inner radius
        nwidths = [widthDict('PositiveActiveMaterial');...
                   widthDict('PositiveCurrentCollector');...
                   widthDict('PositiveActiveMaterial');...
                   widthDict('ElectrolyteSeparator');...
                   widthDict('NegativeActiveMaterial');...
                   widthDict('NegativeCurrentCollector');...
                   widthDict('NegativeActiveMaterial');...
                   widthDict('ElectrolyteSeparator')]; 
        dr = sum(nwidths);dR = rOuter - rInner; 
        % Computed number of windings
        nwindings = ceil(dR/dr);

        tabparams.NegativeElectrode = jsonstruct.NegativeElectrode.CurrentCollector.tabparams;
        tabparams.PositiveElectrode = jsonstruct.PositiveElectrode.CurrentCollector.tabparams;
        
        spiralparams = struct('nwindings'   , nwindings, ...
                              'rInner'      , rInner   , ...
                              'widthDict'   , widthDict, ...
                              'nrDict'      , nrDict   , ...
                              'nas'         , nas      , ...
                              'L'           , L        , ...
                              'nL'          , nL       , ...
                              'tabparams'   , tabparams, ...
                              'angleuniform', true); 
        
        gen = SpiralBatteryGenerator(); 
        paramobj = gen.updateBatteryInputParams(paramobj, spiralparams);

      otherwise
        
        error('Geometry case not recognized')
        
    end

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
        error('interface not yet implemented');
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
    
    step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

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
        tup = 0.1;
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

    % This control is used to set up the schedule
    schedule = struct('control', control, 'step', step); 


    switch initializationSetup
      case "given SOC"
        initstate = model.setupInitialState();
      case "given input"
        error('interface not yet implemented');
      otherwise
        error('initializationSetup not recognized');
    end
    
    
    %% Setup the properties of the nonlinear solver 
    nls = NonLinearSolver();

    % load from json file or use default value. 
    if ~isfield(jsonstruct, 'NonLinearSolver')
        
        % setup default values
        jsonstruct.NonLinearSolver.maxIterations = 10;
        jsonstruct.NonLinearSolver.nonlinearTolerance = 1e-5;
        jsonstruct.NonLinearSolver.verbose = false;
        
        linearSolverSetup.library = "matlab";
        linearSolverSetup.method = "direct";

    else

        linearSolverSetup = jsonstruct.NonLinearSolver.LinearSolver.linearSolverSetup;
        
    end
    
    nls.maxIterations = jsonstruct.NonLinearSolver.maxIterations;
    if isfield(jsonstruct.NonLinearSolver, 'verbose')
        nls.verbose = jsonstruct.NonLinearSolver.verbose;
    else
        nls.verbose = false;
    end
    
    % Change default behavior of nonlinear solver, in case of error
    nls.errorOnFailure = false;
    nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);

    nls.LinearSolver = BatteryLinearSolver('linearSolverSetup', linearSolverSetup);
    
    % Change default tolerance for nonlinear solver
    % For the moment, this is hard-coded here
    if ~isempty(model.Control.Imax)
        model.nonlinearTolerance = 1e-3*model.Control.Imax;
    else
        model.nonlinearTolerance = jsonstruct.NonLinearSolver.nonlinearTolerance;
    end
    
    % Set verbosity
    model.verbose = false;

    if isfield(linearSolverSetup, 'reduction') && linearSolverSetup.reduction.doReduction
        model = model.setupSelectedModel('reduction', linearSolverSetup.reduction);
    end
    
        
    %% Run the simulation
    if isfield(jsonstruct, 'Output') && isfield(jsonstruct.Output, 'saveOutput') && jsonstruct.Output.saveOutput
        saveOptions = jsonstruct.Output.saveOptions;
        dataFolder      = saveOptions.dataFolder;
        name            = saveOptions.name;
        resetSimulation = saveOptions.resetSimulation;
        
        problem = packSimulationProblem(initstate, model, schedule, dataFolder, ...
                                        'Name'           , name, ...
                                        'NonLinearSolver', nls);
        problem.SimulatorSetup.OutputMinisteps = true; 

        if resetSimulation
            %% clear previously computed simulation
            clearPackedSimulatorOutput(problem, 'prompt', false);
        end
        simulatePackedProblem(problem);
        [globvars, states, reports] = getPackedSimulatorOutput(problem);
        
    else
        [~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);
    end

    %% Process output and recover the output voltage and current from the output states.
    ind = cellfun(@(x) not(isempty(x)), states); 
    states = states(ind);
    
    E    = cellfun(@(state) state.Control.E, states); 
    I    = cellfun(@(state) state.Control.I, states);
    time = cellfun(@(state) state.time, states); 

    output = struct('model', model, ...
                    'time' , time , ... % Unit : s
                    'E'    , E    , ... % Unit : V
                    'I'    , I); ... % Unit : A

    output.states = states;
    
    if isfield(jsonstruct, 'Output') && any(strcmp({"energyDensity", "specificEnergy"}, jsonstruct.Output))
        
        % We could fine tuned output
        mass = computeCellMass(model);
        vol = sum(model.G.cells.volumes);
        
        [E, I, energyDensity, specificEnergy, energy] = computeEnergyDensity(E, I, time, vol, mass);

        output.E              = E;
        output.I              = I;
        output.energy         = energy;          % Unit : J
        output.energyDensity  = energyDensity;   % Unit : J/L
        output.specificEnergy = specificEnergy;  % Unit : J/kg
        
    end

end

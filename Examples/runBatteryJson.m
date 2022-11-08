function  output = runBatteryJson(jsonInput)

    if ischar(jsonInput)
        % load file
        jsonstruct = parseBattmoJson(jsonInput);
    elseif isstruct(jsonInput)
        jsonstruct = jsonInput;
    else
        error('jsonInput should be either a filename or a jsonstruct (formatted as output of MATLAB jsondecode)');
    end
    
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

    diffusionModelType = 'simple';
    jsonstruct.(pe).(am).diffusionModelType = diffusionModelType;
    jsonstruct.(ne).(am).diffusionModelType = diffusionModelType;

    jsonstruct.use_particle_diffusion = false;
    
    jsonstruct.use_thermal = false;
    
    paramobj = BatteryInputParams(jsonstruct);

    paramobj.(ne).(am).InterDiffusionCoefficient = 0;
    paramobj.(pe).(am).InterDiffusionCoefficient = 0;

    
    if strcmp(diffusionModelType, 'full')
        paramobj.(ne).(am).(sd).N = 5;
        paramobj.(pe).(am).(sd).N = 5;
    end

    paramobj = paramobj.validateInputParams();

    use_cccv = false;
    if use_cccv
        cccvstruct = struct( 'controlPolicy'     , 'CCCV',  ...
                             'CRate'             , 1         , ...
                             'lowerCutoffVoltage', 2.4       , ...
                             'upperCutoffVoltage', 4.1       , ...
                             'dIdtLimit'         , 0.01      , ...
                             'dEdtLimit'         , 0.01);
        cccvparamobj = CcCvControlModelInputParams(cccvstruct);
        paramobj.Control = cccvparamobj;
    end


    %% Setup the geometry and computational mesh
    % Here, we setup the 1D computational mesh that will be used for the
    % simulation. The required discretization parameters are already included
    % in the class BatteryGenerator1D.

    switch jsonstruct.format

      case '1D'

        gen = BatteryGenerator1D();
        
        xlength = gen.xlength;
        xlength(2) = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        xlength(3) = jsonstruct.Electrolyte.Separator.thickness;
        xlength(4) = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
        gen.xlength = xlength;

        gen.sepnx  = 4;
        gen.nenx   = 4;
        gen.penx   = 4;
        
        gen.faceArea = jsonstruct.faceArea;
        
        % Now, we update the paramobj with the properties of the mesh. 
        paramobj = gen.updateBatteryInputParams(paramobj);

      otherwise
        error('format name not recognized')
    end
    

    %%  Initialize the battery model. 
    % The battery model is initialized by sending paramobj to the Battery class
    % constructor. see :class:`Battery <Battery.Battery>`.

    paramobj.Control.lowerCutoffVoltage = 3;

    model = Battery(paramobj);
    model.AutoDiffBackend= AutoDiffBackend();


    %% Compute the nominal cell capacity and choose a C-Rate
    % The nominal capacity of the cell is calculated from the active materials.
    % This value is then combined with the user-defined C-Rate to set the cell
    % operational current. 

    CRate = model.Control.CRate;

    %% Setup the time step schedule 
    % Smaller time steps are used to ramp up the current from zero to its
    % operational value. Larger time steps are then used for the normal
    % operation.
    switch model.(ctrl).controlPolicy
      case 'CCCV'
        total = 3.5*hour/CRate;
      case 'IEswitch'
        total = 1.4*hour/CRate;
      otherwise
        error('control policy not recognized');
    end

    n  = 50;
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
      otherwise
        error('control policy not recognized');
    end

    % This control is used to set up the schedule
    schedule = struct('control', control, 'step', step); 

    %% Setup the initial state of the model
    % The initial state of the model is setup using the model.setupInitialState() method.

    initstate = model.setupInitialState(); 

    %% Setup the properties of the nonlinear solver 
    nls = NonLinearSolver();

    % Change default maximum iteration number in nonlinear solver
    nls.maxIterations = 10;
    % Change default behavior of nonlinear solver, in case of error
    nls.errorOnFailure = false;
    nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
    % Change default tolerance for nonlinear solver
    model.nonlinearTolerance = 1e-3*model.Control.Imax;
    % Set verbosity
    model.verbose = false;

    %% Run the simulation
    [~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

    %% Process output and recover the output voltage and current from the output states.
    ind = cellfun(@(x) not(isempty(x)), states); 
    states = states(ind);
    
    E    = cellfun(@(state) state.Control.E, states); 
    I    = cellfun(@(state) state.Control.I, states);
    time = cellfun(@(state) state.time, states); 

    mass = computeCellMass(model);
    vol = sum(model.G.cells.volumes);
    
    [E, I, energyDensity, specificEnergy, energy] = computeEnergyDensity(E, I, time, vol, mass);

    output = struct('model'         , model         , ...
                    'E'             , E             , ... % Unit : V
                    'I'             , I             , ... % Unit : A
                    'energy'        , energy        , ... % Unit : Wh
                    'energyDensity' , energyDensity , ... % Unit : Wh/L
                    'specificEnergy', specificEnergy  ... % Unit : Wh/kg
                   );
    output.states = states;

end

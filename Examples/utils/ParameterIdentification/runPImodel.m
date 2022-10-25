function [states, model, schedule, initstate] = runPImodel(jsonstruct, varargin)

% Utility function for running a simple model for parameter identification

    opt = struct('family', '1D', ...
                 'fac', 1, ...
                 'dtFac', 1, ...
                 'chen', false);

    opt = merge_options(opt, varargin{:});

    % Define shorthand names for simplicity.
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    am      = 'ActiveMaterial';
    cc      = 'CurrentCollector';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    itf     = 'Interface';
    sd      = 'SolidDiffusion';
    ctrl    = 'Control';

    % Parse json struct
    paramobj = BatteryInputParams(jsonstruct);

    % Setup the geometry
    switch lower(opt.family)
      case '1d'
        gen = BatteryGenerator1D();
      case '2d'
        gen = BatteryGenerator2D();
      case '3d'
        gen = BatteryGenerator3D();
      case {'spiral', 'coincell'}
        error('Battery geometry family %s is not yet implemented', opt.family);
      otherwise
        error('Battery geometry family %s is not implemented', opt.family);
    end

    if opt.fac > 1
        gen.fac = opt.fac;
        gen = gen.applyResolutionFactors();
    end

    % Update parameters with mesh properties
    paramobj = gen.updateBatteryInputParams(paramobj);

    % Initialize battery model
    model = Battery(paramobj);

    % Set AD backend
    model.AutoDiffBackend = AutoDiffBackend();

    % Get C-rate
    CRate = model.Control.CRate;

    % Setup schedule parameters
    switch model.(ctrl).controlPolicy
      case 'CCCV'
        total = 3.5*hour/CRate;
      case 'IEswitch'
        total = 1.2*hour/CRate;
      otherwise
        error('control policy %s is not recognized', model.(ctrl).controlPolicy);
    end
    n    = 80 * opt.dtFac;
    dt   = total*0.7/n;
    step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

    % Setup control
    switch model.(ctrl).controlPolicy
      case 'IEswitch'
        % Rampup value for the current function (see rampupSwitchControl)
        tup = 0.1;
        srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                    model.Control.Imax, ...
                                                    model.Control.lowerCutoffVoltage);

        % Setup the control by assigning a source and stop function
        control = struct('src', srcfunc, 'IEswitch', true);
      case 'CCCV'
        control = struct('CCCV', true);
      otherwise
        error('control policy %s is not recognized', model.(ctrl).controlPolicy);
    end

    % Setup schedule
    schedule = struct('control', control, 'step', step);

    % Setup initial state
    if opt.chen
        c_ne = 29.866*mol/litre; % initial concentration at negative electrode
        c_pe = 17.038*mol/litre; % initial concentration at positive electrode
        initstate = initStateChen2020(model, c_ne, c_pe);
    else
        initstate = model.setupInitialState();
    end

    % Setup nonlinear solver
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

    % Run the simulation
    [~, states, ~] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

end

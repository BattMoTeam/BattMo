function output = computeCellEnergyGivenCrate(model, CRate, varargin)
% Given a model and a CRate, compute the produced energy
% The output consists of the fields
% - energy
% - dischargeFunction
% - E    % Voltage output (as function of time)
% - I    % Current output (as function of time)
% - time % time output

    opt = struct('computationType'   , 'sharp', ...
                 'lowerCutoffVoltage', []     , ...
                 'cutoffparams'      , []);
    opt = merge_options(opt, varargin{:});

    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    am    = 'ActiveMaterial';
    elyte = 'Electrolyte';
    itf   = 'Interface';
    sd    = 'SolidDiffusion';
    ctrl  = 'Control';
    sep   = 'Separator';

    if isempty(opt.lowerCutoffVoltage)
        lowerCutoffVoltage = model.Control.lowerCutoffVoltage;
    else
        lowerCutoffVoltage = opt.lowerCutoffVoltage;
    end

    assert(strcmp(model.Control.controlPolicy, "CCDischarge"), 'The model should be setup with CCDischarge control');

    capacity = computeCellCapacity(model);
    Imax = (capacity/hour)*CRate;

    model.Control.CRate       = CRate;
    model.Control.Imax        = Imax;
    model.Control.useCVswitch = true;

    step    = model.Control.setupScheduleStep();
    control = model.Control.setupScheduleControl();
    
    schedule = struct('control', control, 'step', step);

    %% Setup the initial state of the model

    state0 = model.setupInitialState();

    nls = NonLinearSolver();
    nls.maxIterations  = 10;
    nls.errorOnFailure = false;

    model.nonlinearTolerance = 1e-5*model.Control.Imax;
    model.verbose = false;

    [~, states, report] = simulateScheduleAD(state0, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

    ind = cellfun(@(x) not(isempty(x)), states);
    states = states(ind);
    E = cellfun(@(x) x.(ctrl).E, states);
    I = cellfun(@(x) x.(ctrl).I, states);
    time = cellfun(@(x) x.time, states);

    switch opt.computationType

      case 'smooth'
        % we use a regularization function

        if isempty(opt.cutoffparams)
            cutoffparams.E0 = lowerCutoffVoltage;
            cutoffparams.alpha = 100;
        else
            cutoffparams = opt.cutoffparams;
        end

        cutoff = @(E) cutoffGeneric(E, cutoffparams);

        energy = sum(I.*E.*cutoff(E).*dt);

        E = E.*cutoff(E);

      case 'sharp'
        % We detect the lowercutoffvoltage point

        ind = find(E <= lowerCutoffVoltage, 1, 'first');
        ind = (1 : ind)';

        I  = I(ind);
        E  = E(ind);
        dt = time(ind);
        dt = diff(dt);
        dt = [time(1); dt];

        energy = sum(I.*E.*dt);

      otherwise

        error('computation type not recognized');

    end

    soc = cumsum(I.*dt)/capacity;

    dischargeFunction = @(s) interp1(soc, E, s, 'linear');

    output = struct('energy'           , energy           , ...
                    'dischargeFunction', dischargeFunction, ...
                    'E'                , E                , ...
                    'I'                , I                , ...
                    'time'             , cumsum(dt));

end


function y = cutoffGeneric(E, cutoffparams)

    E0    = cutoffparams.E0;
    alpha = cutoffparams.alpha;

    y = 0.5 + 1/pi*atan(alpha*(E  - E0));

end

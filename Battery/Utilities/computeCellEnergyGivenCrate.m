function output = computeCellEnergyGivenCrate(model, CRate, varargin)
% Given a model and a CRate, compute the produced energy
% The output consists of the fields
% - energy
% - dischargeFunction
% - E    % Voltage output 
% - I    % Current output
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

    capacity = computeCellCapacity(model);
    Imax = (capacity/hour)*CRate;
    
    model.Control.CRate = CRate;
    model.Control.Imax  = Imax;
    
    assert(strcmp(model.Control.controlPolicy, "IEswitch"), 'The model should be setup with IEswitch control');
    
    totalTime = 1.1*hour/CRate;

    % rampup stage
    n  = 10;
    dt = [];
    dt = [dt; repmat(1e-4, n, 1).*1.5.^(1 : n)'];

    % discharge stage
    n     = 100;
    dt    = [dt; repmat(totalTime/n, n, 1)];

    times = [0; cumsum(dt)];
    dt    = diff(times);
    step  = struct('val', dt, 'control', ones(numel(dt), 1));

    tup = 0.1; % rampup value for the current function, see rampupSwitchControl
    srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                model.Control.Imax, ...
                                                model.Control.lowerCutoffVoltage);

    control = struct('src', srcfunc, 'IEswitch', true);


    schedule = struct('control', control, 'step', step); 

    %% Setup the initial state of the model

    state0 = model.setupInitialState(); 

    nls = NonLinearSolver(); 
    nls.maxIterations  = 10; 
    nls.errorOnFailure = false; 
    
    model.nonlinearTolerance = 1e-5*model.Control.Imax;
    model.verbose = false;

    [wellSols, states, report] = simulateScheduleAD(state0, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

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
        dt = dt(ind);
        
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

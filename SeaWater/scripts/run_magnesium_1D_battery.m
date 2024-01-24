function output = run_magnesium_1D_battery(input, varargin)

    input_default = struct('precipitation', true);

    if ~isempty(input)
        fds = fieldnames(input);
        vals = cellfun(@(fd) input.(fd), fds, 'un', false);
        input = horzcat(fds, vals);
        input = reshape(input', [], 1);
        input = merge_options(input_default, input{:});
    else
        input = input_default;
    end

    opt = struct('runSimulation'  , true , ...
                 'clearSimulation', false, ...
                 'initState'      , []   , ...
                 'directory'      , 'magnesium');
    opt = merge_options(opt, varargin{:});

    if isfield(input, 'directory')
        opt.directory = input.directory;
    end

    %% Load json input
    jsonstruct = parseBattmoJson('SeaWater/json_inputs/magnesium_battery.json');

    if input.precipitation
        jsonstruct.include_precipitation = true;
    else
        jsonstruct.include_precipitation = false;
    end

    inputparams = SeaWaterBatteryInputParams(jsonstruct);

    elyte = 'Electrolyte';
    ct    = 'Cathode';
    ctam  = 'CathodeActiveMaterial';
    an    = 'Anode';
    anam  = 'AnodeActiveMaterial';

    %% Setup geometry

    gen = SeaWaterBatteryGeneratorP2D();
    % We increase the resolution
    gen.fac = 20;
    gen = gen.applyResolutionFactors();

    inputparams = gen.updateBatteryInputParams(inputparams);

    %% Setup model

    model = MagnesiumBattery(inputparams);
    model = model.setupComputationalGraph(); % needed later to evaluate values using evalVarName (we want only to setup it once)

    %% Setup Electrolyte initial content

    pH = 7;
    Cl = 500*mol/meter^3;
    Mg = 1*mol/meter^3;

    totals = [Cl, Mg];
    totalnames = {'Cl', 'Mg'};

    melyte = model.Electrolyte;

    qpdict    = melyte.qpdict;
    spdict    = melyte.spdict;
    logspdict = melyte.logspdict;
    nqp       = melyte.nqp;
    nlogsp    = melyte.nlogsp;

    qpcs = cell(1, nqp);
    pcs = cell(1, nlogsp);

    qpcs{qpdict('Cl')}  = Cl;
    qpcs{qpdict('Mg')}  = Mg;
    qpcs{qpdict('HOH')} = 10^(-pH)*(mol/litre);
    pcs{logspdict('H+')}   = log(10^(-pH)*(mol/litre));
    pcs{logspdict('Cl-')}  = log(qpcs{qpdict('Cl')});
    pcs{logspdict('Mg+2')} = log(qpcs{qpdict('Mg')});

    stateelyteguess.qpcs = qpcs;
    stateelyteguess.pcs  = pcs;

    aqmodel = AqueousModel(melyte, totals, totalnames, pH);
    stateelyte = aqmodel.solveAqueousMixture(stateelyteguess);

    nc = model.(elyte).G.getNumberOfCells();
    state.(elyte).qpcs = cellfun(@(val) val*ones(nc, 1), stateelyte.qpcs, 'uniformoutput', false);
    state.(elyte).pcs  = cellfun(@(val) val*ones(nc, 1), stateelyte.pcs , 'uniformoutput', false);
    for ind = 1 : numel(stateelyte.cs)
        val =  stateelyte.cs{ind};
        if ~isempty(val)
            state.(elyte).cs{ind} = val*ones(nc, 1);
        else
            state.(elyte).cs{ind} = NaN(nc, 1);
        end
    end

    %% setup initial volume fractions

    state.(an).volumeFraction = 0.3*ones(model.(an).G.getNumberOfCells(), 1);
    state.(ct).volumeFraction = 0.8*ones(model.(ct).G.getNumberOfCells(), 1);

    indsolid = model.(elyte).indsolidsp(1);
    state.(elyte).cs{indsolid} = zeros(nc, 1);
    state.(elyte).solidVolumeFraction = zeros(nc, 1);


    %% Setup intial quasi particle total concentration

    state = model.updateElectrolyteVolumeFraction(state);
    vf = state.(elyte).volumeFraction;

    for ind = 1 : nqp
        state.(elyte).qpepscs{ind} = vf.*state.(elyte).qpcs{ind};
    end

    %% initialize potentials

    state = model.initializeTemperature(state);

    state.(elyte).phi = zeros(model.Electrolyte.G.getNumberOfCells(), 1);

    state = model.evalVarName(state, {ctam, 'ENernst'});

    state.(ct).E = state.(ctam).ENernst(end);

    state.(an).E = 0;

    state = model.evalVarName(state, {anam, 'ENernst'});

    coupcells = model.couplingCellDict('Anode-Electrolyte');
    state.(elyte).phi(coupcells(:, 1)) = -state.(anam).ENernst(coupcells(:, 2));

    %% Initialize nucleation value

    state.(elyte).nucleation = zeros(model.(elyte).G.getNumberOfCells(), 1);
    state.(elyte).indicator = ones(model.(elyte).G.getNumberOfCells(), 1);

    initstate = state;

    %% We setup the schedule

    is_prep = model.include_precipitation;

    if is_prep

        dt = [];
        n = 30;
        dt = [dt; repmat(1e-6, n, 1).*1.8.^[1 : n]'];
        dT = dt(end);
        T = 0.03*hour; % roughly activation time
        dt = [dt; repmat(dT, floor(T/dT), 1)];
        T = 150*hour;
        n = 100;
        dt = [dt; repmat(T/n, n, 1)];
        T = 7*hour;
        n = 20;
        dt = [dt; repmat(T/n, n, 1)];

    else

        dt = [];
        n = 30;
        dt = [dt; repmat(1e-5, n, 1).*1.6.^[1 : n]'];
        T = 0.1*hour;
        n = 50;
        dt = [dt; repmat(T/n, n, 1)];

    end


    step = struct('val', dt, 'control', ones(numel(dt), 1));

    ct = 'Cathode';
    % stopFunc = @(model, state, state_prev) (state.(ct).E < 2.0);
    stopFunc = @(model, state, state_prev) (false);

    tup = 0.1;

    % given input current
    inputI = 50;
    % we use a ramp-up function
    srcfunc = @(time) rampupControl(time, tup, inputI);

    control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1);
    schedule = struct('control', control, 'step', step);

    %% We setup the nonlinear solver
    % Setup nonlinear solver
    nls = NonLinearSolver();
    % Change default maximum iteration number in nonlinear solver
    nls.maxIterations   = 15;
    nls.maxTimestepCuts = 10;
    nls.verbose         = true;
    % Change default behavior of nonlinear solver, in case of error
    nls.errorOnFailure = false;

    nls.timeStepSelector = IterationCountTimeStepSelector('targetIterationCount', 10, ...
                                                          'iterationOffSet'     , 1);

    % Change default tolerance for nonlinear solver
    model.nonlinearTolerance = 1e-4;

    model.verbose = true;

    output.model = model;
    output.schedule = schedule;

    % Run simulation

    dopacked = true;
    if dopacked
        simname = md5sum(input);
        dataFolder = opt.directory;
        problem = packSimulationProblem(initstate, model, schedule, dataFolder, 'Name', simname, 'NonLinearSolver', nls);
        output.problem = problem;
        output.dataDirectory = problem.OutputHandlers.states.dataDirectory;
        output.dataFolder = problem.OutputHandlers.states.dataFolder;
        output.input = input;
        inputfilename = fullfile(output.dataDirectory, output.dataFolder, 'input.mat');
        input.dataFolder = simname;
        input.directory = opt.directory;
        if ~opt.runSimulation
            output.input = input;
            try
                [~, states, report] = getPackedSimulatorOutput(problem);
                foundresults = true;
            catch
                foundresults = false;
            end
            if foundresults
                output.states = states;
                fprintf('\nResults of a previous simulation have been found and added to the output\n');
            end
            return
        end
        save(inputfilename, 'input');
        if opt.clearSimulation
            clearPackedSimulatorOutput(problem, 'Prompt', false);
        end
        simulatePackedProblem(problem);
        [~, states, report] = getPackedSimulatorOutput(problem);
        output.states = states;
    else
        [~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', ...
                                                        nls);
        output.states = states;
    end

end

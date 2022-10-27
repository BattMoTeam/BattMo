function obj = leastSquaresEI(model, states, schedule, varargin)

    opt = struct('ComputePartials', false, ...
                 'tStep'          , []   , ...
                 'state'          , []   , ...
                 'from_states'    , true , ...
                 'statesRef'      , []   , ...
                 'relTol'         , 1e-10);
    opt = merge_options(opt, varargin{:});

    statesRef = opt.statesRef;

    dts = schedule.step.val;

    if isempty(opt.tStep) %do all
        time = cumsum(dts);
        numSteps = numel(dts);
        tSteps = (1:numSteps)';
    else
        time = sum(dts(1:(opt.tStep)));
        numSteps = 1;
        dts = dts(opt.tStep);
        tSteps = opt.tStep;
    end

    obj = repmat({[]}, numSteps, 1);

    relTol = opt.relTol;

    for k = 1:numSteps
        t = time(k);
        dt = dts(k);

        % Find states for time t (given by schedule)
        state = findState(t, states, dt, relTol);
        stateRef = findState(t, statesRef, dt, relTol);

        if opt.ComputePartials
            if (opt.from_states)
                state = model.initStateAD(state);
            else
                state = opt.state;
            end
        end

        E = state.Control.E;
        I = state.Control.I;

        Eref = stateRef.Control.E;
        Iref = stateRef.Control.I;

        obj{k} = (E - Eref)^2 * dt + (I - Iref)^2 * dt;

    end

end


function state = findState(tstar, states, dt, relTol)

% Find state with time matching tstar

    times = cellfun(@(x) x.time, states);
    [tdiff, idx] = min(abs(times - tstar));

    assert(tdiff / dt < relTol);

    state = states{idx};

end
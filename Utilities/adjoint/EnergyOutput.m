function obj = EnergyOutput(model, states, schedule, varargin)

    opt     = struct('Price'          , 1.0  , ...
                     'ComputePartials', false, ...
                     'tStep'          , []   , ...
                     'state'          , []   , ...
                     'from_states'    , true);
    opt     = merge_options(opt, varargin{:});
    
    dts   = schedule.step.val;

    tSteps = opt.tStep;
    
    if isempty(tSteps) %do all
        time = 0;
        numSteps = numel(dts);
        tSteps = (1:numSteps)';
    else
        time = sum(dts(1:(opt.tStep-1)));
        numSteps = 1;
        dts = dts(opt.tStep);
    end

    obj = repmat({[]}, numSteps, 1);

    for step = 1:numSteps
        dt = dts(step);
        state = states{tSteps(step)};
        if opt.ComputePartials
            if (opt.from_states) 
                state = model.initStateAD(state);
            else
                state = opt.state;
            end
            E = state.Control.E;
            I = state.Control.I;
        else
            E = state.Control.E;
            I = state.Control.I; 
        end
        obj{step} = opt.Price*I*E*dt;
    end
    
end
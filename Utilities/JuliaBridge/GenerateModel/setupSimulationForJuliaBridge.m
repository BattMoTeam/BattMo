function output = setupSimulationForJuliaBridge(jsonstruct, varargin)

    opt = struct('runSimulation', false);
    opt = merge_options(opt, varargin{:});
    
    output = runBatteryJson(jsonstruct, 'runSimulation', opt.runSimulation);

    model     = output.model;
    schedule  = output.schedule;

    if opt.runSimulation

        output.use_state_ref = true;
        
        states   = output.states;

        %% added all the extra variables on state

        for istate = 1 : numel(states)
            states{istate} = model.addVariables(states{istate});
        end

        output.states = states;
    else

        output.use_state_ref = false;        
        
    end

    model = convertModelForJuliaBridge(model);

    model    = class2data(model);
    schedule = class2data(schedule);

    output.model    = model;
    output.schedule = schedule;
    
end

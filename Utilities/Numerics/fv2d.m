classdef fv2d
% structure to map state into variable y, yp that are needed for ode15i
% This is going to disappear when we switch to own Newton solver
    
    properties

        varsizes
        slots
        
        % Time discretization properties
        % should be moved elsewhere
        
        ti          % Initial time, s
        tf          % Final time, s
        dt          % Time step, s
        tUp         % Ramp up time, s
        tSpan       % Time span vector, s
    
    end
    
    methods
        
        function fv = fv2d(model, state)
            
            varsizes = [];
            varsizes(end + 1) = size(state.elyte.cs{1}, 1);
            varsizes(end + 1) = size(state.elyte.phi, 1);
            varsizes(end + 1) = size(state.ne.am.Li, 1);
            varsizes(end + 1) = size(state.ne.am.phi, 1);
            varsizes(end + 1) = size(state.pe.am.Li, 1);
            varsizes(end + 1) = size(state.pe.am.phi, 1);
            varsizes(end + 1) = size(state.ccne.phi, 1);
            varsizes(end + 1) = size(state.ccpe.phi, 1);
            varsizes(end + 1) = size(state.ccpe.E, 1);

            ind = 1;
            for i = 1 : numel(varsizes)
                slots{i} = ind : (ind + varsizes(i) - 1);
                ind = ind + varsizes(i);
            end
            
            fv.slots = slots;
            
            % Time discretization
            fv.ti = 0;
            fv.tf = 3600*24;
            fv.tf = 360*24;
            fv.dt = 50;
            fv.tUp = 0.1;
            fv.tSpan = (fv.ti : fv.dt : fv.tf);
            
        end

    end
    
end


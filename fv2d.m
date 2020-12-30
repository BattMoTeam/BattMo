classdef fv2d
% structure to map state into variable y, yp that are needed for ode15i
% This is going to disappear when we switch to own Newton solver
    
    properties

        varnames
        varsizes
        
        % State data slots
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
            
            varnames = model.getModelPrimaryVarNames();
            
            ind = 0;
            for i = 1 : numel(varnames)
                val = model.getProp(state, varnames{i});
                varsize =  size(val, 1);
                varsizes{i} = varsize;
                slots{i} = ind + (1 : varsize)';
                ind = ind + varsize;
            end
            
            fv.varnames = varnames;
            fv.varsizes = varsizes;
            fv.slots = slots;

            % Time discretization
            fv.ti = 0;
            fv.tf = 3600*24;
            fv.dt = 10;
            fv.tUp = 0.1;
            fv.tSpan = (fv.ti : fv.dt : fv.tf);
            
        end
        
        function n = varnum(fv)
            n = numel(fv.varnames);
        end
        
        function slot = getSlot(fv, varname)
        
            slots = fv.slots;
            varnames = fv.varnames;
            
            isfound = false;
            for ind = 1 : numel(varnames)
                if varname == varnames{ind}
                    isfound = true;
                    break
                end
            end
            assert(isfound, 'variable name not found');
            
            slot = slots{ind};
        end
        
    end
    
end


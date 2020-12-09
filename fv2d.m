classdef fv2d
    
    properties

        varnames
        varsizes
        varnum
        
        % Time discretization properties
        ti          % Initial time, s
        tf          % Final time, s
        dt          % Time step, s
        tUp         % Ramp up time, s
        tSpan       % Time span vector, s
        
        % State properties
        y
        y0
        yp
        yp0
        
        % State data slots
        slots
        
    end
    
    methods
        function obj = fv2d(model)
            
            compnames = model.compnames
            
            varnames = {};
            slots = {};
            ind = 0;
            
            for icn = 1 : numel(compnames)
                compname = compnames{icn};
                compvarnames = model.(compname).varnames;
                compvarsizes = model.(compname).varsizes;
                compvarnum   = model.(compname).varnum;
                for icv = 1 : compvarnum
                    compvarname = compvarnames{icv};
                    compvarsize = compvarsizes{icv};
                    varnames{end + 1} = fprintf('%s-%s', compname, compvarname);
                    varsizes{end + 1} = compvarsize;
                    slots{end + 1} = ind + (1 : compvarsize)';
                    ind = ind + compvarsize;
                end                
            end
            
        end
    end
    
    function n = varnum(obj)
        n = numel(obj.varnames);
    end
    
    function slot = getSlot(obj, varname)
        
        slots = obj.slots;
        varnames = obj.varnames;
        
        [isok, ind] = ismember(varname, varnames);
        assert(isok, 'variable name not found');
        
        slot = slots{ind};
    end
    
end


classdef fv2d
    
    properties

        varnames
        varsizes
        
        % State data slots
        slots
        
    end
    
    methods
        
        function obj = fv2d(model, state)
            
            varnames = model.getPrimaryVarNames();
            
            ind = 0;
            for i = 1 : numel(varnames)
                val = model.getProp(state, varnames{i});
                varsize =  size(val, 1);
                varsizes{i} = varsize;
                slots{i} = ind + (1 : varsize)';
                ind = ind + varsize;
            end
            
            obj.varnames = varnames;
            obj.varsizes = varsizes;
            obj.slots = slots;

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
    
end


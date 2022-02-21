classdef PropFunction
    
    properties
        
        varname
        
        inputvarnames
        modelnamespace 
        
        fn % function handler

    end
    
    methods
        
        function propfunction = PropFunction(varname, fn, inputvarnames, modelnamespace)
            
            propfunction.varname = varname;
            propfunction.fn = fn;
            propfunction.inputvarnames = inputvarnames;
            propfunction.modelnamespace = modelnamespace;
            
        end
        
        
        function propfunctions = resolveIndex(propfunction)
        % resolve the index for varname (when propfunction.varname is a cell, duplicate the propfunction for each cell entry)
            varname = propfunction.varname;
            varnames = varname.resolveIndex();
            propfunctions = {};
            for ind = 1 : numel(varnames)
                propfunctions{ind} = propfunction;
                propfunctions{ind}.varname = varnames{ind};
            end
            
        end

    end
    
    
end

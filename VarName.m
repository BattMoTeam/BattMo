classdef VarName
    
    properties
        namespace
        name
    end
    
    methods
        function varname = VarName(namespace, varname)
            varname = {};
            name = '';
        end
        
        function name = fullname(varname)
            name = join({varname.namespace{:}, varname.name}, '_');
        end
    end
    
    methods (statics)
        function name = joinvarnames(names)
            name = join({varname.namespace{:}, varname.name}, '_')
        end
    end        
    
end

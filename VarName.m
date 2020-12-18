classdef VarName
    
    properties
        namespace
        name
    end
    
    methods
        function varname = VarName(namespace, name)
            varname.namespace = namespace;
            varname.name = name;
        end
        
        function name = fullname(varname)
            name = join({varname.namespace{:}, varname.name}, '_');
            name = name{1};
        end
    end
    
    methods (Static)
        function name = joinvarnames(names)
            name = join({varname.namespace{:}, varname.name}, '_')
        end
    end        
    
end

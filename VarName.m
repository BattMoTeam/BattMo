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
        
        function name = getfieldname(varname)
            name = join({varname.namespace{:}, varname.name}, '_');
            name = name{1};
        end

        function isnequal = ne(varname, varname1)
            isequal  = varname.eq(varname1);
            if isequal
                isnequal = false;
            else
                isnequal = true;
            end
        end
        
        function isequal = eq(varname, varname1)
            
            namespace = varname.namespace;
            name = varname.name;
            namespace1 = varname1.namespace;
            name1 = varname1.name;
            
            isequal = true;
            
            if ~strcmp(name, name1)
                isequal = false;
                return
            end
            
            if numel(namespace) ~= numel(namespace1)
                isequal = false;
                return
            end
            
            for i = 1 : numel(namespace)
                if ~strcmp(namespace{i}, namespace{i})
                    isequal = false;
                    return
                end
            end
        end
        
    end
    
    methods (Static)
        function name = joinvarnames(names)
            name = join({varname.namespace{:}, varname.name}, '_')
        end
    end        
    
end

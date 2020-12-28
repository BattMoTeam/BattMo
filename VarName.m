classdef VarName
    
    properties
        namespace
        name
        index
    end
    
    methods
        function varname = VarName(namespace, name, index)
            varname.namespace = namespace;
            varname.name = name;
            if nargin > 2
                varname.index = index;
            else
                varname.index = ':';
            end
        end
        
        function name = getfieldname(varname)
            name = join({varname.namespace{:}, varname.name}, '_');
            name = name{1};
            if ~strcmp(varname.index, ':')
                name = sprintf('%s_%d', name, varname.index);
            end
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
            name      = varname.name;
            index     = varname.index;
            
            namespace1 = varname1.namespace;
            name1      = varname1.name;
            index1     = varname1.index;
            
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
            
            if not(ischar(index1) & ischar(index2)) | not(isnumeric(index1) & isisnumeric(index2)) 
                isequal = false;
                return
            end
            
            if (ischar(index1) & ischar(index2)) && ~strcmp(index1, index2)
                isequal = false;
                return
            end
            
            if (isnumeric(index1) & isnumeric(index2)) && (index1 == index2)
                isequal = false;
                return
            end
            
        end
        
    end
    
    methods (Static)
        function name = joinvarnames(names)
            name = join({varname.namespace{:}, varname.name}, '_')
        end
    end        
    
end

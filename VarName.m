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
        
        function isequal = eq(varname1, varname2)
            
            namespace1 = varname1.namespace;
            name1      = varname1.name;
            index1     = varname1.index;
            
            namespace2 = varname2.namespace;
            name2      = varname2.name;
            index2     = varname2.index;
            
            % for the moment, we do not handle the special fields '..', '.'
            namespace = horzcat(namespace1, namespace2);
            assert(all(~strcmp('..', namespace)), 'we do not handle double dots');
            assert(all(~strcmp('.', namespace)), 'we do not handle dot');
                
            isequal = true;
            
            if ~strcmp(name1, name2)
                isequal = false;
                return
            end

            if numel(namespace1) ~= numel(namespace2)
                isequal = false;
                return
            end
            
            for i = 1 : numel(namespace1)
                if ~strcmp(namespace1{i}, namespace2{i})
                    isequal = false;
                    return
                end
            end
            
            if (~ischar(index1) & ischar(index2)) | (ischar(index1) & ~ischar(index2))
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

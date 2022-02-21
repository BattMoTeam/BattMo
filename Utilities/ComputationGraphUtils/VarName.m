classdef VarName
    
    properties
        
        namespace
        name
        index
        dim
        
        
    end
    
    methods
        function varname = VarName(namespace, name, dim, index)
            varname.namespace = namespace;
            varname.name = name;
            
            if nargin > 2
                varname.dim = dim;
            else
                varname.dim = 1;
            end
            
            if nargin > 3
                varname.index = index;
            else
                varname.index = ':';
            end
            
        end
        
        function name = getfieldname(varname)
            name = join({varname.namespace{:}, varname.name}, '_');
            name = name{1};
        end

        function name = getIndexedFieldname(varname)
        % Generate fieldname with the index in bracket (this entry is used in the dictionaries for AdminModel instance)
            name = varname.getfieldname();
            index = varname.index;
            dim = varname.dim;
            if dim > 1
                assert(isnumeric(index) && numel(index) == 1)
                name = sprintf('%s[%d]', name, index);
            end
        end
        
        
        function varnames = resolveIndex(varname)
        % If varname.dim > 1, produce a cell array of VarName, one entry by element in varname.index
            dim = varname.dim;
            if dim > 1
                namespace = varname.namespace;
                name = varname.name;
                index = varname.index;
                
                if ischar(index) && strcmp(':', index)
                    index = [1 : dim];
                end
                
                assert(isnumeric(index), 'index should be numeric');
                varnames = {};
                for ind = 1 : numel(index)
                    varnames{ind} = VarName(namespace, name, dim, index(ind));
                end
            else
                varname.index = 1;
                varnames = {varname};
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

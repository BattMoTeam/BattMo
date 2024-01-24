classdef VarName

    properties (Constant)
        
        strlim = '.' % string delimiter when creating variable name for visualization from namespace and name
        
    end
    
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
        
        function name = getFieldname(varname)
            name = strjoin({varname.namespace{:}, varname.name}, varname.strlim);
        end

        function name = getIndexedFieldname(varname)
        % Generate fieldname with the index in bracket
            name = varname.getFieldname();
            index = varname.index;
            dim = varname.dim;
            if dim > 1
                assert(isnumeric(index) && numel(index) == 1)
                name = sprintf('%s{%d}', name, index);
            end
        end

        function cellname = getPropName(varname)
        % Generate name that is compatible with getProp function in class BaseModel
            cellname = varname.namespace;
            cellname{end + 1} = varname.name;
            index = varname.index;
            dim = varname.dim;
            isok = (isnumeric(index) && numel(index) == 1);
            isok = isok | (ischar(index) && dim == 1);
            assert(isok, 'The variable has multiple index and cannot be processed by getProp');
            if dim > 1
                cellname{end + 1} = index;
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

        function [isequal, compIndices] = compareVarName(varname1, varname2)

            compIndices = [];
            
            namespace1 = varname1.namespace;
            name1      = varname1.name;
            
            namespace2 = varname2.namespace;
            name2      = varname2.name;
            
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
            
            dim1     = varname1.dim;
            dim2     = varname2.dim;

            index1     = varname1.index;
            index2     = varname2.index;

            if (dim1 ~= dim2)
                isequal = false;
                error('There exists two declared variable with the same name and namespace but different dimensions.')
                return
            elseif dim1 == 1
                isequal = true;
                return
            else
                dim = dim1;
            end

            % Same namespace, name and dimension and dimension larger than 1, we process the indices
            if (ischar(index1))
                index1 = (1 : dim)';
            end
            if (ischar(index2))
                index2 = (1 : dim)';
            end
            
            lindex1 = false(dim, 1);
            lindex1(index1) = true;
            lindex2 = false(dim, 1);
            lindex2(index2) = true;

            compIndices.InterInd  = find(lindex1 & lindex2);
            compIndices.OuterInd1 = find(lindex1 & ~lindex2);
            compIndices.OuterInd2 = find(~lindex1 & lindex2);    
            compIndices.index1    = find(index1);
            compIndices.index2    = find(index2);
            
            if all(lindex1 == lindex2)
                isequal = true;
            else
                isequal = false;
            end
            
        end
        
        function isequal = eq(varname1, varname2)
            
            namespace1 = varname1.namespace;
            name1      = varname1.name;
            
            namespace2 = varname2.namespace;
            name2      = varname2.name;
            
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
            
            dim1     = varname1.dim;
            dim2     = varname2.dim;

            index1     = varname1.index;
            index2     = varname2.index;

            if (dim1 ~= dim2)
                isequal = false;
                return
            elseif dim1 == 1
                isequal = true;
                return
            else
                dim = dim1;
            end

            % Same dimension and dimension larger than 1, we have to check indexes
            if (ischar(index1))
                
                index = false(dim, 1);

                index1 = index;
                index1(varname1.index) = true;
                index2 = index;
                index2(varname1.index) = true;
                
                if (~ischar(index1) & ischar(index2)) | (ischar(index1) & ~ischar(index2))
                    isequal = false;
                    return
                end
                
                if (ischar(index1) & ischar(index2)) && ~strcmp(index1, index2)
                    isequal = false;
                    return
                end
                
                if (isnumeric(index1) & isnumeric(index2))

                    index = false()
                    
                    isequal = false;
                    return
                end
                
            end
            
        end

        
    end

end




%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}

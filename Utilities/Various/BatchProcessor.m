classdef BatchProcessor
    properties
        paramnames
    end
    
    methods
        
        function bp = BatchProcessor()
            bp.paramnames = {};
        end
        
        function [bp, simlist] = addElement(bp, simlist, elt)
            simlist{end + 1} = elt;
            fdnames = fieldnames(elt);
            for ind = 1 : numel(fdnames)
                fdname = fdnames{ind};
                bp = bp.addParameterName(fdname);
            end
        end

        function [bp, simlist] = mergeSimLists(bp, simlist, simlist_to_merge)
            for ind = 1 : numel(simlist_to_merge)
                elt = simlist_to_merge{ind};
                [bp, simlist] = bp.addElement(simlist, elt);
            end
        end
        
        
        function bp = addParameterName(bp, paramname)
            if ~ismember(paramname, bp.paramnames)
                bp.paramnames{end + 1} = paramname;
            end        
        end
        
        function ind = getParameterIndex(bp, paramname)
            [~, ind] = ismember(paramname, bp.paramnames);
        end
        
        function [bp, simlist] = modifyElement(bp, simlist, ind, paramname, paramvalue)
            error('not robust')
            simlist{ind}.paramname = paramvalue;
            bp = bp.addParameterName(paramname);
        end        
    
        function [vals, type, dim] = getParameterValue(bp, simlist, paramname)
            bp.assertParam(paramname);
            vals = {};
            types = {};
            dims = {};
            for ind = 1 : numel(simlist)
                s = simlist{ind};
                val  = {};
                type = {};
                dim  = [];
                if isfield(s, paramname)
                    val = s.(paramname);
                    if isnumeric(val)
                        dim = size(val);
                        if all(dim == [1, 1])
                            type = 'scalar';
                        elseif isempty(val)
                            type = 'undef';
                        else
                            type = 'vector';
                        end
                    elseif isa(val, 'char')
                        type = 'char';
                    elseif isa(val, 'logical') & numel(val) == 1
                        type = 'boolean';
                    else
                        error('type not recognized');
                    end
                else
                    type = 'undef';
                end
                vals{end + 1}  = val;
                types{end + 1} = type;
                dims{end + 1}  = dim;
                
            end

            ok = false;

            if all(ismember(types, {'char', 'undef'}))
                ok = true;
                type = 'char';
                dim = [];
            end

            if all(ismember(types, {'boolean', 'undef'}))
                ok = true;
                type = 'boolean';
                dim = [];
            end
            
            if ~isempty(vals) && (isnumeric(vals{1}) | isempty(vals{1}))
                if all(ismember(types, {'scalar', 'undef'}))
                    ok = true;
                    type = 'scalar';
                    dim = 1;
                elseif all(ismember(types, {'vector', 'undef'}))
                    type = 'vector';
                    isdef = ismember(types, 'vector');
                    dims = vertcat(dims{isdef});
                    if all(dims(:, 2) == 1)
                        a = unique(dims(:, 1));
                        if numel(a) == 1
                            ok = true;
                            dim = [dims(1, 1), 1];
                        end
                    elseif all(dims(:, 1) == 1)
                        a = unique(dims(:, 2));
                        if numel(a) == 1
                            ok = true;
                            dim = [1, dims(1, 2)];
                        end
                    end
                end
            end

            if ~ok
                error('values are not consistents');
            end
            
            
        end        

        function T = setupTable(bp, simlist)

            paramnames = bp.paramnames;

            nlist = numel(simlist);
            nparams = numel(paramnames);

            Tc = cell(nlist, nparams);

            for ic  = 1 : nparams
                paramname = paramnames{ic};
                [vals, type, dim] = getParameterValue(bp, simlist, paramname);
                for ir = 1 : numel(simlist)
                    if ismember(type, {'char', 'scalar', 'boolean'})
                        if isempty(vals{ir})
                            Tc{ir, ic} = 'undefined';
                        else
                            Tc{ir, ic} = vals{ir};
                        end
                    elseif strcmp(type, 'vector')
                        if isempty(vals{ir})
                            Tc{ir, ic} = 'undefined';
                        elseif dim(1) == 1
                            Tc{ir, ic} = sprintf('[%s]', strjoin(cellstr(num2str(vals(ir, :)')), ', '));
                        elseif dim(2) == 1
                            Tc{ir, ic} = sprintf('[%s]', strjoin(cellstr(num2str(vals(:, ir))), '; '));
                        else
                            error('problem with dimension.');
                        end
                    end
                end
            end

            T = cell2table(Tc, 'VariableNames', paramnames);
            
        end
        
        function rangevalues = getParameterRanges(bp, simlist, varargin)
            
            if nargin > 2
                paramnames = varargin;
                for iparam = 1 : numel(paramnames)
                    bp.assertParam(paramnames{iparam});
                end                    
            else
                paramnames = bp.paramnames;
            end
            
            rangevalues = [];
            for ind = 1 : numel(paramnames)
                paramname = paramnames{ind};
                [vals, type, dim] = getParameterValue(bp, simlist, paramname);
                ind = cellfun(@(val) ~isempty(val), vals);
                vals = vals(ind);
                if strcmp(type, 'vector')
                    if dim(1) == 1
                        vals = unique(vals, 'rows');
                    elseif dim(2) == 1
                        vals = unique(vals', 'rows');
                        vals = vals';
                    else
                        error('Matrix entry are not coverged');
                    end
                elseif ismember(type, {'scalar', 'boolean'})
                    vals = vertcat(vals{:});
                    vals = unique(vals);
                else
                    vals = unique(vals);                    
                end
                rangevalues.(paramname) = vals;
            end
            
        end
        
        function assertParam(bp, paramname)
            assert(ismember(paramname, bp.paramnames), sprintf('parameter %s not registered', paramname));
        end
        
        function printSimList(bp, simlist, varargin)
            if nargin > 2
                paramnames = varargin;
                for iparam = 1 : numel(paramnames)
                    bp.assertParam(paramnames{iparam});
                end                    
            else
                paramnames = bp.paramnames;
            end
            T = bp.setupTable(simlist);
            T = T(:, paramnames);
            n = size(T, 1);
            take = true(size(T, 2), 1);
            for itake = 1 : numel(take)
                val = T(:, itake);
                if all(strcmp(table2cell(T(:, itake)), repmat('undefined', n ,1)))
                    take(itake) = false;
                end
            end
            T = T(:, take)
                
        end        
    
        function sortedsimlist = sortSimList(bp, simlist, varargin)
            paramname = varargin{1};
            rest = varargin(2 : end);
            if strcmp(paramname, 'inverseOrder')
                sortedsimlist = simlist(end : -1 : 1);
            else
                [vals, type] = getParameterValue(bp, simlist, paramname);
                eind = cellfun(@(val) isempty(val), vals);
                simlist = horzcat(simlist(~eind), simlist(eind));
                vals = vals(~eind);
                if ismember(type, {'scalar', 'boolean'})
                    vals = vertcat(vals{:});
                elseif strcmp(type, 'vector')
                    error('vectors cannot be ordered');
                end
                [~, ind] = sort(vals);
                ne = nnz(~eind);
                sortedsimlist = simlist;
                sortedsimlist(1 : ne) = simlist(ind);
            end
            if ~isempty(rest)
                sortedsimlist = sortSimList(bp, sortedsimlist, rest{:});
            end
        end
        
        function [bp, simlist] = filterMergeSimLists(bp, simlist, simlist_to_merge, varargin)
            
            [bp, simlist_to_merge] = bp.mergeSimLists({}, simlist_to_merge);
            simlist_to_merge = bp.filterSimList(simlist_to_merge, varargin{:});
            [bp, simlist] = bp.mergeSimLists(simlist, simlist_to_merge);
            
        end
        
        
        function filteredsimlist = filterSimList(bp, simlist, varargin)
            assert(mod(numel(varargin), 2) == 0, 'wrong number of argument')
            filteredsimlist = {};
            paramname = varargin{1};
            filter = varargin{2};
            rest = varargin(3 : end);
            for isim = 1 : numel(simlist)
                s = simlist{isim};
                if isfield(s, paramname)
                    paramval = s.(paramname);
                    take = false;
                    if isa(filter, 'numeric') || isa(filter, 'logical')
                        if paramval == filter
                            take = true;
                        end
                    elseif isa(filter, 'char')
                        if strcmp(paramval, filter)
                            take = true;
                        end
                    elseif isa(filter, 'function_handle')
                        if filter(paramval)
                            take = true;
                        end
                    else
                        error('filter type not recognized');
                    end
                    if take == true
                        filteredsimlist{end + 1} = s;
                    end
                end
                if ~isfield(s, paramname) & (isempty(filter) | strcmp(filter, 'undefined'))
                    filteredsimlist{end + 1} = s;
                end
            end
            if ~isempty(rest)
                filteredsimlist = bp.filterSimList(filteredsimlist, rest{:});
            end
        end        
    
    end
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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

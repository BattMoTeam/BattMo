classdef BatchProcessor < Selector

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

    properties
        
        paramnames

    end
    
    methods
        
        function bp = BatchProcessor(varargin)

            bp = bp@Selector();
            
            bp.paramnames = {};
            if nargin > 0
                simlist = varargin{1};
                for isim = 1 : numel(simlist)
                    bp = bp.registerElementParamnames(simlist{isim});
                end
            end
        end

        function bp = registerElementParamnames(bp, elt)
            fdnames = fieldnames(elt);
            for ind = 1 : numel(fdnames)
                fdname = fdnames{ind};
                bp = bp.addParameterName(fdname);
            end            
        end
    
        function [bp, simlist] = addElement(bp, simlist, elt)
            simlist{end + 1} = elt;
            bp = bp.registerElementParamnames(elt);
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

        function [T, singlevalued] = setupTable(bp, simlist)

            paramnames = bp.paramnames;

            nlist = numel(simlist);
            nparams = numel(paramnames);

            Tc = cell(nlist, nparams);

            singlevalued = false(nparams, 1);
            
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
                            Tc{ir, ic} = sprintf('[%s]', strjoin(cellstr(num2str(vals{ir})), '; '));
                        else
                            error('problem with dimension.');
                        end
                    end
                end

                singlevalue = false;
                if nnz(cellfun(@(val) isempty(val), vals)) == numel(vals)
                    singlevalue = true;
                elseif nnz(cellfun(@(val) isempty(val), vals)) > 0
                    singlevalue = false;
                else
                    if ismember(type, {'scalar', 'boolean'}) && numel(unique([vals{:}])) == 1
                        singlevalue = true;
                    elseif ismember(type, {'char'}) && (numel(unique(vals)) == 1)
                        singlevalue = true;
                    end
                end
                singlevalued(ic) = singlevalue;

            end

            singlevalued = find(singlevalued);
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
        
        function print(bp, simlist, varargin)
            
            [T, singlevalued] = bp.setupTable(simlist);
            if nargin > 2 && strcmp(varargin{1}, 'all')
                paramnames = bp.paramnames;
            elseif (nargin > 2 && ismember('varying', varargin)) | (nargin == 2)
                vparamnames = bp.paramnames;
                ind = true(numel(vparamnames), 1);
                ind(singlevalued) = false;
                vparamnames = vparamnames(ind);
                if nargin > 2
                    [lia, loc] = ismember('varying', varargin);
                    ind = true(numel(varargin), 1);
                    ind(loc) = false;
                    aparamnames = varargin(ind);
                    for iparam = 1 : numel(aparamnames)
                        bp.assertParam(aparamnames{iparam});
                    end
                    nv = numel(vparamnames);
                    na = numel(aparamnames);
                    inda = [(1 : loc - 1)'; ((loc + nv) : (na + nv))'];
                    indv = (loc : loc + nv - 1)';
                    paramnames = cell(1, na + nv);
                    paramnames(inda) = aparamnames;
                    paramnames(indv) = vparamnames;
                else
                    paramnames = vparamnames;
                end
            else                
                paramnames = varargin;
                for iparam = 1 : numel(paramnames)
                    bp.assertParam(paramnames{iparam});
                end
            end

            if numel(paramnames) == 0
                if numel(simlist) > 1
                    % we print all the parameters instead of empty table
                    bp.print(simlist, 'all');
                else
                    simlist{1}
                end
            else
                T = T(:, paramnames)
            end
        end        
        
        function sortedsimlist = sort(bp, simlist, varargin)
            paramname = varargin{end};
            rest = varargin(1 : end - 1);

            direction = 'ascend';
            usefunc   = false;
            
            if iscell(paramname)
                func      = paramname{2};
                paramname = paramname{1};
                if isa(func, 'function_handle')
                    usefunc = true;
                elseif ischar(func) && ismember(func, {'ascend', 'descend'})
                    direction = func;
                else
                    error('filter format not recognized');
                end
                
            end
            [vals, type] = getParameterValue(bp, simlist, paramname);
            eind = cellfun(@(val) isempty(val), vals);
            simlist = horzcat(simlist(~eind), simlist(eind));
            vals = vals(~eind);
            if ismember(type, {'scalar', 'boolean'})
                vals = vertcat(vals{:});
            elseif strcmp(type, 'vector')
                error('vectors cannot be ordered');
            end
            
            if usefunc
                vals = func(vals);
            end

            [~, ~, ic] = unique(vals, 'sorted');
            inds = [ic, (1 : numel(vals))'];
            
            switch direction
              case 'ascend'
                [~, ind] = sortrows(inds);
              case 'descend'
                [~, ind] = sortrows(inds, [-1, 2]);
              otherwise
                error('direction not recognized');
            end
            ne = nnz(~eind);
            sortedsimlist = simlist;
            sortedsimlist(1 : ne) = simlist(ind);
            if ~isempty(rest)
                sortedsimlist = sortSimList(bp, sortedsimlist, rest{:});
            end
        end
        
        function [bp, simlist] = filterMergeSimLists(bp, simlist, simlist_to_merge, varargin)
            
            [bp, simlist_to_merge] = bp.mergeSimLists({}, simlist_to_merge);
            simlist_to_merge = bp.filterSimList(simlist_to_merge, varargin{:});
            [bp, simlist] = bp.mergeSimLists(simlist, simlist_to_merge);
            
        end

        function filteredsimlist = filter(bp, simlist, varargin)
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
                    if (isempty(filter) | strcmp(filter, 'undefined'))
                        if isempty(paramval)
                            take = true;
                        end
                    elseif isa(filter, 'numeric') || isa(filter, 'logical')
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


        %%%%%%%%%%%%
        %% Selector class overloaded function
        %%%%%%%%%%%%%%

        function str = selectSelectorToString(slt, selectSelector)
        % Returns the printed form of a 'select' selector

            assert(strcmp(selectSelector{1}, 'select'), 'this is not a select type selector');
            selector_type  = selectSelector{2}{1};
            selector_value = selectSelector{2}{2};
            
            if isa(selector_value, 'function_handle')
                selector_value = func2str(selector_value);
            end
            
            str = sprintf('%s : %s', selector_type, selector_value);
            
        end

        function printSelection(slt, givenset, selection)

            if ~strcmp(selection{1}, 'set')
                selection = slt.parseSelector(givenset, selection);
            end
            
            inds = selection{2};

            simlist = givenset(inds);

            printall = getStructField(slt.interactiveOptions, {'printSelection', 'all'}, true);

            if printall
                slt.print(simlist, 'all');
            else
                slt.print(simlist, 'all');
            end

        end

        function found = find(slt, givenset, selector)
            
            paramname = selector{1};
            filter    = selector{2};

            found = [];
            
            for ielt = 1 : numel(givenset)
                
                elt = givenset{ielt};

                paramnames = fieldnames(elt);

                indparams = regexpSelect(paramnames, paramname);
                
                for iindparams = 1 : numel(indparams)

                    indparam = indparams(iindparams);
                    
                    paramval = elt.(paramnames{indparam});

                    take = false;
                    
                    if (isempty(filter) | strcmp(filter, 'undefined'))
                        if isempty(paramval)
                            take = true;
                        end
                    elseif isa(filter, 'numeric') || isa(filter, 'logical')
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
                        found(end + 1) = ielt;
                    end
                end
                
            end

            found = unique(found);
            
        end
        

    end

    methods(Static)

        function [bp, simlist] = mergeSimList(simlist, simlist_to_merge)
            
            bp = BatchProcessor(simlist);
            for isim = 1 : numel(simlist_to_merge)
                elt = simlist_to_merge{isim};
                [bp, simlist] = bp.addElement(simlist, elt);
            end

        end
        
        function [bp, simlist] = mergeSimLists(simlists)

            simlist = simlists{1};
            simlists = simlists(2 : end);

            while numel(simlists) > 0

                simlist_to_merge = simlists{1};
                simlists = simlists(2 : end);
                
                [bp, simlist] = BatchProcessor.mergeSimList(simlist, simlist_to_merge);
                
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

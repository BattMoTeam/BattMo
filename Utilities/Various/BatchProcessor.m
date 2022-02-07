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
        
        function bp = addParameterName(bp, paramname)
            if ~ismember(paramname, bp.paramnames)
                bp.paramnames{end + 1} = paramname;
            end        
        end
        
        function ind = getParameterIndex(bp, paramname)
            [~, ind] = ismember(paramname, bp.paramnames);
        end
        
        function [bp, simlist] = modifyElement(bp, simlist, ind, paramname, paramvalue)
            error('not robbust')
            simlist{ind}.paramname = paramvalue;
            bp = bp.addParameterName(paramname);
        end        
    
        function vals = getParameterValue(bp, simlist, paramname)
            bp.assertParam(paramname);
            vals = {};
            for ind = 1 : numel(simlist)
                s = simlist{ind};
                if isfield(s, paramname)
                    vals{end + 1} = s.(paramname);
                else
                    vals{end + 1} = {};
                end
            end
            
            if isnumeric(vals{1})
                vals = cell2mat(vals);
            end
            
        end        

        function T = setupTable(bp, simlist)

            paramnames = bp.paramnames;

            nlist = numel(simlist);
            nparams = numel(paramnames);

            Tc = cell(nlist, nparams);

            for ir  = 1 : nlist
                s = simlist{ir};
                for ic  = 1 : nparams
                    paramname = paramnames{ic};
                    if isfield(s, paramname) && ~isempty(s.(paramname))
                        Tc{ir, ic} = s.(paramname); 
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
                vals = getParameterValue(bp, simlist, paramname);
                vals = unique(vals);
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
            T(:, paramnames)
        end        
    
        function sortedsimlist = sortSimList(bp, simlist, varargin)
            paramname = varargin{1};
            rest = varargin(2 : end);
            vals = getParameterValue(bp, simlist, paramname);
            [~, ind] = sort(vals);
            sortedsimlist = simlist(ind);                
            if ~isempty(rest)
                sortedsimlist = sortSimList(bp, sortedsimlist, rest{:});
            end
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
                    if isa(filter, 'numeric')
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
            end
            if ~isempty(rest)
                filteredsimlist = bp.filterSimList(filteredsimlist, rest{:});
            end
        end        
    
    end
    
end



%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}

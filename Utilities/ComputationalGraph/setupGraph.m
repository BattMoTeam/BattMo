function [A, staticprops, resolvedVarNameList, nodenames] = setupGraph2(model, varargin)
    
    opt = struct('resolveIndex', true);
    opt = merge_options(opt, varargin{:});
    
    ss = []; % source nodes
    ts = []; % target nodes
    fs = []; % edge function name
    ps = []; % property function index in property list
    
    varnames  = model.varNameList;
    propfuncs = model.propertyFunctionList;

    nodenames = {};
    resolvedVarNameList = {};
    for ind = 1 : numel(varnames)
        varname = varnames{ind};
        if opt.resolveIndex
            varname_s = varname.resolveIndex();
            for ind = 1 : numel(varname_s)
                fullname = varname_s{ind}.getIndexedFieldname();
                nodenames{end + 1} = fullname;
            end
            resolvedVarNameList = horzcat(resolvedVarNameList, varname_s);
        else
            fullname = varname.getFieldname;
            nodenames{end + 1} = fullname;            
        end
    end
    
    staticprops = {};
    
    for ipropfunc = 1 : numel(propfuncs)

        propfunction = propfuncs{ipropfunc};
        
        varname = propfunction.varname;
        m       = propfunction.modelnamespace;

        inputvarnames = propfunction.inputvarnames;
        
        fullinputvarnames = {};
        
        for i = 1 : numel(inputvarnames)
            inputvarname = inputvarnames{i};
            if opt.resolveIndex
                inputvarnames_s =  inputvarname.resolveIndex();
                for j = 1 : numel(inputvarnames_s)
                    fullinputvarname = inputvarnames_s{j}.getIndexedFieldname();
                    fullinputvarnames = horzcat(fullinputvarnames, {fullinputvarname});
                end
            else
                inputvarname = inputvarname.getFieldname;
                fullinputvarnames = horzcat(fullinputvarnames, {inputvarname});
            end
        end

        fullvarnames = {};
        
        if opt.resolveIndex
            varnames = varname.resolveIndex();
            for ivarname = 1 : numel(varnames)
                fullvarnames{ivarname} = varnames{ivarname}.getIndexedFieldname();
            end
        else
            fullvarnames = {varname.getFieldname};
        end
        
        nv = numel(fullvarnames);
        ni = numel(fullinputvarnames);

        % check that variables have been declared
        allfullvarnames = horzcat(fullvarnames, fullinputvarnames);
        ind = ismember(allfullvarnames, nodenames);
        if any(~ind)
            undefinedvarnames = allfullvarnames(~ind);
            fprintf('The following variables are used in function declaration (function name %s) but have not be declared for themselves before:\n', f);
            for ivar = 1 : numel(undefinedvarnames)
                fprintf('%s\n', undefinedvarnames{ivar});
            end
        end
        
        if ni == 0
            
            staticprop.propind = ipropfunc;
            for inv = 1 : nv
                nodename = fullvarnames{inv};
                staticprop.nodename = nodename;
                staticprops{end + 1} = staticprop;
            end
            
        else

            [~, indv] = ismember(fullvarnames, nodenames);
            [~, indi] = ismember(fullinputvarnames, nodenames);

            clear sourcetbl
            sourcetbl.source = indi';
            sourcetbl = IndexArray(sourcetbl);
            clear targettbl
            targettbl.target = indv';
            targettbl = IndexArray(targettbl);
            
            sourcetargettbl = crossIndexArray(sourcetbl, targettbl, {});
            sourcetargettbl = sortIndexArray(sourcetargettbl, {'source', 'target'});
            
            s = sourcetargettbl.get('source');
            t = sourcetargettbl.get('target');

            p = repmat(ipropfunc, nv*ni, 1);
            
            ss = vertcat(ss, s);
            ts = vertcat(ts, t);
            ps = vertcat(ps, p);
            
        end
        
    end

    n = numel(nodenames);
    
    A = sparse(ss, ts, ps, n, n);
    
end




%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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

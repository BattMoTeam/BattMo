function [g, staticprops, resVarNameList] = setupGraph(model, varargin)
    
    opt = struct('resolveIndex', true);
    opt = merge_options(opt, varargin{:});
    
    g = digraph();
    
    ss = {}; % source nodes
    ts = {}; % target nodes
    fs = {}; % edge function name
    ms = {}; % edge model name    
    ps = {}; % property function index in property list
    
    varnames  = model.varNameList;
    propfuncs = model.propertyFunctionList;

    resVarNameList = {};
    for ind = 1 : numel(varnames)
        varname = varnames{ind};
        if opt.resolveIndex
            varname_s = varname.resolveIndex();
            for ind = 1 : numel(varname_s)
                fullname = varname_s{ind}.getIndexedFieldname();
                g = addnode(g, fullname);
            end
            resVarNameList = horzcat(resVarNameList, varname_s);
        else
            fullname = varname.getFieldname;
            g = addnode(g, fullname);
        end
    end
    
    nodenames = g.Nodes.Variables;
    
    staticprops = {};
    
    for ipropfunc = 1 : numel(propfuncs)

        propfunction = propfuncs{ipropfunc};
        varname = propfunction.varname;
        m = propfunction.modelnamespace;

        f = func2str(propfunction.fn);
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
            indv = rldecode((1 : nv)', ni*ones(nv, 1))';
            indi = repmat((1 : ni), 1, nv);
            
            s = fullinputvarnames(indi);
            t = fullvarnames(indv);
            f = repmat({f}, 1, nv*ni);
            m = repmat({m}, 1, nv*ni);
            p = repmat({ipropfunc}, 1, nv*ni);
            
            ss = horzcat(ss, s);
            ts = horzcat(ts, t);
            fs = horzcat(fs, f);
            ps = horzcat(ps, p);
            ms = horzcat(ms, m);
        end
        
    end
    
    g = addedge(g, ss, ts, [ps{:}]);
    
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

function varname = resolveAlias(model, name, index)
% the function is recursive. The index to varname will be the first index that is detected (if different from ':')
    if nargin < 3
        index = ':';
    end
    [isalias, varname] = model.aliasLookup(name);
    if isalias
        model = model.getAssocModel(varname.namespace);
        if ischar(index) 
            index = varname.index;
        end
        varname = resolveAlias(model, varname.name, index);
    else
        varname = VarName(model.namespace, name);
        varname.index = index;
    end
end

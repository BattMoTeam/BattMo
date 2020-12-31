function [varfieldnames, varindices] = updateVarDict(varfieldnames, varindices, model)
    names = model.names;
    for ind = 1 : numel(names)
        name = names{ind};
        varname = VarName(model.namespace, name);
        fieldname = varname.getfieldname();
        varname = resolveAlias(model, name);
        reffieldname = varname.getfieldname();
        varfieldnames(fieldname) = reffieldname;
        varindices(fieldname) = varname.index;
    end
    aliases = model.aliases;    
    for ind = 1 : numel(aliases)
        name = aliases{ind}{1};
        varname = VarName(model.namespace, name);
        fieldname = varname.getfieldname();
        varname = resolveAlias(model, name);
        reffieldname = varname.getfieldname();
        varfieldnames(fieldname) = reffieldname;
        varindices(fieldname) = varname.index;
    end
    
    if isa(model, 'CompositeModel')
        for ind = 1 : numel(model.SubModels)
            [varfieldnames, varindices] = updateVarDict(varfieldnames, varindices, model.SubModels{ind});
        end
    end
end

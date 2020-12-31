function [propfunctdict, propmodeldict] = updateFunctionDict(propfunctdict, propmodeldict, model)
    
    propfunctions = model.propfunctions;
    
    for ind = 1 : numel(propfunctions)
        
        propfunction = propfunctions{ind};
        name = propfunction.name;
        fn = propfunction.fn;
        fnmodel = model.getAssocModel(propfunction.namespace);
        
        varname = VarName(model.namespace, name);
        fieldname = varname.getfieldname;
        
        propfunctdict(fieldname) = fn;
        
        propmodeldict(fieldname) = fnmodel.getModelFullName();
        
    end
    
    if isa(model, 'CompositeModel')
        for ind = 1 : numel(model.SubModels)
            [propfunctdict, propmodeldict] = updateFunctionDict(propfunctdict, propmodeldict, model.SubModels{ind});
        end
    end
    
end
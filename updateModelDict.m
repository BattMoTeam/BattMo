function models = updateModelDict(models, model)
    
    if model.hasparent
        name = model.getModelFullName;
        models(name) = model;
    else
        name = 'root';
        models(name) = model;
    end
    
    if isa(model, 'CompositeModel')
        for ind = 1 : numel(model.SubModels)
            models = updateModelDict(models, model.SubModels{ind});
        end
    end
    
end

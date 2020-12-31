function models = updateModelDict(models, model)
    
    name = model.getModelFullName;
    models(name) = model;
    
    if isa(model, 'CompositeModel')
        for ind = 1 : numel(model.SubModels)
            models = updateModelDict(models, model.SubModels{ind});
        end
    end
    
end

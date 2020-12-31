function modeldict = updateModelDict(modeldict, model)
    
    name = model.getModelFullName;
    modeldict(name) = model;
    
    if isa(model, 'CompositeModel')
        for ind = 1 : numel(model.SubModels)
            modeldict = updateModelDict(modeldict, model.SubModels{ind});
        end
    end
    
end

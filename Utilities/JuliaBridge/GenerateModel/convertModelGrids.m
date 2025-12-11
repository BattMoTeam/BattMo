function model = convertModelGrids(model)
% Convert the grids from sub-grid representation to a concrete grid with the data needed for a julia simulation.
    
    if isprop(model, 'G') && ~isempty(model.G)
        G = model.G.mrstFormat;
        op.haltTrans = 
        G.operators = op;
        model.G = G;
    end

    submodelnames = model.getSubModelNames();

    if ~isempty(submodelnames)

        for isub = 1 : numel(submodelnames)

            submodelname = submodelnames{isub};

            model.(submodelname) = convertModelGrids(model.(submodelname));

        end

    end
    
end


function model = convertModelForJuliaBridge(model)

    % convert the grids from sub-grids, adding the transmissibilities and other quantities needed in the julia
    % simulation.
    model = convertModelGrids(model);

    model = convertModelForJuliaBridge_(model);
    
end


function model = convertModelForJuliaBridge_(model)

    props = properties(model);

    for iprop = 1 : numel(props)
        prop = props{iprop};
        if isa(model.(prop), 'function_handle')
            model.(prop) = func2str(model.(prop));
        end
    end

    submodelnames = model.getSubModelNames();

    if ~isempty(submodelnames)

        for isub = 1 : numel(submodelnames)

            submodelname = submodelnames{isub};
            model.(submodelname) = convertModelForJuliaBridge_(model.(submodelname));

        end

    end
    
end
    

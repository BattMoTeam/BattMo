function model = convertModelGrids(model)

    if isprop(model, 'G') && ~isempty(model.G)
        G = model.G.mrstFormat;
        rock.perm = ones(G.cells.num, 1);
        rock.poro = ones(G.cells.num, 1);
        op = setupOperatorsTPFA(G, rock);
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


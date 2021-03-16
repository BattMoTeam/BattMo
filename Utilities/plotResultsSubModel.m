function plotResultsSubModel(model, name, states, varargin)
    %% should probably remove all non cell arrays

    plotToolbar(model.(name).G, cellfun(@(x) x.(name), states, 'UniformOutput', false), varargin{:})

end
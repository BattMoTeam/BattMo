 function plotSubModel(model, name,states, varargin)
    %% should proably remove all non cell arrays
    plotToolbar(model.(name).G,cellfun(@(x) x.(name),states,'UniformOutput',false),varargin{:})
 end
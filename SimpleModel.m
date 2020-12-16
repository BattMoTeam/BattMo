classdef SimpleModel < PhysicalModel
    

    methods
        function model = SimpleModel(varargin)
            model = model@PhysicalModel([]);
            model.debugmode = false;
            model = merge_options(model, varargin{:});
        end

        function varnames = getPrimaryVarNames(model)
        % list all the primarty variables that are recognized and can be
        % handled by the model, used by the updateState member function.
            varnames = model.primaryVarNames;
        end

    end

end

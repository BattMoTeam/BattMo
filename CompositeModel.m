classdef CompositeModel < SimpleModel

    properties
        SubModels;
        SubModelNames;
        nSubModels;
        
        isCompositeModel;
    end

    methods
        
        function model = CompositeModel(varargin)
        % The constructor function should be complemented so that the properties
        % SubModels, SubModelNames are defined and the function
        % setupCompositeModel is called.
            model = model@SimpleModel(varargin{:});
        end
        
        function ind = getSubModelInd(model, name)
            ind = strcmp(name, model.SubModelNames);
        end
        
        function submodel = getSubModel(model, name)
            ind = model.getSubModelInd(name);
            submodel = model.SubModels{ind};
        end

        function model = setSubModel(model, name, submodel)
            ind = model.getSubModelInd(name);
            model.SubModels{ind} = submodel;
        end

        function model = setupCompositeModel(model)
            nsubmodels = numel(model.SubModels);
            model.nSubModels     = nsubmodels;
            model.isCompositeModel = true;
        end

        
        function primaryVarNames = getPrimaryVarNames(model)
            nsubmodels = model.nSubModels;
            primaryVarNames = {};
            for i = 1 : nsubmodels
                submodel          = model.SubModels{i};
                submodel          = submodel.setPrimaryVarNames();
                model.SubModels{i} = submodel;

                submodelprimaryvarnames = submodel.getPrimaryVarNames();
                primaryVarNames         = horzcat(primaryVarNames, ...
                                                  submodelprimaryvarnames);
            end
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            nsubmodels = model.nSubModels;

            found = false;
            i = 1;
            while ~found & (i <= nsubmodels)
                submodel = model.SubModels{i};
                [fn, index] = submodel.getVariableField(name);
                if ~isempty(index)
                    found = true;
                end
                i = i + 1;
            end
            assert(found, 'unknown variable');
        end

        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)

            % Due to current use of function updateState, we need to
            % reinitiate the primary variable (this is unfortunate).
            model = model.setPrimaryVarNames();
            nsubmodels = model.nSubModels;
            for i = 1 : nsubmodels
                submodel   = model.SubModels{i};
                [state, ~] = submodel.updateState(state, problem, dx, []);
            end
            report = [];
        end

        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            nsubmodels = model.nSubModels;
            for i = 1 : nsubmodels
                submodel = model.SubModels{i};
                [state, ~] = submodel.updateAfterConvergence(state0, state, ...
                                                             dt, ...
                                                             drivingForces);
            end
            report = [];
        end
        
    end
    

end

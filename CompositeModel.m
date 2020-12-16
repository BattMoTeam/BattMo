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
        % initiateCompositeModel is called.
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
            submodel.useglobalnames = true;
            model.SubModels{ind} = submodel;
        end

        function model = initiateCompositeModel(model)
            nsubmodels = numel(model.SubModels);
            model.nSubModels = nsubmodels;
            model.isCompositeModel = true;
            for ind = 1 : nsubmodels
                model.SubModels{ind}.useglobalnames = true; 
            end
        end
        
        function [globalnames, localnames] = getModelPrimaryVarNames(model)
            
            nsubmodels = model.nSubModels;
            localnames = {};
            for i = 1 : nsubmodels
                submodel = model.SubModels{i};
                submodelvarnames = submodel.getModelPrimaryVarNames();
                localnames = horzcat(localnames, submodelvarnames);
            end
            globalnames = cellfun(@(name) (model.setupGlobalName(name)), localnames, 'uniformoutput', false);
            
        end
        
        function [globalnames, localnames] = getModelVarNames(model)
        % default for compositemodel : fetch all the defined names in the submodels
            nsubmodels = model.nSubModels;
            localnames = {};
            for i = 1 : nsubmodels
                submodel = model.SubModels{i};
                submodelvarnames = submodel.getModelVarNames();
                localnames = horzcat(localnames, submodelvarnames);
            end
            globalnames = cellfun(@(name) (model.setupGlobalName(name)), localnames, 'uniformoutput', false);
        end
        
        
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            waring('to be updated');
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
            waring('to be updated');
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

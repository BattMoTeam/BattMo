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
            model.nSubModels = nsubmodels;
            model.isCompositeModel = true;
        end
        
        function primaryVarNames = getPrimaryVarNames(model)
            nsubmodels = model.nSubModels;
            
            primaryVarNames = model.getModelPrimaryVarNames;
            for i = 1 : nsubmodels
                submodel = model.SubModels{i};
                submodelname = submodel.getModelName();
                submodelprimaryvarnames = submodel.getPrimaryVarNames();
                submodelprimaryvarnames = cellfun(@(str) (sprintf('%s_%s', submodelname, varname)), submodelprimaryvarnames)
                primaryVarNames = horzcat(primaryVarNames, submodelprimaryvarnames);
            end
            
        end

        function varNames = getVarNames(model)
            nsubmodels = model.nSubModels;
            
            varNames = model.getModelVarNames();
            for i = 1 : nsubmodels
                submodel = model.SubModels{i};
                submodelname = submodel.getModelName();
                submodelvarnames = submodel.getPrimaryVarNames();
                submodelvarnames = cellfun(@(str) (sprintf('%s_%s', modelname, submodelname, varname)), submodelvarnames)
                varNames = horzcat(VarNames, submodelvarnames);
            end
            
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

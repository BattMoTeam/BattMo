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
            submodel.usenamespace = true;
            model.SubModels{ind} = submodel;
        end

        function model = initiateCompositeModel(model)
            nsubmodels = numel(model.SubModels);
            model.nSubModels = nsubmodels;
            model.isCompositeModel = true;
            for ind = 1 : nsubmodels
                submodel = model.SubModels{ind};
                submodel.usenamespace = true; 
                
                if model.usenamespace
                    subnamespace = sprintf('%s_%s', model.namespace, submodel.namespace);
                    submodel.namespace = subnamespace;
                end
                
                if isa(submodel, 'CompositeModel')
                    submodel = submodel.initiateCompositeModel();
                end
                model.SubModels{ind} = submodel;
            end
        end
        
        function [namespaces, names] = getModelPrimaryVarNames(model)
        
            nsubmodels = model.nSubModels;
            namespaces = {};
            names = {};
            for i = 1 : nsubmodels
                submodel = model.SubModels{i};
                [namespaces1, names1] = submodel.getModelPrimaryVarNames();
                namespaces = horzcat(namespaces, namespaces1);
                names = horzcat(names, names1);
            end
            
        end
        
        function [namespaces, names] = getModelVarNames(model)
        % default for compositemodel : fetch all the defined names in the submodels
        
            nsubmodels = model.nSubModels;
            namespaces = {};
            names = {};
            for i = 1 : nsubmodels
                submodel = model.SubModels{i};
                [namespaces1, names1] = submodel.getModelVarNames();
                namespaces = horzcat(namespaces, namespaces1);
                names = horzcat(names, names1);
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

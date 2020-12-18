classdef CompositeModel < SimpleModel

    properties
        SubModels;
        SubModelNames;
        nSubModels;
        
        isCompositeModel;
    end

    methods
        
        function model = CompositeModel(name, varargin)
        % The constructor function should be complemented so that the properties
        % SubModels, SubModelNames are defined and the function
        % initiateCompositeModel is called.
            model = model@SimpleModel(name, varargin{:});
            model.isnamespaceroot = false;
        end
        
        function ind = getSubModelInd(model, name)
            ind = strcmp(name, model.SubModelNames);
            if all(ind == 0)
                error('submodel not found');
            end
            ind = find(ind);
            if numel(ind) > 1
                error('submodels have same name');
            end
        end
        
        function submodel = getSubModel(model, name)
            ind = model.getSubModelInd(name);
            submodel = model.SubModels{ind};
        end

        function model = setSubModel(model, submodel, name)
            ind = getSubModelInd(model, name);
            model.SubModels{ind} = submodel;
        end
        
        
        function model = initiateCompositeModel(model)
            nsubmodels = numel(model.SubModels);
            model.nSubModels = nsubmodels;
            model.isCompositeModel = true;
            
            if model.isnamespaceroot
                model.namespace = {};
            end
            
            % Setup the namespaces for all the submodels
            for ind = 1 : nsubmodels
                submodel = model.SubModels{ind};
                submodel.isnamespaceroot = false; 
                submodelname = submodel.getModelName();
                
                if ~model.isnamespaceroot
                    subnamespace = horzcat(model.namespace, {submodelname});
                else
                    subnamespace = {submodelname};
                end
                submodel.namespace = subnamespace;
                
                if isa(submodel, 'CompositeModel')
                    submodel = submodel.initiateCompositeModel();
                end
                
                model.SubModels{ind} = submodel;
            end
            
        end
        
        function varnames = getModelPrimaryVarNames(model)
        
            varnames = model.pvarnames;
            nsubmodels = model.nSubModels;
            for i = 1 : nsubmodels
                submodel = model.SubModels{i};
                varnames1 = submodel.getModelPrimaryVarNames();
                varnames = horzcat(varnames, varnames1);
            end
            
        end
        
        function varnames = getModelVarNames(model)
        % default for compositemodel : fetch all the defined names in the submodels
        
            varnames = model.varnames;
            nsubmodels = model.nSubModels;
            for i = 1 : nsubmodels
                submodel = model.SubModels{i};
                varnames1 = submodel.getModelVarNames();
                varnames = horzcat(varnames, varnames1);
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

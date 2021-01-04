classdef RootModel < CompositeModel
    
    properties
        
        modeldict
        varfieldnames
        varindices
        propfunctdict
        propmodeldict
        
    end
    
    methods
        
        function model = RootModel(name, varargin)
            model = model@CompositeModel(name, varargin{:})
            model.modeldict     = containers.Map;
            model.varfieldnames = containers.Map;
            model.varindices    = containers.Map;
            model.propfunctdict = containers.Map;
            model.propmodeldict = containers.Map;
        end
        
        function model = initiateRootModel(model)
        % set root model
            
            
            % add root in all the submodels
            model = model.addRootModel(model);
          
            modeldict = model.modeldict;
            varfieldnames = model.varfieldnames;
            varindices = model.varindices;
            propfunctdict = model.propfunctdict;
            propmodeldict = model.propmodeldict;
            
            % update model dictionary
            modeldict = updateModelDict(modeldict, model);
            % update variable fieldname dictionary
            [varfieldnames, varindices] = updateVarDict(varfieldnames, varindices, model);
            % update property function dictionary
            [propfunctdict, propmodeldict] = updateFunctionDict(propfunctdict, propmodeldict, model);
            
            model.modeldict = modeldict;
            model.varfieldnames = varfieldnames;
            model.varindices = varindices;
            model.propfunctdict = propfunctdict;
            model.propmodeldict = propmodeldict;
            
        end
        
        
    end
end

classdef SimpleModel < PhysicalModel
    
    properties
        useglobalnames;
    end
    
    
    methods
        function model = SimpleModel(varargin)
            model = model@PhysicalModel([]);
            model.debugmode = false;
            model = merge_options(model, varargin{:});
            model.parentname = '';
        end

        function modelname = getModelName();
        end
        
        function varnames = getPrimaryVarNames(model)
        % list all the primarty variables that are recognized and can be
        % handled by the model, used by the updateState member function.
        end
        
        function varnames = getPrimaryVarNames()
        % for a simple model the primary variable names only consist of those of the model (no child)
            varnames = model.getModelPrimaryVarNames();
        end
        
        function varnames = getVarNames()
        % for a simple model the variable names only consist of those of the model (no child)
            varnames = model.getModelVarNames();
        end
        
        function varnames = getModelPrimaryVarNames(model)
        % Own primary variables of the composite model
        end
        
        function varnames = getModelVarNames(model)
        % Own variables of the composite model (not including the primary ones)
        end

        
        function [fn, index] = getVariableField(model, name, varargin)
            
            primaryVarNames = model.getPrimaryVarNames();
            varNames = model.getVarNames();
            
            varNames = horzcat(primaryVarNames, varNames);
            
            isok = ismember(name, varNames);
            assert(isok, 'unknown variables');
            
            fn = name;
            index = 1;
            
        end
        
        function [globalnames, localnames] = variableNameMapping(model)
        % to know the name of the variable seen from the environment
        % default is to use the model name
            varnames1 = model.getModelVariableNames();
            varnames2 = model.getModelPrimaryVariableNames();
            localnames = horzcat(varnames1, varnames2);
            modelname = model.getModelName();
            
            globalnames = cellfun(@(name) (sprintf('%s_%s', modelname, name)), localnames); 
        end
        

        
    end

end

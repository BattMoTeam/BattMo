classdef SimpleModel < PhysicalModel
    
    properties
        useglobalnames;
    end
    
    
    methods
        function model = SimpleModel(varargin)
            model = model@PhysicalModel([]);
            model = merge_options(model, varargin{:});
        end

        function modelname = getModelName(model);
        % This is a virtual method that should be implemented by the model.
            error('virtual method');
        end
        
        function [globalnames, localnames] = getModelPrimaryVarNames(model)
        % List all the primarty variables that are recognized and can be handled by the model, used by the updateState member
        % function.
            localnames = {};
            globalnames = model.setupGlobalNames(localnames);
        end
        
        
        function [globalnames, localnames] = getModelVarNames(model)
        % for a simple model the variable names only consist of those of the model (no child)
            localnames = {};
            globalnames = model.setupGlobalNames(localnames);
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
           
            [globalnames, localnames] = model.getVarNames();
            
            [isok, ind] = ismember(name, localnames);
            assert(isok, 'unknown variables');
            
            fn = globalnames{ind};
            index = 1;
            
        end
        
        function globalnames = setupGlobalNames(model, localnames)
             globalnames = cellfun(@(name) (model.setupGlobalName(name)), localnames, 'uniformoutput', false);            
        end
        
        function globalname = setupGlobalName(model, localname);
        % Default function to setup global names : we prepend the model name.
            if model.useglobalnames
                modelname = model.getModelName();
                globalname = sprintf('%s_%s', modelname, localname);
            else
                globalname = localname;
            end
        end
                
        function [globalnames, localnames] = getVarNames(model)
        % The variables consists of the model variables and the primary variables (they are handled separately)
        % We can add more variables there and make them recognizable by the model.
            
            [globalnames1, localnames1] = model.getModelVarNames();
            [globalnames2, localnames2] = model.getModelPrimaryVarNames();
            
            globalnames = horzcat(globalnames1, globalnames2);
            localnames = horzcat(localnames1, localnames2);
           
        end
        
    end

end

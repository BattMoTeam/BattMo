classdef SimpleModel < PhysicalModel
    
    properties
        usenamespace;
        namespace;
    end
    
    
    methods
        function model = SimpleModel(varargin)
            model = model@PhysicalModel([]);
            model = merge_options(model, varargin{:});
            model.usenamespace = false; % default choice
            model.namespace = model.getModelName;
        end

        function state = validateState(model, state)
            
            [namespaces, names] = model.getVarNames();
            varnames = model.addNameSpace(namespaces, names);
            
            for i = 1 : numel(varnames)
                varname = varnames{i};
                if ~isfield(state, varname)
                    state.(varname) = [];
                end
            end
            
        end
        
        function modelname = getModelName(model);
        % This is a virtual method that should be implemented by the model.
            error('virtual method');
        end
        
        function [namespaces, names] = getModelPrimaryVarNames(model)
        % List the primary variables. For SimpleModel the variable names only consist of those declared in the model (no child)
            names = {};
            namespaces = {};
        end
        
        
        function [namespaces, names] = getModelVarNames(model)
        % List the variable names (primary variables are handled separately). For SimpleModel the variable names only consist of
        % those declared in the model (no child)
            names = {};
            namespaces = {};
        end
        
        function [namespaces, names] = getVarNames(model)
        % Collect all the variable names (primary and the others). External variable names can be added there
            [namespaces1, names1] = model.getModelPrimaryVarNames;
            [namespaces2, names2] = model.getModelVarNames;
            namespaces = [namespaces1, namespaces2];
            names = [names1, names2];
        end

        function fullnames = addNameSpace(model, namespaces, names)
        % Return a name (string) which identify uniquely the pair (namespace, name). This is handled here by adding "_" at the
        % end of the namespace and join it to the name.
            if iscell(names)
                for ind = 1 : numel(namespaces)
                    fullnames{ind} = model.addNameSpace(namespaces{ind}, names{ind});
                end
            else
                if ~isempty(namespaces)
                    fullnames = sprintf('%s_%s', namespaces, names);
                else
                    fullnames = names;
                end                    
            end
        end
        
        
        function [fn, index] = getVariableField(model, name, varargin)
        % In this function the variable name is associated to a field name (fn) which corresponds to the field where the
        % variable is stored in the state.  See PhysicalModel
            [namespaces, names] = model.getVarNames();
            
            [isok, ind] = ismember(name, names);
            assert(isok, 'unknown variables');
            
            if model.usenamespace
                fn = model.addNameSpace(namespaces{ind}, name);
            else
                fn = name;
            end
            
            index = 1;
            
        end

        function namespaces = assignCurrentNameSpace(model, names)
        % utility function which returns a cell consisting of the current model namespace with the same size as names.
            namespace = model.namespace;
            n = numel(names);
            namespaces = repmat({namespace}, 1, n);
        end
        
    end

end

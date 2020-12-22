classdef SimpleModel < PhysicalModel
    
    properties
        parentmodel
        hasparent

        modelname
        
        % Own variables
        namespace
        
        names
        pnames
        
        aliases
        
        % Variables that are added externally
        varnames
        pvarnames
        
    end
    
    
    methods
        function model = SimpleModel(modelname, varargin)
            model = model@PhysicalModel([]);
            model = merge_options(model, varargin{:});
            model.modelname = modelname;
            model.namespace = {model.getModelName()};
            model.hasparent = false;
            % no variables listed by default
            model.names     = {};
            model.pnames    = {};
            model.aliases   = {};
            model.varnames  = {};
            model.pvarnames = {};
        end

        function state = validateState(model, state)
            
            varnames = model.getVarNames();
            
            for i = 1 : numel(varnames)
                varname = varnames{i};
                name = varname.fullname;
                if ~isfield(state, name)
                    state.(name) = [];
                end
            end
            
        end
        
        function modelname = getModelName(model);
            modelname = model.modelname;
        end
        
        function parentmodel = getParentModel(model)
             parentmodel = model.parentmodel;
        end
        
        function varnames = getModelPrimaryVarNames(model)
        % List the primary variables. For SimpleModel the variable names only consist of those declared in the model (no child)
            varnames1 = model.assignCurrentNameSpace(model.pnames);
            varnames2 = model.pvarnames;
            varnames = horzcat(varnames1, varnames2);
        end
        
        
        function varnames = getModelVarNames(model)
        % List the variable names (primary variables are handled separately). For SimpleModel the variable names only consist of
        % those declared in the model (no child)
            varnames1 = model.assignCurrentNameSpace(model.names);
            varnames2 = model.varnames;
            varnames = horzcat(varnames1, varnames2);
        end
        
        function varnames = getVarNames(model)
        % Collect all the variable names (primary and the others). External variable names can be added there
            varnames1 = model.getModelPrimaryVarNames;
            varnames2 = model.getModelVarNames;
            varnames = horzcat(varnames1, varnames2);
        end

        function names = getFullVarNames(model)
        % Returns full names (i.e. including namespaces)
            varnames = model.getVarNames();
            names = {};
            for i = 1 : numel(varnames)
                names{end + 1} = varnames{i}.fullname;
            end
        end
        
        function fullnames = addNameSpace(model, varnames)
        % Return a name (string) which identify uniquely the pair (namespace, name). This is handled here by adding "_" at the
        % end of the namespace and join it to the name.
            if iscell(names)
                for ind = 1 : numel(namespaces)
                    fullnames{ind} = model.addNameSpace(namespaces{ind}, names{ind});
                end
            else
                fullnames = varnames.join();
            end
        end
        
        function [isalias, varname] = setupVarName(model, name)
        % check if name is an alias
            aliases = model.aliases;
            isalias = false;
            varname = [];
            
            if isempty(aliases)
                return
            end
            
            for ind = 1 : numel(aliases)
                alias = aliases{ind};
                if strcmp(name, alias{1})
                    varname = alias{2};
                    isalias = true;
                    return
                end
            end
        end
            
        function [fn, index] = getVariableField(model, name, varargin)
        % In this function the variable name is associated to a field name (fn) which corresponds to the field where the
        % variable is stored in the state.  See PhysicalModel
            
            % Check if there exist an alias
            [isalias, varname] = model.setupVarName(name);
            
            if isalias & model.hasparent
                % Call alias
                parentmodel = model.getParentModel();
                [fn, index] = parentmodel.getVariableField(varname);
            else
                % Check that name is declared in model.names
                isok = ismember(name, model.names);
                % Otherwise it can be an alias with empty namespace (we do not check for empty name space)
                isok = isok | isalias;
                assert(isok, 'name is not declared/recognized by the model');
                
                % Construct name from namespace
                namespace = model.namespace;
                varname = VarName(namespace, name);
                fn = varname.getfieldname();
                index = 1;
            end
        end
        
        
        function varnames = assignCurrentNameSpace(model, names)
        % utility function which returns a cell consisting of the current model namespace with the same size as names.
            namespace = model.namespace;
            n = numel(names);
            varnames = {};
            for ind = 1 : n
                varnames{end + 1} = VarName(namespace, names{ind});
            end
        end
        
    end

end

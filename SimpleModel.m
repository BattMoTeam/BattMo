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

        varfunctions
    end
    
    
    methods
        function model = SimpleModel(modelname, varargin)
            model = model@PhysicalModel([]);
            model = merge_options(model, varargin{:});
            model.modelname = modelname;
            % By default, the model has no parent so that the namespace is empty
            model.namespace = {};
            model.hasparent = false;
            % no variables listed by default
            model.names     = {};
            model.pnames    = {};
            model.aliases   = {};
        end

        function state = validateState(model, state)
            
            varnames = model.getModelVarNames();
            
            for i = 1 : numel(varnames)
                varname = varnames{i};
                name = varname.getfieldname;
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
         
        function submodel = getAssocModel(model, name)
            if isa(name, 'char')
                if strcmp(name, '..')
                    submodel = model.parentmodel;
                elseif strcmp(name, '.')
                    submodel = model;
                end
            elseif isa(name, 'cell') && numel(name) > 1
                firstparentname = name{1};
                name = name{2 : end};
                submodel = model.getAssocModel(firstparentname);
                submodel = submodel.getAssocModel(name);
            elseif isa(name, 'cell') && numel(name) == 1
                name = name{1};
                submodel = model.getAssocModel(name);
            else
                error('name type not recognized');
            end
        end
        
        function varnames = getModelPrimaryVarNames(model)
        % List the primary variables. For SimpleModel the variable names only consist of those declared in the model (no child)
            varnames = model.assignCurrentNameSpace(model.pnames);
        end
        
        
        function varnames = getModelVarNames(model)
        % List the variable names (primary variables are handled separately). For SimpleModel the variable names only consist of
        % those declared in the model (no child)
            varnames = model.assignCurrentNameSpace(model.names);
        end

        function names = getFullVarNames(model)
        % Returns full names (i.e. including namespaces)
            varnames = model.getModelVarNames();
            names = {};
            for i = 1 : numel(varnames)
                names{end + 1} = varnames{i}.getfieldname;
            end
        end
        
        function names = getFullPrimaryVarNames(model)
        % Returns full names (i.e. including namespaces)
            varnames = model.getModelPrimaryVarNames();
            names = {};
            for i = 1 : numel(varnames)
                names{end + 1} = varnames{i}.getfieldname;
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

        function model = addAlias(model, alias)
            aliases = model.aliases;
            aliases{end + 1} = alias;
            model.aliases = aliases;
        end
        
        
        function model = setVarFunction(model, varfunction)
        % If there exist a name in varfunctions that correspond to varfunction, overwrite this entry with
        % varfunction. Otherwise, add varfunction to varfunctions.
            isdefined = false;
            name = varfunction{1};
            varfunctions = model.varfunctions;
            for ind = 1 : numel(varfunctions)
                if strcmp(name, varfunctions{ind}{1})
                    isdefined = true;
                    break
                end
            end
            
            if isdefined
                varfunctions{ind} = varfunction;
            else
                varfunctions{end + 1} = varfunction;
            end
            
            model.varfunctions = varfunctions;
            
        end
        
        
        function [isalias, varname] = aliasLookup(model, name)
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
            
        
        function [val, state] = getUpdatedProp(model, state, name)
            state = model.updateProp(state, name);
            val = model.getProp(state, name);
        end
        
        
        function state = updateProp(model, state, name)
            val = model.getProp(state, name);
            
            % check for vector variables
            isupdated = false;
            
            if ~isempty(val) 
                isupdated = true;
            end
            
            % Check for cell variables
            if iscell(val)
                isupdated = true;
                for ind = 1 : numel(val)
                    if isempty(val{ind})
                        isupdated = false;
                        break
                    end
                end
            end
            
            if isupdated
                return
            end
            
            % We look for an updating function for the variable
            [isalias, varname] = model.aliasLookup(name);
            if isalias
                fnmodel = model.getAssocModel(varname.namespace);
                state  = fnmodel.updateProp(state, varname.name);
            end
            
            % find the updating property
            varfunctions = model.varfunctions;
            updatefn = [];
            for ind = 1 : numel(varfunctions)
                varfunction = varfunctions{ind};
                if strcmp(varfunction{1}, name)
                    updatefn = varfunction{2}{1};
                    fnmodelname = varfunction{2}{2};
                    fnmodel = model.getAssocModel(fnmodelname);
                    break
                end
            end
            state = updatefn(fnmodel, state);
            
        end

        function [fn, index] = getVariableField(model, name, throwError, varargin)
        % In this function the variable name is associated to a field name (fn) which corresponds to the field where the
        % variable is stored in the state.  See PhysicalModel
            
            opt = struct('index', []);
            opt = merge_options(opt, varargin{:});
            if isempty(opt.index)
                index = ':';
            else
                index = opt.index;
            end
            
            if isa(name, 'VarName')
                varname = name;
                varmodel = model.getAssocModel(varname.namespace);
                [fn, index] = varmodel.getVariableField(varname.name, true, 'index', varname.index);
            elseif iscell(name)
                % syntaxic sugar (do not need to setup VarName)
                varname = VarName({name{1 : end - 1}}, name{end});
                [fn, index] = model.getVariableField(varname);
            else
                % Check if there exist an alias
                [isalias, varname] = model.aliasLookup(name);

                if isalias
                    [fn, index] = model.getVariableField(varname);
                else
                    % Check that name is declared
                    isok = ismember(name, model.names);
                    % Otherwise it can be an alias refering to an external field 
                    assert(isok, 'name is not declared/recognized by the model');
                    
                    % Construct name from namespace
                    namespace = model.namespace;
                    varname = VarName(namespace, name);
                    fn = varname.getfieldname();
                    index = index;
                end
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

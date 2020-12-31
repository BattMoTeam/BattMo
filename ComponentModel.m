classdef ComponentModel < PhysicalModel
    
    properties
        parentmodel
        hasparent

        modelname
        
        % Own variables
        namespace
        
        names
        pnames
        
        aliases

        propfunctions
    end
    
    
    methods
        function model = ComponentModel(modelname, varargin)
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

        function state = initiateState(model, state)
            varnames = model.assignCurrentNameSpace(model.names);
            for i = 1 : numel(varnames)
                varname = varnames{i};
                name = varname.getfieldname();
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

        function name = getModelFullName(model)
            if model.hasparent
                name = join(model.namespace, '_');
                name = name{1};
            else
                name = 'root';
            end
        end
        
        function submodel = getAssocModel(model, name)
            if isempty(name)
                submodel = model;
            elseif isa(name, 'char')
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
        % List the primary variables. For ComponentModel the variable names only consist of those declared in the model (no child)
            varnames = model.assignCurrentNameSpace(model.pnames);
        end
        
        
        function varnames = getModelVarNames(model)
        % List the variable names (primary variables are handled separately). For ComponentModel the variable names only consist of
        % those declared in the model (no child)
            varnames = model.assignCurrentNameSpace(model.names);
        end

        function names = getVarFieldNames(model)
        % Returns full names (i.e. including namespaces)
            varnames = model.getModelVarNames();
            names = {};
            for i = 1 : numel(varnames)
                names{end + 1} = varnames{i}.getfieldname;
            end
        end
        
        function names = getPrimaryVarFieldNames(model)
        % Returns full names (i.e. including namespaces)
            varnames = model.getModelPrimaryVarNames();
            names = {};
            for i = 1 : numel(varnames)
                names{end + 1} = varnames{i}.getfieldname;
            end
        end
        
        function model = addAlias(model, alias)
            aliases = model.aliases;
            aliases{end + 1} = alias;
            model.aliases = aliases;
        end
        
        
        function model = setPropFunction(model, propfunction)
        % If there exist a name in propfunctions that correspond to propfunction, overwrite this entry with
        % propfunction. Otherwise, add propfunction to propfunctions.
            isdefined = false;
            name = propfunction.name;
            propfunctions = model.propfunctions;
            for ind = 1 : numel(propfunctions)
                if strcmp(name, propfunctions{ind}.name)
                    isdefined = true;
                    break
                end
            end
            
            if isdefined
                propfunctions{ind} = propfunction;
            else
                propfunctions{end + 1} = propfunction;
            end
            
            model.propfunctions = propfunctions;
            
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
                return
            end
            
            % find the updating property
            propfunctions = model.propfunctions;
            updatefn = [];
            fnmodel = [];
            for ind = 1 : numel(propfunctions)
                propfunction = propfunctions{ind};
                if strcmp(propfunction.name, name)
                    updatefn = propfunction.fn;
                    fnmodelname = propfunction.namespace;
                    fnmodel = model.getAssocModel(fnmodelname);
                    break
                end
            end
            
            assert(~isempty(updatefn), 'property is empty and no property function is provided.')
            state = updatefn(fnmodel, state);
            
        end

        function [fn, index] = getVariableField(model, name, index)
        % In this function the variable name is associated to a field name (fn) which corresponds to the field where the
        % variable is stored in the state.  See PhysicalModel
            
            if nargin < 3
                index = ':';
            end
            
            if isa(name, 'VarName')
                varname = name;
                varmodel = model.getAssocModel(varname.namespace);
                [fn, index] = varmodel.getVariableField(varname.name, varname.index);
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

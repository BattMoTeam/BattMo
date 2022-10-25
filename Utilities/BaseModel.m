classdef BaseModel < PhysicalModel

    properties

        propertyFunctionList
        varNameList
        subModelNameList

        % used only for a root model to cleanup output of some of the functions in ComputationalGraphTool
        staticVarNameList
        
    end
        
    methods

        function model = BaseModel()

            model = model@PhysicalModel([]);
            
            model.propertyFunctionList = {};
            model.varNameList          = {};
            model.subModelNameList     = {};
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = model.registerSubModels();
            
        end

        
        function model = registerVarName(model, varname)
            if isa(varname, 'VarName')
                model.varNameList = mergeList(model.varNameList, {varname});
            elseif isa(varname, 'char')
                varname = VarName({}, varname);
                model = model.registerVarName(varname);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                model = model.registerVarName(varname);
            else
                error('varname not recognized');
            end
        end
        
        function model = registerVarNames(model, varnames)
            for ivarnames = 1 : numel(varnames)
                model = model.registerVarName(varnames{ivarnames});
            end
        end


        function model = registerStaticVarName(model, varname)
            if isa(varname, 'VarName')
                model.staticVarNameList = mergeList(model.staticVarNameList, {varname});
            elseif isa(varname, 'char')
                varname = VarName({}, varname);
                model = model.registerStaticVarName(varname);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                model = model.registerStaticVarName(varname);
            else
                error('varname not recognized');
            end
        end
        
        function model = registerStaticVarNames(model, varnames)
            for ivarnames = 1 : numel(varnames)
                model = model.registerStaticVarName(varnames{ivarnames});
            end
        end

        
        function model = removeVarNames(model, varnames)
            for ivar = 1 : numel(varnames)
                model = model.removeVarName(varnames{ivar});
            end
        end
        
        function model = removeVarName(model, varname)

            if isa(varname, 'char')
                varname = VarName({}, varname);
                model = model.removeVarName(varname);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                model = model.removeVarName(varname);
            elseif isa(varname, 'VarName')
                % remove from varNameList
                varnames = model.varNameList;
                nvars = numel(varnames);
                keep = true(nvars, 1);
                for ivar = 1 : nvars
                    if varnames{ivar} == varname
                        keep(ivar) = false;
                        break
                    end
                end
                model.varNameList = varnames(keep);
                
                % remove from propertyFunctionList if varname occurs in output and input variables
                propfuncs = model.propertyFunctionList;
                nprops = numel(propfuncs);
                keep = true(nprops, 1);
                for iprop = 1 : nprops
                    propfunc = propfuncs{iprop};
                    if propfunc.varname == varname
                        keep(iprop) = false;
                        break
                    end
                    inputvarnames = propfunc.inputvarnames;
                    for iinput = 1 : numel(inputvarnames)
                        if inputvarnames{iinput} == varname
                            keep(iprop) = false;
                            break
                        end
                    end
                end
                model.propertyFunctionList = propfuncs(keep);                

            else
                error('varname not recognized');
            end
            
            
        end
        
        
        function model = registerPropFunction(model, propfunc)
            if isa(propfunc, 'PropFunction')
                model.propertyFunctionList = mergeList(model.propertyFunctionList, {propfunc});
            elseif isa(propfunc, 'cell')
                assert(numel(propfunc) == 3, 'format of propfunc not recognized');
                varname = propfunc{1};
                if isa(varname, 'VarName')
                    % ok
                elseif isa(varname, 'char')
                    varname = VarName({}, varname);
                elseif isa(varname, 'cell')
                    varname = VarName(varname(1 : end - 1), varname{end});
                else
                    error('format of propfunc not recognized');
                end
                fn = propfunc{2};
                assert(isa(propfunc{3}, 'cell'), 'format of propfunc not recognized');
                inputvarnames = cell(1, numel(propfunc{3}));
                for iinputvarnames = 1 : numel(propfunc{3})
                    inputvarname = propfunc{3}{iinputvarnames};
                    if isa(inputvarname, 'VarName')
                        % ok
                    elseif isa(inputvarname, 'char')
                        inputvarname = VarName({}, inputvarname);
                    elseif isa(inputvarname, 'cell')
                        inputvarname = VarName(inputvarname(1 : end - 1), inputvarname{end});
                    else
                        error('format not recognized')
                    end
                    inputvarnames{iinputvarnames} = inputvarname;
                end
                modelnamespace = {};
                propfunc = PropFunction(varname, fn, inputvarnames, modelnamespace);
                model = model.registerPropFunction(propfunc);
            else
                error('propfunc not recognized');
            end                
        end

                
        function stateAD = initStateAD(model, state)
        % initialize a new cleaned-up state with AD variables
            
            pnames  = model.getPrimaryVariables();
            vars = cell(numel(pnames),1);
            for i=1:numel(pnames)
                vars{i} = model.getProp(state,pnames{i});
            end
            % Get the AD state for this model           
            [vars{:}] = model.AutoDiffBackend.initVariablesAD(vars{:});
            stateAD =struct();
            for i=1:numel(pnames)
               stateAD = model.setNewProp(stateAD, pnames{i}, vars{i});
            end

            stateAD = model.addStaticVariables(stateAD, state);
            
        end 
        
        
        function model = registerSubModels(model)
            
            propfuncs = model.propertyFunctionList;
            varnames = model.varNameList;

            if ~isempty(model.subModelNameList)
                
                submodelnames = model.subModelNameList;
                
            else
                
                props = propertynames(model);
                submodelnames = {};
                
                for iprops = 1 : numel(props)
                    prop = props{iprops};
                    if isa(model.(prop), 'BaseModel')
                        submodelnames{end + 1} = prop;
                    end
                end
                model.subModelNameList = submodelnames;
                
            end
            
            for isub = 1 : numel(submodelnames)

                submodelname = submodelnames{isub};
                submodel = model.(submodelname);
                
                submodel = registerVarAndPropfuncNames(submodel);

                % Register the variable names from the submodel after adding the model name in the name space

                subvarnames = submodel.varNameList;
                
                for isubvar = 1 : numel(subvarnames)
                    subvarname = subvarnames{isubvar};
                    subvarname.namespace = {submodelname, subvarname.namespace{:}};
                    subvarnames{isubvar} = subvarname;
                end
                
                varnames = mergeList(varnames, subvarnames);
                
                % Register the property functions from the submodel after adding the model name in the name space
                
                subpropfuncs = submodel.propertyFunctionList;

                for isubpropfunc = 1 : numel(subpropfuncs)

                    subpropfunc = subpropfuncs{isubpropfunc};
                    
                    subpropfunc.varname.namespace = {submodelname, subpropfunc.varname.namespace{:}};
                    subpropfunc.modelnamespace = {submodelname, subpropfunc.modelnamespace{:}};
                    
                    subinputvarnames = subpropfunc.inputvarnames;
                    for isubinput = 1 : numel(subinputvarnames)
                        subinputvarname = subinputvarnames{isubinput};
                        subinputvarname.namespace = {submodelname, subinputvarname.namespace{:}};
                        subinputvarnames{isubinput} = subinputvarname;
                    end
                    subpropfunc.inputvarnames = subinputvarnames;
                    
                    subpropfuncs{isubpropfunc} = subpropfunc;
                end               
                
                propfuncs = mergeList(propfuncs, subpropfuncs);
                
            end
            
             model.varNameList = varnames;
             model.propertyFunctionList = propfuncs;
            
        end
        
        
        function state = setProp(model, state, names, val)
            if iscell(names) & (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                state.(name) = model.setProp(state.(name), names, val);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    state{name} = val;
                else
                    state.(name) = val;
                end
            else
                error('format not recognized');
            end
        end
        
        
        function state = setNewProp(model, state, names, val)
            if iscell(names) & (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                if(not(isfield(state,name)))
                    state.(name)=struct();
                end
                state.(name) = model.setNewProp(state.(name), names, ...
                                                               val);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    if(not(iscell(state)))
                        state = {};
                    end
                    state{name} = val;
                else
                    state.(name) = val;
                end
            else
                error('format not recognized');
            end
        end

        
        function var = getProp(model, state, names)
            if iscell(names) && (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                var = model.getProp(state.(name), names);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    var = state{name};
                else
                    var = state.(name);
                end
            else
                error('format not recognized');
            end
        end
        
        function submod = getSubmodel(model, names)
            submod = model.(names{1});
            for i=2:numel(names)
                submod = submod.(names{i});
            end
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            p = model.getPrimaryVariables();

            for i = 1 : numel(dx)
                val = model.getProp(state, p{i});
                val = val + dx{i};
                state = model.setProp(state, p{i}, val);
            end
           
            report = [];
            
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

            [state, report] = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
            p = model.getPrimaryVariables();
            cleanState = [];
            for ind = 1 : numel(p)
                cleanState = model.copyProp(cleanState, state, p{ind});
            end

            cleanState = model.addStaticVariables(cleanState, state);
            
            state = cleanState;
            report = [];
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state)
        % function to add static variables (not AD) on the cleanState, called in updateAfterConvergence and
        % initStateAD. Time is added by default here
            cleanState.time = state.time;
            
        end
        
        function state = copyProp(model, state, refState, names)
            
            if iscell(names) & (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                if isfield(state, name)
                    state.(name) = model.copyProp(state.(name), refState.(name), names);
                else
                    state.(name) = model.copyProp([], refState.(name), names);
                end
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    state{name} = refState{name};
                else
                    state.(name) = refState.(name);
                end
            else
                error('format not recognized');
            end
        end

        
        function scale = getScales()
            scale = [];
            error();
        end
        
        function [state, report] = updateStateNew(model, state, problem, dx, drivingForces)
            
            scales = model.getScales();
            for i = 1:numel(problem.primaryVariables)
                p = problem.primaryVariables{i};
                % Update the state
                scale=model.getProp(scales,p);
                if(isempty(scale))
                     state = model.updateStateFromIncrement(state, dx{i}, ...
                                                            problem, p);
                else
                    state = model.updateStateFromIncrement(state, dx{i}, ...
                                                           problem, p, ...
                                                           scale.relchangemax, ...
                                                           scale.abschangemax);
                    val = model.getProp(state, p);
                    val = max(val,scale.min);
                    val = min(val,scale.max);
                    state = model.setProp(state, p, val);
                end               
            end
            report = []
        end
        
        function state = reduceState(model, state, removeContainers)
            state = value(state, false);
        end
        
        function outputvars = extractGlobalVariables(model, states)
            ns = numel(states);
            outputvars = cell(1, ns);
        end
        
    end
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}

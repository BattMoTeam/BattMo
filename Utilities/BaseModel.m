classdef BaseModel < PhysicalModel

    properties
        propertyFunctionList
        varNameList
    end
        
    methods

        function model = BaseModel()
            model = model@PhysicalModel([]);
            model.propertyFunctionList = {};
            model.varNameList = {};
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
                    elseif isa(varname, 'iscell')
                        varname = VarName(varname(1 : end - 1), varname{end});
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
        
        
        function model = registerSubModels(model, submodelnames)
        % submodels is a list of submodel given as submodel and name. They should have been already assigned before this
        % function is called
            
            propfuncs = model.propertyFunctionList;
            varnames = model.varNameList;
            
            for isub = 1 : numel(submodelnames)

                submodelname = submodelnames{isub};
                submodel = model.(submodelname);
                
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
            cleanState.time = state.time;
            
            state = cleanState;
            report = [];
            
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
        
        function [state, report] = updateStateNew(model, state, problem, ...
                                                  dx, drivingForces)
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
            state = value(state);
        end
        
    end
    
end



%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}

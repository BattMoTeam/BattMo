classdef BaseModel < PhysicalModel

    properties

        propertyFunctionList % List of property functions (PropFunction instances) for the model
        varNameList          % List of variable names  (VarName instances) for the model
        subModelNameList     % List of submodels for the model

        % The variables listed in extraVarNameList are elements in VarNameList that are not used in assembly of the
        % residuals, those has been typically introduced for post-processing.
        extraVarNameList
       
        computationalGraph


        %% The following properties are set when a model is equiped for simulation. It is then a root model (as opposed to the sub-models)
        % See method equipModelForComputation

        isRootSimulationModel % boolean if the model is meant to be run for simulation
        
        funcCallList
        primaryVarNames
        equationVarNames
        equationNames
        equationTypes

        scalings % cell array. Each element is a cell pair. The first element in the pair is the variable name using the
                 % syntax of the getProp methods. The second element is the scaling coefficient. See methods applyScaling
                
    end
        
    methods

        function model = BaseModel()

            model = model@PhysicalModel([]);
            
            model.propertyFunctionList = {};
            model.varNameList          = {};
            model.subModelNameList     = {};
            model.extraVarNameList     = {};
            
        end

        function model = registerVarAndPropfuncNames(model)

        % The function will be overload by the model. First, the submodels are setup. Then, the varNameList and
        % propertyFunctionList and others lists are populated using the methods
        %
        % - registerVarName 
        % - registerVarNames
        % - registerPropFunction
        % - setAsStaticVarName
        % - setAsStaticVarNames
        % - setAsExtraVarName
        % - setAsExtraVarNames
        %
        % This function can be used to modify the same list but for the submodels.
        %
        % Then, we can also use the methods
        % - removeVarName
        % - removeVarNames
        % to remove variables defined in the submodels. The property function using those will then be also removed.
        %
            model = model.registerSubModels();
            
        end

        function model = setupComputationalGraph(model)

            model.computationalGraph = ComputationalGraph(model);
            
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


        function model = setAsStaticVarName(model, varname)
        % Register a static variable name. A static variable is a variable that is set up at the beginning of the
        % assembly and does not depend on the primary variables.

            if isa(varname, 'VarName')
                fn = [];
                inputnames = {};
                model = model.registerPropFunction({varname, fn, inputnames});
            elseif isa(varname, 'char')
                varname = VarName({}, varname);
                model = model.setAsStaticVarName(varname);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                model = model.setAsStaticVarName(varname);
            else
                error('varname not recognized');
            end
        end
        
        function model = setAsStaticVarNames(model, varnames)
        % Register several static variable name using setAsStaticVarName

            for ivarnames = 1 : numel(varnames)
                model = model.setAsStaticVarName(varnames{ivarnames});
            end
            
        end
        
        function model = setAsExtraVarName(model, varname)
            
            if isa(varname, 'char')
                varname = VarName({}, varname);
                model = model.setAsExtraVarName(varname);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                model = model.setAsExtraVarName(varname);
            elseif isa(varname, 'VarName')
                model.extraVarNameList = mergeList(model.extraVarNameList, {varname});
            end
        end
        
        function model = setAsExtraVarNames(model, varnames)
            for ivar = 1 : numel(varnames)
                model = model.setAsExtraVarName(varnames{ivar});
            end
        end
        
        function model = removeVarNames(model, varnames)
        % Remove several variable name using removeVarName
            for ivar = 1 : numel(varnames)
                model = model.removeVarName(varnames{ivar});
            end
        end
        
        function model = removeVarName(model, varname)
        % Remove a variable name. It will also remove the property functions where the variable name has been us either
        % in the input or output.
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
                    [found, keep(ivar), varnames{ivar}] = BaseModel.extractVarName(varname, varnames{ivar});
                    if found
                        break;
                    end
                end
                model.varNameList = varnames(keep);
                
                % remove from propertyFunctionList if varname occurs in output and input variables
                propfuncs = model.propertyFunctionList;
                nprops = numel(propfuncs);
                keep = true(nprops, 1);
                for iprop = 1 : nprops
                    [~, keep(iprop),  propfuncs{iprop}.varname] = BaseModel.extractVarName(varname, propfuncs{iprop}.varname);
                    inputvarnames = propfuncs{iprop}.inputvarnames;
                    for iinput = 1 : numel(inputvarnames)
                        [~, keepprop,  inputvarnames{iinput}] = BaseModel.extractVarName(varname, inputvarnames{iinput});
                        keep(iprop) = keepprop & keep(iprop);
                    end
                    propfunctions{iprop}.inputvarnames = inputvarnames;
                end
                
                model.propertyFunctionList = propfuncs(keep);                

            else
                error('varname not recognized');
            end
            
        end

        
        function model = unsetAsExtraVarName(model, varname)

            model.extraVarNameList = BaseModel.removeVarNameFromList(varname, model.extraVarNameList);

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
                if iscell(fn)
                    functionCallSetupFn = fn{2};
                    fn = fn{1};
                else
                    functionCallSetupFn = [];
                end
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

                propfunc = PropFunction(varname, fn, inputvarnames, modelnamespace, functionCallSetupFn);

                model = model.registerPropFunction(propfunc);
            else
                error('propfunc not recognized');
            end                
        end

                
        function stateAD = initStateAD(model, state)
        % initialize a new cleaned-up state with AD variables
            
            pnames  = model.getPrimaryVariableNames();
            vars = cell(numel(pnames), 1);
            for i = 1:numel(pnames)
                vars{i} = model.getProp(state, pnames{i});
            end
            % Get the AD state for this model           
            [vars{:}] = model.AutoDiffBackend.initVariablesAD(vars{:});
            stateAD = struct();
            for i = 1:numel(pnames)
               stateAD = model.setNewProp(stateAD, pnames{i}, vars{i});
            end

            stateAD = model.addStaticVariables(stateAD, state);
            
        end 
        
        function submodelnames = getSubModelNames(model)

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

            end

        end
        
        function model = registerSubModels(model)

            propfuncs     = model.propertyFunctionList;
            varnames      = model.varNameList;
            extravarnames = model.extraVarNameList;
            
            submodelnames = model.getSubModelNames();
            
            for isub = 1 : numel(submodelnames)

                submodelname = submodelnames{isub};
                submodel = model.(submodelname);

                if isa(submodel, 'BaseModel')
                    
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

                    % Register the static variables
                    subextravarnames = submodel.extraVarNameList;
                    
                    for isubvar = 1 : numel(subextravarnames)
                        subaddvarname = subextravarnames{isubvar};
                        subaddvarname.namespace = {submodelname, subaddvarname.namespace{:}};
                        subextravarnames{isubvar} = subaddvarname;
                    end
                    
                    extravarnames = mergeList(extravarnames, subextravarnames);

                end
                
            end
            
             model.varNameList          = varnames;
             model.propertyFunctionList = propfuncs;
             model.extraVarNameList     = extravarnames;
             
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
                if not(isfield(state,name))
                    state.(name) = struct();
                end
                state.(name) = model.setNewProp(state.(name), names, ...
                                                               val);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    if not(iscell(state))
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

            p = model.getPrimaryVariableNames();

            for i = 1 : numel(dx)
                val = model.getProp(state, p{i});
                val = val + dx{i};
                state = model.setProp(state, p{i}, val);
            end
           
            report = [];
            
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

            [state, report] = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
            p = model.getPrimaryVariableNames();
            cleanState = [];
            for ind = 1 : numel(p)
                cleanState = model.copyProp(cleanState, state, p{ind});
            end

            cleanState = model.addStaticVariables(cleanState, state);
            cleanState = model.addVariablesAfterConvergence(cleanState, state);
            
            state = cleanState;
            report = [];
            
        end

        function newstate = addVariablesAfterConvergence(model, newstate, state)
        % Function called in updateAfterConvergence
            
            submodelnames = model.getSubModelNames();
            
            for isub = 1 : numel(submodelnames)

                submodelname = submodelnames{isub};

                if isfield(state, submodelname)
                    
                    if ~isfield(newstate, submodelname)
                        newstate.(submodelname) = [];
                    end

                    if isa(model.(submodelname), 'BaseModel')
                        newstate.(submodelname) = model.(submodelname).addVariablesAfterConvergence(newstate.(submodelname), state.(submodelname));
                    end
                    
                end

            end
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state)
        % function to add static variables (not AD) on the cleanState, called in updateAfterConvergence and
        % initStateAD. Time is added by default here (when it exists)
            if isfield(state, 'time')
                cleanState.time = state.time;
            end

            submodelnames = model.getSubModelNames();
            
            for isub = 1 : numel(submodelnames)

                submodelname = submodelnames{isub};

                if isfield(state, submodelname)
                    
                    if ~isfield(cleanState, submodelname)
                        cleanState.(submodelname) = [];
                    end

                    if isa(model.(submodelname), 'BaseModel')
                        cleanState.(submodelname) = model.(submodelname).addStaticVariables(cleanState.(submodelname), state.(submodelname));
                    end 
                end

            end
            
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
        
        function state = reduceState(model, state, removeContainers)
            state = value(state, false);
        end
        
        function outputvars = extractGlobalVariables(model, states)
            ns = numel(states);
            outputvars = cell(1, ns);
        end
        
        function state = evalVarName(model, state, varname, extravars)
        % varname is valid input to method ComputationalGraph.getPropFunctionCallList

            cg = model.computationalGraph;

            if isempty(cg)
                fprintf('The computational graph has not been set up in model so that we compute set it up now, but it is a time coslty operation\n')
                cg = ComputationalGraph(model);
            end

            if iscell(varname)
                varname = VarName(varname(1 : end - 1), varname{end});
            end
            
            funcCallList = cg.getPropFunctionCallList(varname);

            if (nargin > 3) && ~isempty(extravars)
                for ivar = 1 : numel(extravars)
                    extravar = extravars{ivar};
                    varname  = extravar{1};
                    varvalue = extravar{2};
                    eval(sprintf('%s = varvalue;', varname));
                end
            end
            
            funcCall = strjoin(funcCallList, '');
            eval(funcCall);
            
        end

        %% Methods used when the model is used as root model for a simulation. Then, the model is equipped for simulation
        %
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            opts = struct('ResOnly', false, 'iteration', 0, 'reverseMode', false);
            opts = merge_options(opts, varargin{:});
            
            if (not(opts.ResOnly) && not(opts.reverseMode))
                state = model.initStateAD(state);
            elseif opts.reverseMode
                dispif(mrstVerbose, 'No AD initialization in equation old style')
                state0 = model.initStateAD(state0);
            else
                assert(opts.ResOnly);
            end
            
            %% We call the assembly equations ordered from the graph

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            state = model.applyScaling(state);
            
            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end
            
            names       = model.equationNames;
            types       = model.equationTypes;
            primaryVars = model.primaryVarNames;
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function primaryvarnames = getPrimaryVariableNames(model)

            primaryvarnames = model.primaryVarNames;

        end
        
        function state = applyScaling(model, state)
            
            if ~isempty(model.scalings)

                scalings = model.scalings;

                for iscal = 1 : numel(scalings)

                    scaling = scalings{iscal};
                    name = scaling{1};
                    coef = scaling{2};

                    val = model.getProp(state, name);
                    val = 1./coef.*val;

                    state = model.setProp(state, name, val);
                    
                end
                
            end
            
        end
                
        function model = equipModelForComputation(model, varargin)

            opt = struct('shortNames', []);
            opt = merge_options(opt, varargin{:});

            model = model.setupComputationalGraph();

            cg = model.computationalGraph();
            
            model.funcCallList     = cg.getOrderedFunctionCallList();
            model.primaryVarNames  = cg.getPrimaryVariableNames();
            model.equationVarNames = cg.getEquationVariableNames();

            function str = shortenName(name)
                [found, ind] = ismember(name, opt.shortNames(:, 1));
                if found
                    str = opt.shortNames{ind, 2};
                else
                    str = name;
                end
            end
            
            function str = setupName(varname)

                if isnumeric(varname{end})
                    varname{end} = sprintf('%d', varname{end});
                end

                if ~isempty(opt.shortNames)
                    varname = cellfun(@(elt) shortenName(elt), varname, 'uniformoutput', false);
                end

                str = strjoin(varname, '_');
                
            end
            
            model.equationNames = cellfun(@(varname) setupName(varname), model.equationVarNames, 'uniformoutput', false);
            model.equationTypes = repmat({'cell'}, 1, numel(model.equationNames));
            
        end

        function jsonstruct = exportParams(model)

            submodelnames = model.getSubModelNames();

            if numel(submodelnames) > 0
                for isubmodel = 1 : numel(submodelnames)

                    submodelname = submodelnames{isubmodel};

                    jsonstruct.(submodelname) = model.(submodelname).exportParams();
                    
                end
            else
                jsonstruct = [];
            end
        end

        function cgit = getComputationalGrapInteractiveTool(model)
        % setup and retrieve the computational graph interactive tool
            if isempty(model.computationalGraph)
                model = model.setupComputationalGraph();
            end
            cgit = ComputationalGraphInteractiveTool(model.computationalGraph);

        end

        function cgit = cgit(model)
            cgit = model.getComputationalGrapInteractiveTool();
        end

        function G = grid(model)
        % Shorcut to retrieve grid in defaut MRST format, which can be used for plotting

            if isa(model.G, 'GenericGrid')
                G = model.G.mrstFormat();
            else
                G = model.G;
            end
            
        end

        function model = validateModel(model, varargin)

        % By default, discard validateModel from MRST
            
        end

        function inds = getRangePrimaryVariable(model, adsample, varname)

            ivar = model.getIndexPrimaryVariable(varname);

            ss = cellfun(@(jac) size(jac, 2), adsample.jac);
            ss = cumsum(ss);
            ss = [1, ss(1 : end - 1) + 1; ...
                  ss];
            inds = ss(:, ivar);

        end
        
        function ind = getIndexPrimaryVariable(model, varname)

            primvarnames = model.getPrimaryVariableNames();

            for ivar = 1 : numel(primvarnames)
                isequal = ImpedanceBattery.compareVarName(varname, primvarnames{ivar});
                if isequal
                    ind = ivar;
                    return
                end
            end

            error('primary variable not found');
            
        end
                
    end

    methods(Static)

        function [found, keep, varname] = extractVarName(rvarname, varname)

        % This function supports index, meaning, it handles the case where the rvarname and varname has same name but
        % different indices where the indices of rvarname make a subset of those of varname.
            
            [isequal, compIndices] = compareVarName(rvarname, varname);
            if isequal
                keep  = false;
                found = true;
            elseif isempty(compIndices)
                keep  = true;
                found = false;
            elseif ~isempty(compIndices.InterInd)
                found = true;
                if isempty(compIndices.OuterInd2)
                    keep = false;
                else
                    keep = true;
                    varname.index = compIndices.OuterInd2;
                end
            else
                found = false;
                keep  = true;
            end
        end
        

        function varnames = removeVarNameFromList(varname, varnames)

            if isa(varname, 'char')
                varname = VarName({}, varname);
                varnames = BaseModel.removeVarNameFromList(varname, varnames);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                varnames = BaseModel.removeVarNameFromList(varname, varnames);
            elseif isa(varname, 'VarName')
                % remove from varnames
                nvars = numel(varnames);
                keep = true(nvars, 1);
                for ivar = 1 : nvars
                    [found, keep(ivar), varnames{ivar}] = BaseModel.extractVarName(varname, varnames{ivar});
                    if found
                        break;
                    end
                end
                varnames = varnames(keep);
            else
                error('varname not recognized');
            end
            
        end
        
    end
    
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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

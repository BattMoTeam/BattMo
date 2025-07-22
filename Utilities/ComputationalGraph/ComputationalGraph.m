classdef ComputationalGraph
%
% A model is essentially characterised by its computational graph. The functions to setup the graph of a model are given
% in the BaseModel class.
%
    
    properties (SetAccess = private)

        model

    end

    properties
        
        adjencyMatrix    % Adjency matrix for the computational graph
                         % - column index      : output variable index (as in cgt.varNameList)
                         % - row index         : input variable (as in cgt.varNameList)
                         % - coefficient value : property function index as in model.propertyFunctionList
        varNameList      % List of variable names (object of class VarName).
                         % They have been processed model.varNameList so that there exist separate entries for each index in the case of variable names with indices.
        nodenames        % The node names in the graph are unique string representations of the variable names in varNameList.
        staticprops      % list of the static properties (property functions that take no input). It is computed in setupComputationalGraph.
                         % Each element is a struct with the  the fields
                         % - nodename (element in nodenames)
                         % - propind (index of the property function in the list model.propertyFunctionList)
                         % - varnameind (index of the variable name in varNameList)
        extraVarNameInds % index relative with respect to varNameList
        isok             % The graph is a proper directed acyclic graph that can be used for computation.

        %% model graph (for visualization only, i.e. not used in assembly)
        modelnames         % List of strings given the model names (contructed using getSubModelNames from BaseModel)
        modelAdjencyMatrix % model adjency matrix

    end

    methods

        function cg = ComputationalGraph(model)

            if isempty(model.propertyFunctionList)
                model = model.registerVarAndPropfuncNames();
            end
            
            cg.model = model;
            [A, staticprops, varNameList, nodenames] = setupGraph(model);

            % In adjacency matrix A,
            % - column index      : output variable index (as in cg.varNameList)
            % - row index         : input variable (as in cg.varNameList)
            % - coefficient value : property function index as in model.propertyFunctionList

            cg.adjencyMatrix = A;
            cg.varNameList   = varNameList;
            cg.nodenames     = nodenames;
            cg.staticprops   = staticprops;
            cg.isok        = false;
            if size(A, 1) == numel(varNameList)
                cg = cg.setupComputationalGraph();
            else
                fprintf('\nThe graph could not be ordered properly due a mismatch in the variable declarations\nFix that before using the graph in computations\n');
            end

        end

        function cg = setupComputationalGraph(cg)

            A            = cg.adjencyMatrix;
            varNameList  = cg.varNameList;
            nodenames    = cg.nodenames;
            staticprops  = cg.staticprops;

            p = compatible_topological_order(A);

            if isempty(p)
                fprintf('The graph contains cycles. It implies that some variables cannot be evaluated.\n')
                cg.isok = false;
                return
            end

            varNameList = varNameList(p);
            A           = A(p, p);
            nodenames   = nodenames(p);

            for istat = 1 : numel(staticprops)
                nodename = staticprops{istat}.nodename;
                [isok, varnameind] = ismember(nodename, nodenames);
                if ~isok
                    fprintf('The static variable %s has not been found in list\n', nodename);
                    cg.isok = false;
                    return
                end

                staticprops{istat}.varnameind = varnameind;
            end

            cg.varNameList   = varNameList;
            cg.adjencyMatrix = A;
            cg.nodenames     = nodenames;
            cg.staticprops   = staticprops;

            extraVarNames = cg.model.extraVarNameList;
            extraVarNameInds = [];
            for ivar = 1 : numel(extraVarNames);
                extraVarName = extraVarNames{ivar};
                extraVarName_s = extraVarName.resolveIndex();
                for iivar = 1 : numel(extraVarName_s)
                    extraVarNameInds(end + 1) = cg.getVarNameIndex(extraVarName_s{iivar});
                end
            end
            cg.extraVarNameInds = extraVarNameInds;

            cg.isok = true;
        end

        function [varnames, varnameinds, propfuncinds, levels] = getDependencyList(cg, varname, direction)
        % dependency direction can be
        % - 'upwards' for upwards in the graph
        % - 'downwards' for downwards in the graph
            
            A = cg.adjencyMatrix;

            % for cell-valued variable pick-up a valid index
            if (varname.dim > 1) && (ischar(varname.index))
                varname.index = 1;
            elseif isnumeric(varname.index) && (numel(varname.index) > 1)
                varname.index = varname.index(1);
            end

            varnameind = cg.getVarNameIndex(varname);

            switch direction
              case 'upwards'
                [depvarnameinds, propfuncinds, propvarnameinds, propdeplevels, rootinds, rootdeplevels] = getDependencyVarNameInds(varnameind, A);
              case 'downwards'
                [depvarnameinds, propfuncinds, propvarnameinds, propdeplevels, rootinds, rootdeplevels] = getDependencyVarNameInds(varnameind, A');
              otherwise
                error('direction not recognized')
            end
            
            varnameinds = depvarnameinds;
            levels      = [rootdeplevels; propdeplevels];
            varnames    = cg.varNameList(varnameinds);

        end

        function varnameind = findVarName(cg, varname)

            nodenames = cg.nodenames;

            if isa(varname, 'VarName')
                nodename = ComputationalGraph.getNodeName(varname);
                varnameind = findVarName(nodename);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                varnameind = findVarName(varname);
            elseif isa(varname, 'char')
                varnameinds = regexpSelect(cg.nodenames, varname);
                if numel(varnameinds) > 1
                    fprintf('Several variables found:\n\n')
                    for ivar = 1 : numel(varnameinds)
                        fprintf('%d : %s\n', ivar, nodenames{varnameinds(ivar)});
                    end
                    ivar = input('Pick one:');
                    varnameind = varnameinds(ivar);
                else
                    varnameind = varnameinds;
                end
            else
                error('input varname not recognized');
            end

        end
        
        function [propfuncs, propfuncinds, varnameinds] = getPropFunctionList(cg, propfunc)

        % Get the list of property functions and the corresponding indices (with respect to
        % cg.model.propertyFunctionList) that should be call so that the variable of propfunc get updated.

            A           = cg.adjencyMatrix;
            staticprops = cg.staticprops;

            varname = propfunc.varname;

            % for cell-valued variable pick-up a valid index
            if (varname.dim > 1) && (ischar(varname.index))
                varname.index = 1;
            elseif isnumeric(varname.index) && (numel(varname.index) > 1)
                varname.index = varname.index(1);
            end

            varnameind = cg.getVarNameIndex(varname);

            [~, propfuncinds, propvarnameinds, ~, rootinds, ~] = getDependencyVarNameInds(varnameind, A);

            [rootinds, rootpropinds] = cg.findStaticVarNameInds(rootinds);

            if ~isempty(rootinds)
                propfuncinds = [rootpropinds; propfuncinds];
                varnameinds  = [rootinds; propvarnameinds];
            else
                varnameinds = propvarnameinds;
            end

            propfuncs = cg.model.propertyFunctionList(propfuncinds);

        end

        function [staticinds, staticpropinds] = findStaticVarNameInds(cg, varnameinds)

        % Given a set list of index of varnameinds (index with respect to cg.varNameList), extract and return the indices that
        % correspond to static variables.

            staticprops = cg.staticprops;

            allstaticinds     = cellfun(@(staticprop) staticprop.varnameind, staticprops)';
            allstaticpropinds = cellfun(@(staticprop) staticprop.propind, staticprops)';

            [isok, ia] = ismember(varnameinds, allstaticinds);

            if any(isok)
                staticpropinds = allstaticpropinds(ia(isok), 1);
                staticinds     = allstaticinds(ia(isok), 1);
            else
                staticpropinds = {};
                staticinds = {};
            end

        end

        function funcCallList = getPropFunctionCallList(cg, propfunc)
        % input propfunc is either
        % - an instance of PropFunction
        % - a valid input for findPropFunction, that is, either
        %     - a VarName instance
        %     - a cell which then uses shortcuts for VarName (see implementation in findPropFunction)
        %     - a string giving a regexp. It will be used to select varnames by the string name
        %   In this case, findPropFunction is first run to obtain a list of propfunctions
        % setup the list of function call (as list of cell of strings) that will be run to update the property function propfunc.

            if isa(propfunc, 'PropFunction')
                propfuncs = cg.getPropFunctionList(propfunc);
            elseif isa(propfunc, 'VarName')
                % We set up a propfunction with only the varname
                propfunc = PropFunction(propfunc, [], [], []);
                propfuncs = cg.getPropFunctionList(propfunc);
            else

                varname = propfunc;
                selectedpropfuncs = cg.findPropFunction(propfunc);
                if isa(selectedpropfuncs, 'PropFunction')
                    selectedpropfuncs = {selectedpropfuncs};
                end
                propfuncs   = {};
                for isel = 1 : numel(selectedpropfuncs)
                    deppropfuncs = cg.getPropFunctionList(selectedpropfuncs{isel});
                    propfuncs   = horzcat(propfuncs, deppropfuncs);
                end
            end

            funcCallList = {};

            for iprop = 1 : numel(propfuncs)

                propfunc = propfuncs{iprop};
                fncallstr = propfunc.functionCallSetupFn(propfunc);
                if ~isempty(fncallstr)
                    funcCallList{end + 1} = fncallstr;
                end

            end

        end


        function nodestaticnames = getStaticVarNames(cg)

        % add the variables that are updated with no input arguments
            propfunctions = cg.model.propertyFunctionList;
            varnames = {};
            for iprop = 1 : numel(propfunctions)
                propfunction = propfunctions{iprop};
                if isempty(propfunction.inputvarnames)
                    varnames{end + 1} = propfunction.varname;
                end
            end

            nodestaticnames = cg.getNodeNames(varnames);

        end

        function [propfuncs, propfuncinds] = findPropFunction(cg, varname)
        % The input varnames can be either:
        % - a VarName instance
        % - a cell which then uses shortcuts for VarName (see implementation below)
        % - a string giving a regexp. It will be used to select varnames by the string name
        %
        % The function returns a list of PropFunction that updates the variable given by varname.

            nodenames = cg.nodenames;
            A         = cg.adjencyMatrix;
            model     = cg.model;

            if isa(varname, 'VarName')
                varnameinds = cg.getVarNameIndex(varname);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                [propfuncs, propfuncinds] = cg.findPropFunction(varname);
                return
            elseif isa(varname, 'char')
                selectednodenames = varname;
                varnameinds = regexpSelect(cg.nodenames, varname);
            else
                error('input type not recognized');
            end

            propfuncinds = max(A(:, varnameinds), [], 1)';
            staticinds = varnameinds(propfuncinds == 0);
            propfuncinds = propfuncinds(propfuncinds > 0); % remove the zero elements

            [staticinds, staticpropinds] = cg.findStaticVarNameInds(staticinds);

            if ~isempty(staticinds)
                propfuncinds = [staticpropinds; propfuncinds];
            end

            propfuncinds = unique(propfuncinds);
            
            propfuncs = model.propertyFunctionList(propfuncinds);

            if numel(propfuncs) == 1
                propfuncs = propfuncs{1};
            end

        end

        function varnameind = getVarNameIndex(cg, varname)
        % Given a varname, find its index in cg.varNameList
            for varnameind = 1 : numel(cg.varNameList)
                if varname.compareVarName(cg.varNameList{varnameind})
                    return
                end
            end

            fprintf('Variable not found\n');
            varnameind = [];


        end

        function cg = setupModelGraph(cg, varargin)


            function g = addModel(g, model, modelname, namespace, parentnodename)

                namecomponents = horzcat(namespace, {modelname});
                nodename = strjoin(namecomponents, '.');

                g = addnode(g, nodename);

                if ~isempty(parentnodename)
                    g = addedge(g, parentnodename, nodename);
                end
                submodelnames = model.getSubModelNames();
                namespace = horzcat(namespace, modelname);

                for isubmodel = 1 : numel(submodelnames)

                    submodelname = submodelnames{isubmodel};
                    submodel = model.(submodelname);

                    g = addModel(g, submodel, submodelname, namespace, nodename);

                end

            end

            g = digraph();
            g = addModel(g, cg.model, 'model', {}, []);

            cg.modelnames = g.Nodes.Variables;
            cg.modelAdjencyMatrix = adjacency(g, 'weighted');

        end

        function primvarnames = getPrimaryVariableNames(cg)
        % Return the primary variables, which are defined as the root variables and not declared or recognized as static.

            A = cg.adjencyMatrix;
            nodenames = cg.nodenames;

            ind = find(all(A == 0, 1));

            rootnames = nodenames(ind);
            staticnames = cg.getStaticVarNames();
            indstatic = ismember(rootnames, staticnames);
            ind = ind(~indstatic);
            primVarNameList = cg.varNameList(ind);

            primvarnames = cellfun(@(varname) varname.getPropName(), primVarNameList, 'uniformoutput', false);

        end

        function eqvarnames = getEquationVariableNames(cg)
        % Return the equation variables, which are defined as the tail variables and not declared as extravarnames.

            A                = cg.adjencyMatrix;
            nodenames        = cg.nodenames;
            extravarnameinds = cg.extraVarNameInds;

            eqinds = find(all(A' == 0, 1));
            isadded  = ismember(eqinds, extravarnameinds);
            eqinds = eqinds(~isadded);

            equationVarNameList = cg.varNameList(eqinds);

            eqvarnames = cellfun(@(varname) varname.getPropName(), equationVarNameList, 'uniformoutput', false);

        end


        function funcCallList = getOrderedFunctionCallList(cg, varargin)
        % Return a list of strings that corresponds to all the function calls ordered in the right order.

            opt = struct('removeExtraVariables', true);
            opt = merge_options(opt, varargin{:});

            A                = cg.adjencyMatrix;
            staticprops      = cg.staticprops;
            extraVarNameInds = cg.extraVarNameInds;

            funcCallList = {};

            if ~isempty(staticprops)

                for istatic = 1 : numel(staticprops)

                    iprop = staticprops{istatic}.propind;
                    propfunction = cg.model.propertyFunctionList{iprop};
                    fncallstr = propfunction.getFunctionCallString();
                    if ~isempty(fncallstr)
                        funcCallList{end + 1} = fncallstr;
                    end

                end

            end

            if opt.removeExtraVariables
                nA = size(A, 2);
                ind = true(nA, 1);
                ind(extraVarNameInds) = false;
                A = A(ind, ind);
            end

            for ind = 1 : size(A, 2)
                iprop = full(A(:, ind));
                iprop = unique(iprop(iprop>0));
                if ~isempty(iprop)
                    assert(numel(iprop) == 1, 'There should be only one value for a given row');
                    propfunction = cg.model.propertyFunctionList{iprop};
                    fncallstr = propfunction.getFunctionCallString();
                    if ~isempty(fncallstr)
                        funcCallList{end + 1} = fncallstr;
                    end
                end
            end

            % we remove duplicates (the same function can be used for several output variables)
            [~, ia, ic] = unique(funcCallList, 'first');

            ia = sort(ia); % Not sure if it is necessary, but in this way we make sure the ordering for which  funcCallList was built is respected
            funcCallList = funcCallList(ia);

        end
        
    end

    methods (Static)

        function nodenames = getNodeName(varname)
        % convert variable names to graph node names
        % TODO : we should enforce that the same function is used in setupGraph (important!)
            nodenames = {};
            varname_s = varname.resolveIndex();
            for ind = 1 : numel(varname_s)
                nodenames{end + 1} = varname_s{ind}.getIndexedFieldname();
            end

        end

        function nodenames = getNodeNames(varnames)
        % convert variable names to graph node names
        % TODO : we should enforce that the same function is used in setupGraph (important!)
            nodenames = {};
            for ind = 1 : numel(varnames)
                varname = varnames{ind};
                nodenames = horzcat(nodenames, ...
                                    ComputationalGraph.getNodeName(varname));
            end

        end

        function printHeader(headertxt, n)
        % Minor utility function used in this class to print a header with a number
            str = sprintf('\n%d %s', n, headertxt);
            fprintf('%s:\n', str);
            fprintf('%s\n', repmat('-', length(str), 1));
        end
    end

end




%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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

classdef ComputationalGraphTool

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
        
        function cgt = ComputationalGraphTool(model)
            if isempty(model.propertyFunctionList)
                model = model.registerVarAndPropfuncNames();
            end
            cgt.model = model;
            [A, staticprops, varNameList, nodenames] = setupGraph(model);
            
            % In adjacency matrix A, 
            % - column index      : output variable index (as in cgt.varNameList)
            % - row index         : input variable (as in cgt.varNameList)
            % - coefficient value : property function index as in model.propertyFunctionList

            cgt.adjencyMatrix = A;
            cgt.varNameList   = varNameList;
            cgt.nodenames     = nodenames;
            cgt.staticprops   = staticprops;
            cgt.isok        = false;
            if size(A, 1) == numel(varNameList)
                try
                    require('matlab_bgl');
                catch
                    mrstModule add matlab_bgl
                end
                cgt = cgt.setupComputationalGraph();
            else
                fprintf('\nThe graph could not be ordered properly due a mismatch in the variable declarations\nFix that before using the graph in computations\n');
            end

        end

        function cgt = setupComputationalGraph(cgt)

            A            = cgt.adjencyMatrix;
            varNameList  = cgt.varNameList;
            nodenames    = cgt.nodenames;
            staticprops  = cgt.staticprops;
            
            p = compatible_topological_order(A);

            if isempty(p)
                fprintf('The graph contains cycles. It implies that some variables cannot be evaluated.\n')
                cgt.isok = false;
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
                    cgt.isok = false;
                    return
                end
                
                staticprops{istat}.varnameind = varnameind;
            end

            cgt.varNameList   = varNameList;
            cgt.adjencyMatrix = A;
            cgt.nodenames     = nodenames;
            cgt.staticprops   = staticprops;

            extraVarNames = cgt.model.extraVarNameList;
            extraVarNameInds = [];
            for ivar = 1 : numel(extraVarNames);
                extraVarNameInds(end + 1) = cgt.getVarNameIndex(extraVarNames{ivar});
            end
            cgt.extraVarNameInds = extraVarNameInds;

            cgt.isok = true;
        end

        function [varnames, varnameinds, propfuncinds, distance] = getDependencyList(cgt, varname)

            A = cgt.adjencyMatrix;
            
            % for cell-valued variable pick-up a valid index
            if (varname.dim > 1) && (ischar(varname.index))
                varname.index = 1;
            elseif isnumeric(varname.index) && (numel(varname.index) > 1)
                varname.index = varname.index(1);
            end
            
            varnameind = cgt.getVarNameIndex(varname);

            distance = dfs(A, varnameind);
            
            varnameinds = find(distance >= 0);
            distance    = distance(distance >= 0);
            [distance, ia] = sort(distance);
            varnameinds = varnameinds(ia);

            propfuncinds = max(A(:, varnameinds), [], 1);

            varnames = cgt.varNameList(varnameinds);
            
        end
        
        function printDependencyList(cgt, varname)
        % input varname is either
        %     - a VarName instance
        %     - a cell which then uses shortcuts for VarName (see implementation below)
        %     - a string giving a regexp. It will be used to select varnames by the string name

            nodenames = cgt.nodenames;

            if isa(varname, 'VarName')
                % ok
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
            elseif isa(varname, 'char')
                varnameind = regexpSelect(cgt.nodenames, varname);
                if numel(varnameind) > 1
                    fprintf('Several variables found:\n\n')
                    for ivar = 1 : numel(varnameind)
                        fprintf('%s\n', nodenames{varnameind(ivar)});
                    end
                    fprintf('\npick a single one\n')
                    return
                end
                varname = cgt.varNameList{varnameind};
            else
                error('input varname not recognized');
            end

            [varnames, varnameinds, propfuncinds, distance] = cgt.getDependencyList(varname);

            for ivar = 1 : numel(varnameinds)
                varnameind = varnameinds(ivar);
                fprintf('%s (%d)\n', nodenames{varnameind}, distance(ivar));
            end

        end

        function [propfuncs, propfuncinds, varnameinds] = getPropFunctionList(cgt, propfunc)

        % Get the list of property functions and the corresponding indices (with respect to
        % cgt.model.propertyFunctionList) that should be call so that the variable of propfunc get updated.
            
            A           = cgt.adjencyMatrix;
            staticprops = cgt.staticprops;
            
            varname = propfunc.varname;

            % for cell-valued variable pick-up a valid index
            if (varname.dim > 1) && (ischar(varname.index))
                varname.index = 1;
            elseif isnumeric(varname.index) && (numel(varname.index) > 1)
                varname.index = varname.index(1);
            end
            
            varnameind = cgt.getVarNameIndex(varname);

            [~, propfuncinds, propvarnameinds, staticinds] = getDependencyVarNameInds2(varnameind, A);
            
            [staticinds, staticpropinds] = cgt.findStaticVarNameInds(staticinds);

            if ~isempty(staticinds)
                propfuncinds = [staticpropinds; propfuncinds];
                varnameinds  = [staticinds; propvarnameinds];
            end

            propfuncs = cgt.model.propertyFunctionList(propfuncinds);
            
        end

        function [staticinds, staticpropinds] = findStaticVarNameInds(cgt, varnameinds)

        % Given a set list of index of varnameinds (index with respect to cgt.varNameList), extract and return the indices that
        % correspond to static variables.

            staticprops = cgt.staticprops;
           
            allstaticinds     = cellfun(@(staticprop) staticprop.varnameind, staticprops);
            allstaticpropinds = cellfun(@(staticprop) staticprop.propind, staticprops);
            
            [isok, ia] = ismember(varnameinds, allstaticinds);
            
            if any(isok)
                staticpropinds = allstaticpropinds(ia(isok));
                staticinds = allstaticinds(ia(isok))'; %
            else
                staticpropinds = {};
                staticinds = {};
            end
            
        end
        
        function funcCallList = getPropFunctionCallList(cgt, propfunc)
        % input propfunc is either
        % - an instance of PropFunction
        % - a valid input for findPropFunction, that is, either
        %     - a VarName instance
        %     - a cell which then uses shortcuts for VarName (see implementation in findPropFunction)
        %     - a string giving a regexp. It will be used to select varnames by the string name
        %   In this case, findPropFunction is first run to obtain a list of propfunctions
        % setup the list of function call (as list of cell of strings) that will be run to update the property function propfunc.

            if isa(propfunc, 'PropFunction')
                propfuncs = cgt.getPropFunctionList(propfunc);
            elseif isa(propfunc, 'VarName')
                % We set up a propfunction with only the varname
                propfunc = PropFunction(propfunc, [], [], []);
                propfuncs = cgt.getPropFunctionList(propfunc);                
            else
                
                varname = propfunc;
                selectedpropfuncs = cgt.findPropFunction(propfunc);
                if isa(selectedpropfuncs, 'PropFunction')
                    selectedpropfuncs = {selectedpropfuncs};
                end
                propfuncs   = {};
                for isel = 1 : numel(selectedpropfuncs)
                    deppropfuncs = cgt.getPropFunctionList(selectedpropfuncs{isel});
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

        function printPropFunctionCallList(cgt, propfunc, varargin)
        % input propfunc is either
        % - an instance of PropFunction
        % - a valid input for findPropFunction, that is, either
        %     - a VarName instance
        %     - a cell which then uses shortcuts for VarName (see implementation below)
        %     - a string giving a regexp. It will be used to select varnames by the string name
        %   In this case, findPropFunction is first run to obtain the property function. It should return a unique element, otherwise we get a warning with a list of returned elements
        % Print the list of function call (as string) that will update the property function propfunc.


            opt = struct('fullSignature', false);
            opt = merge_options(opt, varargin{:});

            if ~isa(propfunc, 'PropFunction')

                varname = propfunc;

                propfunc = cgt.findPropFunction(varname);
                
                if isempty(propfunc)
                    fprintf('No property matching regexp has been found\n');
                    return
                end            
                
                if numel(propfunc) > 1
                    fprintf('Several property functions are matching\n\n');
                    cgt.printPropFunction(varname);
                    return
                end
                
            end
            
            nodenames = cgt.nodenames;
            
            [propfuncs, ~, varnameinds] = cgt.getPropFunctionList(propfunc);
            
            funcCallList = {};

            strls = arrayfun(@(varnameind) strlength(nodenames{varnameind}), varnameinds);
            strl = max(strls);
            
            for iprop = 1 : numel(propfuncs)

                propfunc = propfuncs{iprop};
                varnameind = varnameinds(iprop);

                if opt.fullSignature
                    fncallstr = propfunc.functionCallSetupFn(propfunc);
                    varstr = sprintf('state.%s', nodenames{varnameind});
                    varstr = sprintf('%-*s', strl + 6, varstr);
                    fprintf('%s <- %s\n', varstr, fncallstr);
                else
                    inputvarnames = propfunc.inputvarnames;
                    strs = {};
                    for ivar = 1 : numel(inputvarnames)
                        inputvarname = inputvarnames{ivar};
                        strs = horzcat(strs, cgt.getNodeName(inputvarname));
                    end
                    strs = strjoin(strs, ', ');
                    varstr = sprintf('%-*s', strl + 6, nodenames{varnameind});
                    fprintf('%s [%s]\n', varstr, strs);
                end
                
            end

        end

        function nodestaticnames = getStaticVarNames(cgt)

            % add the variables that are updated with no input arguments
            propfunctions = cgt.model.propertyFunctionList;
            varnames = {};
            for iprop = 1 : numel(propfunctions)
                propfunction = propfunctions{iprop};
                if isempty(propfunction.inputvarnames)
                    varnames{end + 1} = propfunction.varname;
                end
            end
            
            nodestaticnames = cgt.getNodeNames(varnames);
            
        end
        
        function [propfuncs, propfuncinds] = findPropFunction(cgt, varname)
        % The input varnames can be either:
        % - a VarName instance
        % - a cell which then uses shortcuts for VarName (see implementation below)
        % - a string giving a regexp. It will be used to select varnames by the string name
        %
        % The function returns a list of PropFunction that updates the variable given by varname.
            
            nodenames = cgt.nodenames;
            A         = cgt.adjencyMatrix;
            model     = cgt.model;

            if isa(varname, 'VarName')
                varnameinds = cgt.getVarNameIndex(varname);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                propfuncs = cgt.findPropFunction(varname);
                return
            elseif isa(varname, 'char')
                selectednodenames = varname;
                varnameinds = regexpSelect(cgt.nodenames, varname);
            else
                error('input type not recognized');
            end
            
            propfuncinds = max(A(:, varnameinds), [], 1);
            staticinds = varnameinds(propfuncinds == 0);
            propfuncinds = propfuncinds(propfuncinds > 0); % remove the zero elements

            [staticinds, staticpropinds] = cgt.findStaticVarNameInds(staticinds);

            if ~isempty(staticinds)
                propfuncinds = [staticpropinds, propfuncinds];
            end

            propfuncs = model.propertyFunctionList(propfuncinds);

            if numel(propfuncs) == 1
                propfuncs = propfuncs{1};
            end
            
        end
        

        function openPropFunction(cgt, nodename)
        % Open in editor the place where the variable that matches the regexp nodename is updated. The regexp nodename
        % should return a unique match

            propfunc = cgt.findPropFunction(nodename);

            if isempty(propfunc)
                fprintf('No property matching regexp has been found\n');
                return
            end
            
            if numel(propfunc) > 1
                fprintf('Several property functions are matching\n\n');
                cgt.printPropFunction(nodename);
                return
            end

            fn = propfunc.fn;
            mn = propfunc.modelnamespace;
            mn = strjoin(mn, '.');
            fnname = func2str(fn);
            fnname = regexp(fnname, "\.(.*)", 'tokens');
            fnname = fnname{1}{1};
            fnname = sprintf('%s([])', fnname);
            fnname = strjoin({mn, fnname}, '.');
            fnname = sprintf('cgt.model.%s', fnname);

            try
                eval(fnname)
            catch ME
                fprintf('%s\n', fnname);
                stack = ME.stack;
                file = stack.file;
                lineNum = stack.line;
                editor = settings().matlab.editor.OtherEditor.ActiveValue;
                if contains(editor, 'emacs')
                    cmd = sprintf('%s +%d %s &', editor, lineNum, file);
                    system(cmd);
                else
                    matlab.desktop.editor.openAndGoToLine(file, lineNum);
                end
            end

        end

        function printVarNames(cgt, nodename)
        % Print the variable(s) that match the regexp nodename

            nodenames = cgt.nodenames;

            if nargin < 2
                nodename = '.';
            end
            inds = regexpSelect(cgt.nodenames, nodename);
            nodenames = nodenames(inds);
            for ind = 1 : numel(nodenames)
                fprintf('%s\n', nodenames{ind});
            end
            
        end

        
        function printPropFunction(cgt, nodename)
        % Print property function including output variable name, function name and input variable names for the variable(s)
        % that match the regexp nodename
            
            propfuncs = cgt.findPropFunction(nodename);

            if numel(propfuncs) == 1
                propfuncs = {propfuncs};
            end

            for iprop = 1 : numel(propfuncs)

                propfunc      = propfuncs{iprop};
                varname       = propfunc.varname;
                varname_s     = varname.resolveIndex();
                outputvarstrs = {};
                
                fncallstr = propfunc.getFunctionCallString();

                if ~isempty(fncallstr)
                    for ind = 1 : numel(varname_s)
                        fullname = varname_s{ind}.getIndexedFieldname();
                        outputvarstrs{end + 1} = sprintf('state.%s', fullname);
                    end
                    outputvarstr = strjoin(outputvarstrs, {', '});
                    fprintf('%s <-', outputvarstr);
                    fprintf(' (%s) ', fncallstr);
                    inputvarnames = propfunc.inputvarnames;
                    if isempty(inputvarnames)
                        fprintf('[no state field is used]\n');
                    else
                        fprintf(' <- ', fncallstr);
                        inputvarstrs = {};
                        for ind = 1 : numel(inputvarnames)
                            varname = inputvarnames{ind};
                            inputvarstrs{end + 1} = sprintf('state.%s', varname.getFieldname());
                        end
                        inputvarstr = strjoin(inputvarstrs, {', '});
                        fprintf('[%s]\n', inputvarstr);
                    end
                end
            end
            
        end

        function varnameind = getVarNameIndex(cgt, varname)
        % Given a varname, find its index in cgt.varNameList
            for varnameind = 1 : numel(cgt.varNameList)
                if varname.compareVarName(cgt.varNameList{varnameind})
                    return
                end
            end

            fprintf('Variable not found\n');
            varnameind = [];

            
        end

        function cgt = setupModelGraph(cgt, varargin)


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
            g = addModel(g, cgt.model, 'model', {}, []);

            cgt.modelnames = g.Nodes.Variables;
            cgt.modelAdjencyMatrix = adjacency(g, 'weighted');
            
        end

        function cgt2 = printModelNames(cgt, modelname)
            
            if isempty(cgt.modelnames)
                cgt = cgt.setupModelGraph();
                if nargout > 0
                    cgt2 = cgt;
                end
            end

            if nargin < 2
                modelname = '.';
            end
            
            modelnames = cgt.modelnames;
            inds = regexpSelect(modelnames, modelname);
            modelnames = modelnames(inds);

            for imodel = 1 : numel(modelnames)
                fprintf('%s\n', modelnames{imodel});
            end
            
        end

        
        function printRootVariables(cgt)
        % Print the root variables in computational graph 

            A = cgt.adjencyMatrix;
            nodenames = cgt.nodenames;
            
            %% print root variables after removing the variables that were declared as static in the model

            rootnames = nodenames(all(A == 0, 1));
            staticnames = cgt.getStaticVarNames();

            fprintf('Root variables \n');

            ind = ismember(rootnames, staticnames);
            rootnames(~ind)

            if nnz(ind)
                fprintf('static variables \n');
                rootnames(ind)
            end

        end

        function printTailVariables(cgt, nodename)
        % Print the tail variables of the computational graph
            
            A                = cgt.adjencyMatrix;
            nodenames        = cgt.nodenames;
            extravarnameinds = cgt.extraVarNameInds;
            
            tailinds = find(all(A' == 0, 1));
            isadded  = ismember(tailinds, extravarnameinds);
            tailinds = tailinds(~isadded);

            if nargin > 1
                inds = cgt.regexpVarNameSelect(nodename);
                tailinds = intersect(inds, tailinds);
            end
            
            fprintf('Tail variables \n');
            nodenames(tailinds)

            if nargin == 1
                fprintf('Extra Variables (do not enter in residual evaluation)\n')
                nodenames(extravarnameinds)
            end
            
        end        

        function printDetachedVariables(cgt)
        % Print the "detached" variables, which are the variables that are not connected to the graph. This is specially useful
        % in debugging because such variables should be eliminated from the final graph.
            
            A = cgt.adjencyMatrix;
            nodenames = cgt.nodenames;
            
            fprintf('Detached variables \n');
            ind1 = all(A == 0, 1);
            ind2 = all(A' == 0, 1);
            nodenames(ind1&ind2)
            
        end

        function primvarnames = getPrimaryVariableNames(cgt)
        % Return the primary variables, which are defined as the root variables and not declared or recognized as static.
            
            A = cgt.adjencyMatrix;
            nodenames = cgt.nodenames;

            ind = find(all(A == 0, 1));
            
            rootnames = nodenames(ind);
            staticnames = cgt.getStaticVarNames();
            indstatic = ismember(rootnames, staticnames);
            ind = ind(~indstatic);
            primVarNameList = cgt.varNameList(ind);
            
            primvarnames = cellfun(@(varname) varname.getPropName(), primVarNameList, 'uniformoutput', false);

        end

        function eqvarnames = getEquationVariableNames(cgt)
        % Return the equation variables, which are defined as the tail variables and not declared as extravarnames.
            
            A                = cgt.adjencyMatrix;
            nodenames        = cgt.nodenames;
            extravarnameinds = cgt.extraVarNameInds;

            eqinds = find(all(A' == 0, 1));
            isadded  = ismember(eqinds, extravarnameinds);
            eqinds = eqinds(~isadded);

            equationVarNameList = cgt.varNameList(eqinds);
            
            eqvarnames = cellfun(@(varname) varname.getPropName(), equationVarNameList, 'uniformoutput', false);

        end

        
        function funcCallList = getOrderedFunctionCallList(cgt, varargin)
        % Return a list of strings that corresponds to all the function calls ordered in the right order.

            opt = struct('removeExtraVariables', true);
            opt = merge_options(opt, varargin{:});
            
            A                = cgt.adjencyMatrix;
            staticprops      = cgt.staticprops;
            extraVarNameInds = cgt.extraVarNameInds;
            
            funcCallList = {};

            if ~isempty(staticprops)
                
                for istatic = 1 : numel(staticprops)

                    iprop = staticprops{istatic}.propind;
                    propfunction = cgt.model.propertyFunctionList{iprop};
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
                    propfunction = cgt.model.propertyFunctionList{iprop};
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

        function printOrderedFunctionCallList(cgt)
        % Print the function calls ordered in the right order.
            
            funcCallList = cgt.getOrderedFunctionCallList();
            
            fprintf('Function call list\n');
            for ind = 1 : numel(funcCallList)
                fprintf('%s\n', funcCallList{ind});
            end

        end

        function printSubModelNames(cgt, model, parents)

            if nargin < 2
                model = cgt.model;
                parents = {};
            end

            modelnames = model.getSubModelNames();
            for imodelname = 1 : numel(modelnames)
                modelname = modelnames{imodelname};

                subparents = horzcat(parents, {modelname});
                fprintf('%s\n', strjoin(subparents, '.'))
                submodel = model.(modelname);

                cgt.printSubModelNames(submodel, subparents);

            end
            
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
                                    ComputationalGraphTool.getNodeName(varname));
            end
        
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
